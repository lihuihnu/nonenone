/**
 * @file petsc_io.cpp
 * @brief PETSc 向量 CSV 输入输出接口的实现。
 */
#include <common/petsc_io.hpp>

#include <petscviewer.h>

#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace MPMC::petsc
{

namespace
{

/**
 * @brief 将文件名或用户给定名称规范为合法 MATLAB 变量名。
 *
 * MATLAB 标识符必须以字母开头，后续只保留字母、数字和下划线。
 * 文件路径和扩展名不属于变量名，例如
 * `results/phase_state_step_60.m -> phase_state_step_60`。
 */
std::string matlabVariableName(
    const std::string &filename,
    const std::string &requestedName)
{
    std::string name = requestedName;
    if (name.empty())
    {
        const std::size_t slash = filename.find_last_of("/\\");
        const std::size_t begin = slash == std::string::npos ? 0 : slash + 1;
        const std::size_t dot = filename.find_last_of('.');
        const std::size_t end = dot == std::string::npos || dot < begin
            ? filename.size()
            : dot;
        name = filename.substr(begin, end - begin);
    }

    if (name.empty())
        name = "petsc_object";

    for (char &ch : name)
    {
        const unsigned char value = static_cast<unsigned char>(ch);
        if (!std::isalnum(value) && ch != '_')
            ch = '_';
    }

    if (!std::isalpha(static_cast<unsigned char>(name.front())))
        name.insert(0, "v_");

    return name;
}

template <class PetscHandle, class ViewFunction>
PetscErrorCode saveAsciiMatlabObject(
    PetscHandle object,
    const std::string &filename,
    const std::string &variableName,
    const char *nullErrorMessage,
    ViewFunction view)
{
    PetscFunctionBegin;

    PetscCheck(
        object != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "%s",
        nullErrorMessage);

    PetscObject petscObject = reinterpret_cast<PetscObject>(object);
    const MPI_Comm comm = PetscObjectComm(petscObject);
    const std::string matlabName = matlabVariableName(filename, variableName);

    const char *oldNamePointer = nullptr;
    PetscCall(PetscObjectGetName(petscObject, &oldNamePointer));
    const std::string oldName = oldNamePointer != nullptr ? oldNamePointer : "";
    PetscCall(PetscObjectSetName(petscObject, matlabName.c_str()));

    PetscViewer viewer = nullptr;
    PetscCall(PetscViewerASCIIOpen(comm, filename.c_str(), &viewer));
    PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(view(object, viewer));
    PetscCall(PetscViewerPopFormat(viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    // 输出变量名只属于该文件，不改变求解过程中 PETSc 对象的长期名称。
    PetscCall(PetscObjectSetName(petscObject, oldName.c_str()));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

PetscErrorCode saveVectorAsciiMatlab(
    Vec vector,
    const std::string &filename,
    const std::string &variableName)
{
    return saveAsciiMatlabObject(
        vector,
        filename,
        variableName,
        "saveVectorAsciiMatlab requires a valid Vec.",
        [](Vec value, PetscViewer viewer) { return VecView(value, viewer); });
}

PetscErrorCode saveMatrixAsciiMatlab(
    Mat matrix,
    const std::string &filename,
    const std::string &variableName)
{
    return saveAsciiMatlabObject(
        matrix,
        filename,
        variableName,
        "saveMatrixAsciiMatlab requires a valid Mat.",
        [](Mat value, PetscViewer viewer) { return MatView(value, viewer); });
}


namespace
{

std::string csvField(const std::string &value)
{
    bool quote = false;
    for (char ch : value)
    {
        if (ch == ',' || ch == '"' || ch == '\n' || ch == '\r')
        {
            quote = true;
            break;
        }
    }

    if (!quote)
        return value;

    std::string result;
    result.reserve(value.size() + 2);
    result.push_back('"');
    for (char ch : value)
    {
        if (ch == '"')
            result.push_back('"');
        result.push_back(ch);
    }
    result.push_back('"');
    return result;
}

} // namespace

PetscErrorCode saveCellVectorCsv(
    Vec vector,
    PetscInt dofPerCell,
    const std::string &filename,
    const std::vector<std::string> &columnNames)
{
    PetscFunctionBegin;

    PetscCheck(
        vector != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "saveCellVectorCsv requires a valid Vec.");
    PetscCheck(
        dofPerCell > 0,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_OUTOFRANGE,
        "saveCellVectorCsv requires dofPerCell > 0.");
    PetscCheck(
        columnNames.empty() ||
            columnNames.size() == static_cast<std::size_t>(dofPerCell),
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_SIZ,
        "CSV column-name count must match dofPerCell.");

    const MPI_Comm comm =
        PetscObjectComm(reinterpret_cast<PetscObject>(vector));

    PetscInt globalSize = 0;
    PetscCall(VecGetSize(vector, &globalSize));
    PetscCheck(
        globalSize % dofPerCell == 0,
        comm,
        PETSC_ERR_ARG_SIZ,
        "Vector size must be divisible by dofPerCell.");

    PetscMPIInt rank = 0;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    VecScatter scatter = nullptr;
    Vec sequential = nullptr;
    PetscCall(VecScatterCreateToZero(vector, &scatter, &sequential));
    PetscCall(VecScatterBegin(
        scatter, vector, sequential, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(
        scatter, vector, sequential, INSERT_VALUES, SCATTER_FORWARD));

    if (rank == 0)
    {
        const PetscScalar *values = nullptr;
        PetscCall(VecGetArrayRead(sequential, &values));

        std::ofstream file(filename, std::ios::out | std::ios::trunc);
        if (!file.is_open())
        {
            PetscCall(VecRestoreArrayRead(sequential, &values));
            PetscCall(VecScatterDestroy(&scatter));
            PetscCall(VecDestroy(&sequential));
            SETERRQ(
                comm,
                PETSC_ERR_FILE_OPEN,
                "Unable to open CSV output file '%s'.",
                filename.c_str());
        }

        file << "input_index";
        for (PetscInt component = 0; component < dofPerCell; ++component)
        {
            std::string name;
            if (columnNames.empty())
            {
                name = "value_" + std::to_string(component);
            }
            else
            {
                name = columnNames[static_cast<std::size_t>(component)];
            }
            file << ',' << csvField(name);
        }
        file << '\n';

        file << std::setprecision(17) << std::scientific;
        const PetscInt cells = globalSize / dofPerCell;
        for (PetscInt cell = 0; cell < cells; ++cell)
        {
            file << cell;
            for (PetscInt component = 0; component < dofPerCell; ++component)
            {
                const PetscInt index = cell * dofPerCell + component;
                file << ',' << PetscRealPart(values[index]);
            }
            file << '\n';
        }

        file.close();
        PetscCall(VecRestoreArrayRead(sequential, &values));
    }

    PetscCall(VecScatterDestroy(&scatter));
    PetscCall(VecDestroy(&sequential));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode saveVectorBinary(
    Vec vector,
    const std::string &filename)
{
    PetscFunctionBegin;

    PetscCheck(
        vector != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "saveVectorBinary requires a valid Vec.");

    const MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                vector));

    PetscViewer viewer = nullptr;

    PetscCall(
        PetscViewerBinaryOpen(
            comm,
            filename.c_str(),
            FILE_MODE_WRITE,
            &viewer));

    PetscCall(
        VecView(
            vector,
            viewer));

    PetscCall(
        PetscViewerDestroy(
            &viewer));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode loadVectorBinary(
    Vec vector,
    const std::string &filename)
{
    PetscFunctionBegin;

    PetscCheck(
        vector != nullptr,
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_NULL,
        "loadVectorBinary requires a valid Vec.");

    const MPI_Comm comm =
        PetscObjectComm(
            reinterpret_cast<PetscObject>(
                vector));

    PetscViewer viewer = nullptr;

    PetscCall(
        PetscViewerBinaryOpen(
            comm,
            filename.c_str(),
            FILE_MODE_READ,
            &viewer));

    PetscCall(
        VecLoad(
            vector,
            viewer));

    PetscCall(
        PetscViewerDestroy(
            &viewer));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace MPMC::petsc
