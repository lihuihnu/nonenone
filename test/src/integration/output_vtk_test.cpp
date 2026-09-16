/**
 * @file output_vtk_test.cpp
 * @brief 集成测试：验证 `output_vtk_test` 涉及的模块组合、PETSc/网格或输出链路。
 */
#include <output/vtk/cpgrid_vtk_writer.hpp>

#include <cpgrid/cpgrid.hpp>

#include <petscsys.h>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct VtkTag
{
    static constexpr int numVars_ = 1;
};

using Grid =
    MPMC::CpGrid<VtkTag>;

void run(
    const std::string &meshDirectory)
{
    MPMC::Mesh mesh(
        meshDirectory);

    mesh.prepareForUse();

    Grid grid(mesh);

    if (mesh.rank() == 0)
    {
        std::vector<double>
            currentIdField(
                mesh.cellCount(),
                0.0);

        std::vector<double>
            storageField(
                mesh.cellCount(),
                0.0);

        for (const auto &cell :
             mesh.allCells())
        {
            currentIdField[
                static_cast<std::size_t>(
                    cell.id())] =
                static_cast<double>(
                    cell.id());

            storageField[
                static_cast<std::size_t>(
                    cell.storageIndex())] =
                static_cast<double>(
                    cell.storageIndex());
        }

        MPMC::CpGridVtkWriter<Grid>
            writer(grid);

        writer.addCellScalar(
            "current_id_value",
            currentIdField,
            MPMC::VtkCellOrdering::CurrentId);

        writer.addCellScalar(
            "storage_index_value",
            storageField,
            MPMC::VtkCellOrdering::InputIndex);

        writer.write(
            "output_vtk_test.vtu");

        if (!std::filesystem::exists(
                "output_vtk_test.vtu"))
        {
            throw std::runtime_error(
                "VTK writer did not create output_vtk_test.vtu.");
        }

        std::filesystem::remove(
            "output_vtk_test.vtu");
    }

    PetscCallMPIAbort(
        mesh.communicator(),
        MPI_Barrier(
            mesh.communicator()));

    PetscPrintf(
        mesh.communicator(),
        "Output VTK validation: ALL PASS\n");
}

} // namespace

int main(
    int argc,
    char **argv)
{
    PetscInitialize(
        &argc,
        &argv,
        nullptr,
        nullptr);

    int status = 0;

    try
    {
        char meshDirectory[
            PETSC_MAX_PATH_LEN] = {};

        PetscBool found =
            PETSC_FALSE;

        PetscCallAbort(
            PETSC_COMM_WORLD,
            PetscOptionsGetString(
                nullptr,
                nullptr,
                "-mesh_dir",
                meshDirectory,
                sizeof(meshDirectory),
                &found));

        if (!found)
        {
            throw std::invalid_argument(
                "Please provide -mesh_dir <MRST-export-directory>.");
        }

        {
            run(
                std::string(
                    meshDirectory));
        }
    }
    catch (const std::exception &error)
    {
        PetscPrintf(
            PETSC_COMM_WORLD,
            "Output VTK test failed: %s\n",
            error.what());

        status = 1;
    }

    PetscFinalize();
    return status;
}
