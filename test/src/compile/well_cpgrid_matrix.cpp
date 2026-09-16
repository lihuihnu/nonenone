/**
 * @file well_cpgrid_matrix.cpp
 * @brief 编译矩阵测试：验证 `well_cpgrid_matrix` 所覆盖模板/后端组合能够实例化。
 */
#include <indices/indices.hpp>
#include <natural/petsc/well_runtime.hpp>
#include <well/well.hpp>

#include <type_traits>

using Config = MPMC::CompositionalModelConfig<6, true, true>;
using Indices = MPMC::ADIndices<Config>;

static_assert(std::is_same_v<
    MPMC::NaturalPerforation<Indices>,
    MPMC::WellPerforation<PetscInt>>);

static_assert(std::is_same_v<
    MPMC::NaturalWell<Indices>,
    MPMC::WellSpecification<Indices, PetscInt>>);

int main()
{
    MPMC::NaturalWell<Indices> well;
    well.type = MPMC::WellType::Producer;
    well.control = MPMC::WellControl::Bhp;
    well.target = 10.0e6;
    well.bhpCellId = 0;
    well.perforations = {{0, 1.0e-12}};
    well.validate(1);
    return 0;
}
