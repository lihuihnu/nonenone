/**
 * @file output_cpgrid_matrix.cpp
 * @brief 编译矩阵测试：验证 `output_cpgrid_matrix` 所覆盖模板/后端组合能够实例化。
 */
#include <output/petsc/output_petsc.hpp>

#include <cpgrid/cpgrid.hpp>
#include <indices/indices.hpp>

namespace
{

struct DummyTag
{
    static constexpr int numVars_ = 15;
};

using Grid =
    MPMC::CpGrid<DummyTag>;

using Base =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, false, false>>;

using Dissolution =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            true, false, false>>;

using Land =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, false, true>>;

using Adsorption =
    MPMC::ScalarIndices<
        MPMC::CompositionalModelConfig<
            6, true, true,
            false, true, false>>;

} // namespace

/*
 * Explicit instantiations must appear in a namespace that encloses the
 * namespace where CpGridSaver is declared.  MPMC is declared directly
 * in the global namespace, so these instantiations belong here, not inside
 * the anonymous namespace above.
 */
template class MPMC::CpGridSaver<
    Base,
    Grid>;

template class MPMC::CpGridSaver<
    Dissolution,
    Grid>;

template class MPMC::CpGridSaver<
    Land,
    Grid>;

template class MPMC::CpGridSaver<
    Adsorption,
    Grid>;
