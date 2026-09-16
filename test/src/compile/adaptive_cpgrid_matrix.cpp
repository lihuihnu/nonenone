/**
 * @file adaptive_cpgrid_matrix.cpp
 * @brief 编译矩阵测试：验证 `adaptive_cpgrid_matrix` 所覆盖模板/后端组合能够实例化。
 */
#include <adaptive_timestep/core/stepper.hpp>
#include <adaptive_timestep/natural/natural_petsc_backend.hpp>

#include <natural/petsc/natural_cpgrid_runtime.hpp>
#include <natural/numerics.hpp>
#include <indices/indices.hpp>

namespace
{

template <class Indices>
using Grid = MPMC::CpGridCore;

template <class Indices>
using Runtime = MPMC::NaturalCpGridRuntime<
    Indices,
    Grid<Indices>>;

using BaseIndices = MPMC::ADIndices<
    MPMC::CompositionalModelConfig<6, true, true, false, false, false>>;

using DissolutionIndices = MPMC::ADIndices<
    MPMC::CompositionalModelConfig<6, true, true, true, false, false>>;

using LandIndices = MPMC::ADIndices<
    MPMC::CompositionalModelConfig<6, true, true, false, false, true>>;

using AdsorptionIndices = MPMC::ADIndices<
    MPMC::CompositionalModelConfig<6, true, true, false, true, false>>;

using BaseBackend = MPMC::NaturalAdaptiveBackend<Runtime<BaseIndices>>;
using DissolutionBackend = MPMC::NaturalAdaptiveBackend<Runtime<DissolutionIndices>>;
using LandBackend = MPMC::NaturalAdaptiveBackend<Runtime<LandIndices>>;
using AdsorptionBackend = MPMC::NaturalAdaptiveBackend<Runtime<AdsorptionIndices>>;

static_assert(sizeof(BaseBackend) > 0);
static_assert(sizeof(DissolutionBackend) > 0);
static_assert(sizeof(LandBackend) > 0);
static_assert(sizeof(AdsorptionBackend) > 0);

} // namespace

// Explicit instantiations belong in global namespace because MPMC is declared
// directly in global namespace.
template class MPMC::AdaptiveTimeStepper<BaseBackend>;
template class MPMC::AdaptiveTimeStepper<DissolutionBackend>;
template class MPMC::AdaptiveTimeStepper<LandBackend>;
template class MPMC::AdaptiveTimeStepper<AdsorptionBackend>;
