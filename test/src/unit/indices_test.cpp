/**
 * @file indices_test.cpp
 * @brief 单元测试：验证 `indices` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>

#include <ad/Evaluation.hpp>

#include <array>
#include <iostream>
#include <type_traits>

namespace
{

using ThreePhaseSixComponent =
    MPMC::CompositionalModelConfig<
        6,      // components
        true,   // water
        true,   // well unknown
        false,  // aqueous CO2 dissolution
        false,  // adsorption
        false>; // Land trapping

using ThreePhaseSixComponentIndices =
    MPMC::ADIndices<
        ThreePhaseSixComponent>;

using ThreePhaseSixComponentDissolved =
    MPMC::CompositionalModelConfig<
        6,
        true,
        true,
        true,
        false,
        false>;

using ThreePhaseSixComponentDissolvedIndices =
    MPMC::ADIndices<
        ThreePhaseSixComponentDissolved>;

using TwoPhaseSixComponent =
    MPMC::CompositionalModelConfig<
        6,
        false,
        true>;

using TwoPhaseSixComponentIndices =
    MPMC::ADIndices<
        TwoPhaseSixComponent>;

using SectorTwentyComponent =
    MPMC::CompositionalModelConfig<
        20,
        true,
        true>;

using SectorTwentyComponentIndices =
    MPMC::ADIndices<
        SectorTwentyComponent>;

using DaqingDissolved =
    MPMC::CompositionalModelConfig<
        8,
        true,
        true,
        true,
        false,
        true>;

using DaqingDissolvedIndices =
    MPMC::ADIndices<
        DaqingDissolved>;

using AdsorptionAndTrapping =
    MPMC::CompositionalModelConfig<
        6,
        true,
        true,
        false,
        true,
        true>;

using AdsorptionAndTrappingIndices =
    MPMC::ADIndices<
        AdsorptionAndTrapping>;

using ScalarConfig =
    MPMC::ScalarIndices<
        ThreePhaseSixComponent>;

using FullyCompositionalFourComponent =
    MPMC::CompositionalModelConfig<
        4, true, false, false, false, false,
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using FullyCompositionalFourComponentIndices =
    MPMC::ADIndices<FullyCompositionalFourComponent>;

using FullyCompositionalFourComponentWell =
    MPMC::CompositionalModelConfig<
        4, true, true, false, false, false,
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using FullyCompositionalFourComponentWellIndices =
    MPMC::ADIndices<FullyCompositionalFourComponentWell>;

template <class LhsArray, class RhsArray>
constexpr bool arraysEqual(
    const LhsArray &lhs,
    const RhsArray &rhs)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }

    for (std::size_t i = 0;
         i < lhs.size();
         ++i)
    {
        if (lhs[i] != rhs[i])
        {
            return false;
        }
    }

    return true;
}

// -------------------------------------------------------------------------
// 3p6c base layout: must remain numerically identical to the original code.
// -------------------------------------------------------------------------
static_assert(
    ThreePhaseSixComponentIndices::
        numPrimaryVariables == 15);

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::pressure == 0);

static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            Primary::liquidComposition,
        std::array<int, 5>{
            1, 2, 3, 4, 5}));

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::waterSaturation == 6);

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::wellPressure == 7);

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::liquidSaturation == 8);

static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            Primary::vaporComposition,
        std::array<int, 5>{
            9, 10, 11, 12, 13}));

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::vaporSaturation == 14);

static_assert(
    ThreePhaseSixComponentIndices::
        Primary::aqueousCO2MoleFraction ==
    ThreePhaseSixComponentIndices::
        disabledIndex);

// Equations
static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            Equation::massConservation,
        std::array<int, 6>{
            0, 1, 2, 3, 4, 5}));

static_assert(
    ThreePhaseSixComponentIndices::
        Equation::waterConservation == 6);

static_assert(
    ThreePhaseSixComponentIndices::
        Equation::wellControl == 7);

static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            Equation::fugacity,
        std::array<int, 6>{
            8, 9, 10, 11, 12, 13}));

static_assert(
    ThreePhaseSixComponentIndices::
        Equation::volumeClosure == 14);

// Phase state
static_assert(
    ThreePhaseSixComponentIndices::
        PhaseState::flag == 0);

static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            PhaseState::equilibriumRatio,
        std::array<int, 6>{
            1, 2, 3, 4, 5, 6}));

static_assert(
    arraysEqual(
        ThreePhaseSixComponentIndices::
            PhaseState::overallComposition,
        std::array<int, 6>{
            7, 8, 9, 10, 11, 12}));

static_assert(
    ThreePhaseSixComponentIndices::
        PhaseState::liquidFraction == 13);

static_assert(
    ThreePhaseSixComponentIndices::
        PhaseState::liquidCompressibilityFactor == 14);

static_assert(
    ThreePhaseSixComponentIndices::
        PhaseState::vaporCompressibilityFactor == 15);

static_assert(
    ThreePhaseSixComponentIndices::
        numPhaseStateVariables == 16);

// -------------------------------------------------------------------------
// Aqueous CO2 dissolution adds exactly one primary and one equation.
// -------------------------------------------------------------------------
static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        Primary::baseCount == 15);

static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        Primary::aqueousCO2MoleFraction == 15);

static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        Equation::aqueousCO2Equilibrium == 15);

static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        numPrimaryVariables == 16);

static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        numEquations == 16);

// -------------------------------------------------------------------------
// No-water case.
// -------------------------------------------------------------------------
static_assert(
    TwoPhaseSixComponentIndices::
        numPhases == 2);

static_assert(
    TwoPhaseSixComponentIndices::
        Phase::water ==
    TwoPhaseSixComponentIndices::
        disabledIndex);

static_assert(
    TwoPhaseSixComponentIndices::
        Primary::waterSaturation ==
    TwoPhaseSixComponentIndices::
        disabledIndex);

static_assert(
    TwoPhaseSixComponentIndices::
        Equation::waterConservation ==
    TwoPhaseSixComponentIndices::
        disabledIndex);

static_assert(
    TwoPhaseSixComponentIndices::
        Primary::wellPressure == 6);

static_assert(
    TwoPhaseSixComponentIndices::
        Primary::liquidSaturation == 7);

static_assert(
    TwoPhaseSixComponentIndices::
        Primary::vaporSaturation == 13);

static_assert(
    TwoPhaseSixComponentIndices::
        numPrimaryVariables == 14);

// -------------------------------------------------------------------------
// 20-component 3-phase model naturally produces 43 AD derivatives.
// This replaces the historical need for Evaluation43.hpp.
// -------------------------------------------------------------------------
static_assert(
    SectorTwentyComponentIndices::
        numPrimaryVariables == 43);

static_assert(
    SectorTwentyComponentIndices::
        ValueType::size() == 43);

// -------------------------------------------------------------------------
// 8-component water/well/dissolution model naturally produces 20 vars.
// -------------------------------------------------------------------------
static_assert(
    DaqingDissolvedIndices::
        numPrimaryVariables == 20);

static_assert(
    DaqingDissolvedIndices::
        Primary::aqueousCO2MoleFraction == 19);

static_assert(
    DaqingDissolvedIndices::
        hasLandTrapping);

// -------------------------------------------------------------------------
// Adsorption/Land are physics capabilities, not extra algebraic unknowns.
// -------------------------------------------------------------------------
static_assert(
    AdsorptionAndTrappingIndices::
        hasAdsorption);

static_assert(
    AdsorptionAndTrappingIndices::
        hasLandTrapping);

static_assert(
    AdsorptionAndTrappingIndices::
        numPrimaryVariables == 15);

// -------------------------------------------------------------------------
// Fully compositional O/G/W layout: H2O is one of N common components.
// N=4 gives p + 3*(N-1) compositions + 3 saturations = 13 primaries.
// The phase-state block is mask + Kg/o[N] + Kw/o[N] + z[N] + beta[3] + Z[3] + hysteresis-mask.
// -------------------------------------------------------------------------
static_assert(FullyCompositionalFourComponentIndices::fullyCompositionalThreePhase);
static_assert(FullyCompositionalFourComponentIndices::numThermodynamicPhases == 3);
static_assert(!FullyCompositionalFourComponentIndices::hasIndependentWaterConservation);
static_assert(FullyCompositionalFourComponentIndices::numPrimaryVariables == 13);
static_assert(FullyCompositionalFourComponentIndices::numEquations == 13);
static_assert(FullyCompositionalFourComponentIndices::numPhaseStateVariables == 20);
static_assert(FullyCompositionalFourComponentIndices::Primary::liquidCompositionBegin == 1);
static_assert(FullyCompositionalFourComponentIndices::Primary::vaporCompositionBegin == 4);
static_assert(FullyCompositionalFourComponentIndices::Primary::aqueousCompositionBegin == 7);
static_assert(FullyCompositionalFourComponentIndices::Primary::liquidSaturation == 10);
static_assert(FullyCompositionalFourComponentIndices::Primary::vaporSaturation == 11);
static_assert(FullyCompositionalFourComponentIndices::Primary::waterSaturation == 12);
static_assert(FullyCompositionalFourComponentIndices::Equation::fugacityBegin == 4);
static_assert(FullyCompositionalFourComponentIndices::Equation::aqueousFugacityBegin == 8);
static_assert(FullyCompositionalFourComponentIndices::Equation::volumeClosure == 12);
static_assert(FullyCompositionalFourComponentWellIndices::numPrimaryVariables == 14);
static_assert(FullyCompositionalFourComponentWellIndices::Equation::wellControl == 13);

// -------------------------------------------------------------------------
// Scalar/AD numerical policy is separate from physical config.
// -------------------------------------------------------------------------
static_assert(
    std::is_same_v<
        ScalarConfig::ValueType,
        double>);

static_assert(
    !ScalarConfig::
        useAutomaticDifferentiation);

static_assert(
    ThreePhaseSixComponentIndices::
        useAutomaticDifferentiation);

static_assert(
    ThreePhaseSixComponentIndices::
        ValueType::size() ==
    ThreePhaseSixComponentIndices::
        numPrimaryVariables);

// System is square for every tested formulation.
static_assert(
    ThreePhaseSixComponentIndices::
        numPrimaryVariables ==
    ThreePhaseSixComponentIndices::
        numEquations);

static_assert(
    ThreePhaseSixComponentDissolvedIndices::
        numPrimaryVariables ==
    ThreePhaseSixComponentDissolvedIndices::
        numEquations);

static_assert(
    TwoPhaseSixComponentIndices::
        numPrimaryVariables ==
    TwoPhaseSixComponentIndices::
        numEquations);

} // namespace

int main()
{
    std::cout
        << "Indices validation: ALL PASS\n";

    std::cout
        << "3p6c vars = "
        << ThreePhaseSixComponentIndices::
               numPrimaryVariables
        << '\n';

    std::cout
        << "3p6c+dissolution vars = "
        << ThreePhaseSixComponentDissolvedIndices::
               numPrimaryVariables
        << '\n';

    std::cout
        << "2p6c vars = "
        << TwoPhaseSixComponentIndices::
               numPrimaryVariables
        << '\n';

    std::cout
        << "20-component vars = "
        << SectorTwentyComponentIndices::
               numPrimaryVariables
        << '\n';

    std::cout
        << "8-component+dissolution vars = "
        << DaqingDissolvedIndices::
               numPrimaryVariables
        << '\n';

    return 0;
}
