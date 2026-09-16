/**
 * @file natural_core_test.cpp
 * @brief 单元测试：验证 `natural_core` 的核心语义、边界条件和回归行为。
 */
#include <indices/indices.hpp>
#include <natural/assembly/cell_residual.hpp>
#include <natural/assembly/local_equations.hpp>
#include <natural/fluid_system.hpp>
#include <natural/physics/accumulation.hpp>
#include <natural/physics/face_flux.hpp>
#include <natural/petsc/natural_layout_registration.hpp>
#include <natural/petsc/natural_scaling.hpp>
#include <natural/physics/land_trapping.hpp>
#include <natural/physics/well_source.hpp>
#include <natural/state/cell_property_evaluator.hpp>
#include <natural/state/newton_limiter.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
using BaseConfig = MPMC::CompositionalModelConfig<6, true, true, false, false, false>;
using Base = MPMC::ScalarIndices<BaseConfig>;
using FullConfig = MPMC::CompositionalModelConfig<6, true, true, true, true, true>;
using Full = MPMC::ScalarIndices<FullConfig>;
using AdBase = MPMC::ADIndices<BaseConfig>;
using AdFull = MPMC::ADIndices<FullConfig>;
using FullyCompositionalConfig = MPMC::CompositionalModelConfig<
    6, true, true, false, false, false,
    MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase>;
using FullyCompositional = MPMC::ScalarIndices<FullyCompositionalConfig>;

void require(bool value, const std::string &message)
{
    if (!value) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
    {
        std::ostringstream details;
        details << message << ": actual=" << std::setprecision(17) << actual
                << ", expected=" << expected
                << ", abs_diff=" << std::abs(actual - expected);
        throw std::runtime_error(details.str());
    }
}

template <class Indices>
MPMC::CompositionalMixture<Indices> makeMixture()
{
    return MPMC::CompositionalMixture<Indices>(
        {189.515, 304.2, 387.607, 597.497, 698.515, 875.0},
        {4580011.59, 7386592.50, 4095515.97, 3345244.875, 1768374.5625, 1169006.79},
        {9.97012032965401e-5, 9.26344713533338e-5, 0.000217076707259486,
         0.000381162235869935, 0.000721410148917871, 0.00113570073874421},
        {0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141},
        {0.0161594, 0.04401, 0.0455725, 0.11774, 0.248827, 0.48152},
        {{0.0, 0.00070981, 0.00077754, 0.0100, 0.0110, 0.0110},
         {0.00070981, 0.0, 0.1500, 0.1500, 0.1500, 0.1500},
         {0.00077754, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0100, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0},
         {0.0110, 0.1500, 0.0, 0.0, 0.0, 0.0}});
}

template <class Indices>
auto makeEos()
{
    return MPMC::CubicEquationOfState<Indices>(
        0.4572355, 0.0779691, makeMixture<Indices>(), 1,
        2.414213562373095, -0.414213562373095);
}

void testEOSAndProperties()
{
    auto eos = makeEos<Base>();
    const std::array<double, 6> x{
        0.3246914, 0.0128351, 0.2278401, 0.2606985, 0.1134144, 0.0605206};
    const std::array<double, 6> y{
        0.808671, 0.025284, 0.1487798, 0.0175878, 0.0006759, 1.51e-06};
    constexpr double p = 15.0e6;
    constexpr double T = 387.45;

    const auto liquid = eos.phaseResult(p, T, x, true);
    const auto vapor = eos.phaseResult(p, T, y, false);
    near(liquid.compressibility, 0.66880987807546899, 2e-12, "PR liquid Z regression");
    near(vapor.compressibility, 0.84105046763526459, 2e-12, "PR vapor Z regression");

    const std::array<double, 6> expectedLiquidF{
        11148938.7289865, 285739.262764987, 958802.358789237, 32132.229549948, 222.962965571357, 0.0303098653818211};
    const std::array<double, 6> expectedVaporF{
        11364784.8024743, 289457.590456904, 964022.546010914, 31595.5016348791, 198.090892606485, 0.0143840907404893};
    for (int i = 0; i < 6; ++i)
    {
        near(liquid.fugacity[static_cast<std::size_t>(i)], expectedLiquidF[static_cast<std::size_t>(i)], 2e-8,
             "PR liquid fugacity regression");
        near(vapor.fugacity[static_cast<std::size_t>(i)], expectedVaporF[static_cast<std::size_t>(i)], 2e-8,
             "PR vapor fugacity regression");
    }

    MPMC::CompositionalPropertyModel<Base> properties(eos.mixture());
    const double rhoL = properties.density(p, x, liquid.compressibility, T);
    const double rhoV = properties.density(p, y, vapor.compressibility, T);
    near(rhoL, 725.8104330709524, 2e-12, "liquid density regression");
    near(rhoV, 128.44402682548275, 2e-12, "vapor density regression");

    const double muL = properties.viscosity(p, x, liquid.compressibility, T);
    const double muV = properties.viscosity(p, y, vapor.compressibility, T);
    near(muL, 1.2351230861591054e-4, 3e-12, "LBC liquid viscosity regression");
    near(muV, 1.8313005681915054e-5, 3e-12, "LBC vapor viscosity regression");

    const auto mass = eos.massFractions(x);
    near(mass[0], 0.05032822, 2e-7, "mass fraction conversion");
    near(mass[5], 0.27953300, 2e-7, "mass fraction conversion heavy component");
}

void testAqueousCO2AdsorptionAndLand()
{
    MPMC::AqueousCO2Model aqueous(387.45, 0.04401, 0.01801528, 0.0);
    near(aqueous.waterSaturationPressure(), 165354.25421818954, 2e-12, "water saturation pressure");
    near(aqueous.referenceHenryConstant(), 551569549.4768968, 2e-12, "reference Henry constant");
    near(aqueous.henryConstant(15.0e6), 654690656.0402844, 2e-12, "pressure corrected Henry constant");
    near(aqueous.massFraction(0.01), 0.02408177995182462, 2e-12, "aqueous CO2 mass fraction");
    near(aqueous.fugacity(15.0e6, 0.01), 6546906.560402844, 2e-12, "aqueous CO2 fugacity");

    MPMC::CompetitiveLangmuirAdsorption<6> adsorption(
        {1.468e-3, 2.860e-3, 0.0, 0.0, 0.0, 0.0},
        {1.471e-7, 2.472e-7, 0.0, 0.0, 0.0, 0.0});
    const std::array<double, 6> y{
        0.808671, 0.025284, 0.1487798, 0.0175878, 0.0006759, 1.51e-06};
    const auto theta = adsorption.adsorbedVolume(15.0e6, y);
    near(theta[0], 0.00091011892411852391, 2e-9, "competitive adsorption component 0");
    near(theta[1], 9.3163936054927677e-05, 2e-9, "competitive adsorption CO2");

    MPMC::LandTrappingModel land(2.0);
    near(land.residualGasAtMaximum(0.6), 0.2727272727272727, 2e-14, "Land Sgr");
    near(land.trappedGasSaturation(0.3, 0.6), 0.16879529857838577, 2e-14, "Land trapped saturation");
    near(land.trappedGasSaturation(0.6, 0.6), 0.0, 1e-14, "Land drainage branch");
}

void testAccumulationAndFlux()
{
    MPMC::CellProperties<Full, double> p{};
    p.porosity = 0.25;
    p.density[Full::Phase::liquid] = 800.0;
    p.density[Full::Phase::vapor] = 100.0;
    p.density[Full::Phase::water] = 1000.0;
    p.saturation[Full::Phase::liquid] = 0.5;
    p.saturation[Full::Phase::vapor] = 0.2;
    p.saturation[Full::Phase::water] = 0.3;
    p.massFraction[Full::Phase::liquid][0] = 0.7;
    p.massFraction[Full::Phase::liquid][1] = 0.3;
    p.massFraction[Full::Phase::vapor][0] = 0.1;
    p.massFraction[Full::Phase::vapor][1] = 0.9;
    p.aqueousCO2MassFraction = 0.02;

    const auto a = MPMC::computeFluidAccumulation<Full>(p, 1);
    near(a.componentMass[0], 0.25 * (800.0 * 0.5 * 0.7 + 100.0 * 0.2 * 0.1), 1e-14, "component accumulation");
    near(a.componentMass[1], 0.25 * (800.0 * 0.5 * 0.3 + 100.0 * 0.2 * 0.9) + 0.25 * 1000.0 * 0.3 * 0.02,
         1e-14, "CO2 accumulation including aqueous mass");
    near(a.waterMass, 0.25 * 1000.0 * 0.3 * 0.98, 1e-14, "H2O accumulation");

    MPMC::CellState<Full, double> left{}, right{};
    left.pressure = 100.0;
    right.pressure = 90.0;
    left.liquidSaturation = right.liquidSaturation = 1.0;
    MPMC::CellProperties<Full, double> lp{}, rp{};
    lp.saturation[Full::Phase::liquid] = rp.saturation[Full::Phase::liquid] = 1.0;
    lp.density[Full::Phase::liquid] = rp.density[Full::Phase::liquid] = 100.0;
    lp.mobility[Full::Phase::liquid] = 0.5;
    rp.mobility[Full::Phase::liquid] = 0.1;
    lp.massFraction[Full::Phase::liquid][0] = 0.7;
    lp.massFraction[Full::Phase::liquid][1] = 0.3;
    rp.massFraction[Full::Phase::liquid][0] = 0.2;
    rp.massFraction[Full::Phase::liquid][1] = 0.8;
    const auto flux = MPMC::computeFaceMassFlux<Full>(
        left, lp, right, rp, 2.0, 0.0, 1.0, 1.0, 0, 1, 1);
    near(flux.component[0], 700.0, 1e-14, "TPFA upwind component 0");
    near(flux.component[1], 300.0, 1e-14, "TPFA upwind component 1");

    MPMC::NaturalScalingOptions scaling;
    scaling.enabled = true;
    scaling.validate();
    bool invalidScalingRejected = false;
    try
    {
        scaling.pressureScale = 0.0;
        scaling.validate();
    }
    catch (const std::invalid_argument &)
    {
        invalidScalingRejected = true;
    }
    require(invalidScalingRejected,
            "Natural consistent scaling must reject non-positive characteristic scales");
}

void testWellAndLimiter()
{
    std::array<double, Base::numPhases> pressure{100.0, 100.0, 100.0};
    std::array<double, Base::numPhases> density{800.0, 100.0, 1000.0};
    std::array<double, Base::numPhases> mobility{0.1, 0.2, 0.3};
    std::array<double, Base::numPhases> surfaceDensity{800.0, 2.0, 1000.0};
    std::array<double, Base::numPhases> phaseFraction{1.0, 0.0, 0.0};
    std::array<double, Base::numComponents> injectionComposition{1.0, 0, 0, 0, 0, 0};
    std::array<std::array<double, Base::numComponents>, Base::numPhases> massFraction{};
    massFraction[Base::Phase::liquid][0] = 0.75;
    massFraction[Base::Phase::liquid][1] = 0.25;
    massFraction[Base::Phase::vapor][0] = 0.1;
    massFraction[Base::Phase::vapor][1] = 0.9;

    const auto well = MPMC::computePerforationWellSource<Base>(
        MPMC::WellType::Producer, 2.0, 90.0, pressure, density, mobility,
        surfaceDensity, phaseFraction, injectionComposition, massFraction);
    near(well.surfacePhaseRate[Base::Phase::liquid], -2.0, 1e-14, "producer oil SC rate");
    near(well.surfacePhaseRate[Base::Phase::vapor], -200.0, 1e-14, "producer gas SC rate");
    near(well.reservoirPhaseRate[Base::Phase::vapor], -4.0, 1e-14, "producer gas RC rate");
    near(well.phaseMassRate[Base::Phase::water], -6000.0, 1e-14, "producer water mass rate");
    near(well.componentMassSource[0], -1240.0, 1e-14, "producer component 0 mass source");
    near(well.componentMassSource[1], -760.0, 1e-14, "producer component 1 mass source");

    // Dissolution branch: produced aqueous phase carries the cell's dissolved CO2
    // mass fraction.  With total water mass rate -6000 and XwCO2=0.02,
    // H2O source is -5880 and dissolved CO2 contributes another -120 to
    // component 1.  This test also forces the constexpr dissolution branch to
    // compile under -Werror.
    std::array<double, Full::numPhases> fullPressure{100.0, 100.0, 100.0};
    std::array<double, Full::numPhases> fullDensity{800.0, 100.0, 1000.0};
    std::array<double, Full::numPhases> fullMobility{0.1, 0.2, 0.3};
    std::array<double, Full::numPhases> fullSurfaceDensity{800.0, 2.0, 1000.0};
    std::array<double, Full::numPhases> fullPhaseFraction{1.0, 0.0, 0.0};
    std::array<double, Full::numComponents> fullInjectionComposition{1.0, 0, 0, 0, 0, 0};
    std::array<std::array<double, Full::numComponents>, Full::numPhases> fullMassFraction{};
    fullMassFraction[Full::Phase::liquid][0] = 0.75;
    fullMassFraction[Full::Phase::liquid][1] = 0.25;
    fullMassFraction[Full::Phase::vapor][0] = 0.1;
    fullMassFraction[Full::Phase::vapor][1] = 0.9;

    const auto dissolvedWell = MPMC::computePerforationWellSource<Full>(
        MPMC::WellType::Producer, 2.0, 90.0,
        fullPressure, fullDensity, fullMobility,
        fullSurfaceDensity, fullPhaseFraction,
        fullInjectionComposition, fullMassFraction,
        0.02, 1);

    near(dissolvedWell.waterMassSource, -5880.0, 1e-14,
         "producer H2O source with dissolved CO2");
    near(dissolvedWell.componentMassSource[0], -1240.0, 1e-14,
         "producer component 0 source with dissolution");
    near(dissolvedWell.componentMassSource[1], -880.0, 1e-14,
         "producer CO2 source including produced aqueous CO2");

    // RATE 控制井方程必须严格使用物理 residual 的真实导数，不能再叠加
    // 与 residual 不对应的固定 BHP 对角项。这里采用与三 EOS 3-D CO2 注入井
    // 同量级的 WI / mobility / pressure drop，对 AD 的 dq/dBHP 做中心差分回归。
    using Eval = AdBase::ValueType;
    std::array<Eval, AdBase::numPhases> adPressure{
        Eval(6.0e6), Eval(6.0e6), Eval(6.0e6)};
    std::array<Eval, AdBase::numPhases> adDensity{
        Eval(560.0), Eval(210.0), Eval(857.0)};
    std::array<Eval, AdBase::numPhases> adMobility{
        Eval(10.0), Eval(4.0e4), Eval(17.0)};
    std::array<double, AdBase::numPhases> adSurfaceDensity{620.0, 1.8, 998.0};
    std::array<double, AdBase::numPhases> adInjectionPhaseFraction{0.0, 1.0, 0.0};
    std::array<double, AdBase::numComponents> adInjectionComposition{};
    adInjectionComposition[1] = 1.0;
    std::array<std::array<Eval, AdBase::numComponents>, AdBase::numPhases>
        adMassFraction{};

    constexpr double testWI = 7.0e-12;
    constexpr double testBhp = 6.1e6;
    const Eval adBhp = Eval::createVariable(
        testBhp, AdBase::Primary::wellPressure);
    const auto adInjector = MPMC::computePerforationWellSource<AdBase>(
        MPMC::WellType::Injector,
        testWI,
        adBhp,
        adPressure,
        adDensity,
        adMobility,
        adSurfaceDensity,
        adInjectionPhaseFraction,
        adInjectionComposition,
        adMassFraction);
    const Eval adRate = MPMC::selectControlledRate<AdBase>(
        MPMC::WellControl::TotalRate, adInjector.surfacePhaseRate);
    const Eval adReservoirRate = MPMC::selectControlledRate<AdBase>(
        MPMC::WellControl::ReservoirTotalRate,
        adInjector.surfacePhaseRate, adInjector.reservoirPhaseRate);
    Eval expectedReservoirRate = 0.0;
    for (const auto &rate : adInjector.reservoirPhaseRate)
        expectedReservoirRate += rate;
    near(MPMC::scalarValue(adReservoirRate),
         MPMC::scalarValue(expectedReservoirRate), 1e-14,
         "reservoir-rate control must select in-situ volume rather than surface volume");

    // G8I runtime 复用 in-place 穿孔 scratch。其值与 AD 导数必须与传统
    // 值返回 API 完全一致，并且同一 scratch 连续覆盖时不能残留上一次结果。
    MPMC::PerforationWellResult<AdBase, Eval> adScratchResult;
    MPMC::PerforationWellWorkspace<AdBase, Eval> adScratchWorkspace;
    MPMC::computePerforationWellSourceInto<AdBase>(
        adScratchResult,
        adScratchWorkspace,
        MPMC::WellType::Injector,
        testWI,
        adBhp,
        adPressure,
        adDensity,
        adMobility,
        adSurfaceDensity,
        adInjectionPhaseFraction,
        adInjectionComposition,
        adMassFraction,
        Eval(0.0),
        -1);
    const Eval adScratchRate = MPMC::selectControlledRate<AdBase>(
        MPMC::WellControl::TotalRate, adScratchResult.surfacePhaseRate);
    near(adScratchRate.value(), adRate.value(), 1e-14,
         "in-place well scratch must preserve RATE value");
    near(adScratchRate.derivative(AdBase::Primary::wellPressure),
         adRate.derivative(AdBase::Primary::wellPressure), 1e-14,
         "in-place well scratch must preserve RATE derivative");

    const Eval shiftedBhp = Eval::createVariable(
        testBhp + 1234.0, AdBase::Primary::wellPressure);
    const auto shiftedReference = MPMC::computePerforationWellSource<AdBase>(
        MPMC::WellType::Injector, testWI, shiftedBhp,
        adPressure, adDensity, adMobility, adSurfaceDensity,
        adInjectionPhaseFraction, adInjectionComposition, adMassFraction);
    MPMC::computePerforationWellSourceInto<AdBase>(
        adScratchResult, adScratchWorkspace,
        MPMC::WellType::Injector, testWI, shiftedBhp,
        adPressure, adDensity, adMobility, adSurfaceDensity,
        adInjectionPhaseFraction, adInjectionComposition, adMassFraction,
        Eval(0.0), -1);
    const Eval shiftedScratchRate = MPMC::selectControlledRate<AdBase>(
        MPMC::WellControl::TotalRate, adScratchResult.surfacePhaseRate);
    const Eval shiftedReferenceRate = MPMC::selectControlledRate<AdBase>(
        MPMC::WellControl::TotalRate, shiftedReference.surfacePhaseRate);
    near(shiftedScratchRate.value(), shiftedReferenceRate.value(), 1e-14,
         "reused well scratch must overwrite prior RATE value");
    near(shiftedScratchRate.derivative(AdBase::Primary::wellPressure),
         shiftedReferenceRate.derivative(AdBase::Primary::wellPressure), 1e-14,
         "reused well scratch must overwrite prior RATE derivative");

    const double adDqDbhp =
        adRate.derivative(AdBase::Primary::wellPressure);

    const auto rateAtBhp = [&](double bhpValue)
    {
        std::array<double, Base::numPhases> p{6.0e6, 6.0e6, 6.0e6};
        std::array<double, Base::numPhases> rho{560.0, 210.0, 857.0};
        std::array<double, Base::numPhases> mob{10.0, 4.0e4, 17.0};
        std::array<double, Base::numPhases> rhoSc{620.0, 1.8, 998.0};
        std::array<double, Base::numPhases> injPhase{0.0, 1.0, 0.0};
        std::array<double, Base::numComponents> injComp{};
        injComp[1] = 1.0;
        std::array<std::array<double, Base::numComponents>, Base::numPhases> mf{};
        const auto injector = MPMC::computePerforationWellSource<Base>(
            MPMC::WellType::Injector,
            testWI,
            bhpValue,
            p,
            rho,
            mob,
            rhoSc,
            injPhase,
            injComp,
            mf);
        return MPMC::selectControlledRate<Base>(
            MPMC::WellControl::TotalRate, injector.surfacePhaseRate);
    };

    constexpr double h = 1.0;
    const double fdDqDbhp =
        (rateAtBhp(testBhp + h) - rateAtBhp(testBhp - h)) / (2.0 * h);
    require(adDqDbhp > 0.0, "injector RATE derivative with respect to BHP must be positive");
    near(adDqDbhp, fdDqDbhp, 1e-10,
         "rate-control AD Jacobian must match finite difference");

    // 每个普通单元的 dummy well-pressure 行默认带 +1 对角；RATE 井代表
    // 单元必须在插入真实 q-target Jacobian 前把这个占位行完全清零。
    MPMC::CellResidualBlock<AdBase, Eval> wellBlock{};
    wellBlock.value[static_cast<std::size_t>(AdBase::Equation::wellControl)] =
        Eval::createVariable(testBhp, AdBase::Primary::wellPressure) - testBhp;
    near(wellBlock.value[static_cast<std::size_t>(AdBase::Equation::wellControl)]
             .derivative(AdBase::Primary::wellPressure),
         1.0, 1e-14, "default dummy well-control row must carry identity derivative");
    MPMC::clearRateWellControlPlaceholder<AdBase>(wellBlock);
    near(wellBlock.value[static_cast<std::size_t>(AdBase::Equation::wellControl)]
             .derivative(AdBase::Primary::wellPressure),
         0.0, 1e-14, "RATE representative must remove dummy well-control identity row");

    std::array<double, Full::numPrimaryVariables> delta{};
    std::array<double, Full::numPrimaryVariables> state{};
    state[Full::Primary::pressure] = 100.0;
    delta[Full::Primary::pressure] = 50.0;
    delta[Full::Primary::liquidSaturation] = 0.2;
    delta[Full::Primary::vaporSaturation] = -0.2;
    delta[Full::Primary::aqueousCO2MoleFraction] = 0.04;
    MPMC::limitNaturalNewtonIncrement<Full>(delta, state, MPMC::HydrocarbonPhaseState::TwoPhase);
    near(delta[Full::Primary::pressure], 25.0, 1e-14, "pressure Newton limiter");
    near(delta[Full::Primary::liquidSaturation], 0.1, 1e-14, "saturation Newton limiter");
    near(delta[Full::Primary::aqueousCO2MoleFraction], 0.02, 1e-14, "aqueous Newton limiter");

    // Fully-compositional regression: a pathological composition correction in
    // one phase must not freeze saturation or the other two composition blocks.
    // The pre-v59 limiter used one global thermodynamic scale, so the O(100)
    // oil-composition correction below would shrink every thermodynamic update
    // by O(1e-3), reproducing the near-phase-boundary residual plateau seen in
    // the three-EOS HPC logs.
    std::array<double, FullyCompositional::numPrimaryVariables> threePhaseDelta{};
    std::array<double, FullyCompositional::numPrimaryVariables> threePhaseState{};
    threePhaseState[FullyCompositional::Primary::pressure] = 100.0;
    threePhaseDelta[FullyCompositional::Primary::liquidSaturation] = 0.05;
    threePhaseDelta[FullyCompositional::Primary::vaporSaturation] = -0.04;
    threePhaseDelta[FullyCompositional::Primary::waterSaturation] = -0.01;
    threePhaseDelta[static_cast<std::size_t>(
        FullyCompositional::Primary::liquidComposition[0])] = 100.0;
    threePhaseDelta[static_cast<std::size_t>(
        FullyCompositional::Primary::vaporComposition[0])] = 0.02;
    threePhaseDelta[static_cast<std::size_t>(
        FullyCompositional::Primary::waterComposition[0])] = -0.03;

    MPMC::limitNaturalNewtonIncrement<FullyCompositional>(
        threePhaseDelta, threePhaseState, MPMC::PhasePresence::all());

    near(threePhaseDelta[FullyCompositional::Primary::liquidSaturation], 0.05, 1e-14,
         "trace-phase composition must not freeze oil saturation update");
    near(threePhaseDelta[FullyCompositional::Primary::vaporSaturation], -0.04, 1e-14,
         "trace-phase composition must not freeze gas saturation update");
    near(threePhaseDelta[FullyCompositional::Primary::waterSaturation], -0.01, 1e-14,
         "trace-phase composition must not freeze water saturation update");
    near(threePhaseDelta[static_cast<std::size_t>(
             FullyCompositional::Primary::liquidComposition[0])],
         0.1, 1e-14, "oil composition block must still be safely limited");
    near(threePhaseDelta[static_cast<std::size_t>(
             FullyCompositional::Primary::vaporComposition[0])],
         0.02, 1e-14, "gas composition block should keep its own Newton scale");
    near(threePhaseDelta[static_cast<std::size_t>(
             FullyCompositional::Primary::waterComposition[0])],
         -0.03, 1e-14, "water composition block should keep its own Newton scale");

    // Saturation closure direction is still limited as one block.
    threePhaseDelta.fill(0.0);
    threePhaseDelta[FullyCompositional::Primary::liquidSaturation] = 0.2;
    threePhaseDelta[FullyCompositional::Primary::vaporSaturation] = -0.12;
    threePhaseDelta[FullyCompositional::Primary::waterSaturation] = -0.08;
    MPMC::limitNaturalNewtonIncrement<FullyCompositional>(
        threePhaseDelta, threePhaseState, MPMC::PhasePresence::all());
    near(threePhaseDelta[FullyCompositional::Primary::liquidSaturation], 0.1, 1e-14,
         "three-phase saturation limiter maximum change");
    near(threePhaseDelta[FullyCompositional::Primary::vaporSaturation], -0.06, 1e-14,
         "three-phase saturation limiter preserves relative direction");
    near(threePhaseDelta[FullyCompositional::Primary::waterSaturation], -0.04, 1e-14,
         "three-phase saturation limiter preserves closure direction");

    // The H2O/CO2 aqueous closure guard must include the dependent last
    // component and scale only the water-composition Newton block.
    threePhaseState.fill(0.0);
    threePhaseDelta.fill(0.0);
    threePhaseState[static_cast<std::size_t>(
        FullyCompositional::Primary::waterComposition[0])] = 0.9;
    threePhaseState[static_cast<std::size_t>(
        FullyCompositional::Primary::waterComposition[1])] = 0.09995;
    threePhaseDelta[static_cast<std::size_t>(
        FullyCompositional::Primary::waterComposition[0])] = -0.01;
    threePhaseDelta[FullyCompositional::Primary::liquidSaturation] = 0.03;
    MPMC::limitAqueousCompositionNewtonIncrement<FullyCompositional>(
        threePhaseDelta, threePhaseState, 0, 1, 1.0e-4);
    near(threePhaseDelta[static_cast<std::size_t>(
             FullyCompositional::Primary::waterComposition[0])],
         -4.999999e-5, 5e-8,
         "aqueous domain limiter must bound the dependent unsupported component");
    near(threePhaseDelta[FullyCompositional::Primary::liquidSaturation], 0.03, 1e-14,
         "aqueous domain limiter must not alter saturation updates");
}


void testNaturalLayoutRegistration()
{
    struct RecordingBackend final
    {
        std::array<int, 3> calls{};
        int count{0};

        void registerLayout(int dof)
        {
            calls[static_cast<std::size_t>(count++)] = dof;
        }

        void registerAdAuxiliaryLayout(int dof)
        {
            calls[static_cast<std::size_t>(count++)] = -dof;
        }
    };

    RecordingBackend scalar;
    MPMC::registerNaturalGridLayouts<Base>(scalar);
    require(scalar.count == 2, "scalar Natural layout registration count");
    require(scalar.calls[0] == Base::numPrimaryVariables,
            "scalar Natural primary layout registration");
    require(scalar.calls[1] == Base::numPhaseStateVariables,
            "scalar Natural phase-state layout registration");

    RecordingBackend ad;
    MPMC::registerNaturalGridLayouts<AdBase>(ad);
    require(ad.count == 3, "AD Natural layout registration count");
    require(ad.calls[0] == AdBase::numPrimaryVariables,
            "AD Natural primary layout registration");
    require(ad.calls[1] == -(1 + AdBase::numPrimaryVariables),
            "AD Natural auxiliary layout registration");
    require(ad.calls[2] == AdBase::numPhaseStateVariables,
            "AD Natural phase-state layout registration");

    struct LegacyBackend final
    {
        std::array<int, 2> calls{};
        int count{0};

        void registerLayout(int dof)
        {
            calls[static_cast<std::size_t>(count++)] = dof;
        }
    };

    LegacyBackend legacy;
    MPMC::registerNaturalGridLayouts<AdBase>(legacy);
    require(legacy.count == 2,
            "legacy backend must not require AD auxiliary hook");
    require(legacy.calls[0] == AdBase::numPrimaryVariables &&
                legacy.calls[1] == AdBase::numPhaseStateVariables,
            "legacy backend Natural layout registration compatibility");
}

void testAD()
{
    auto eos = makeEos<AdFull>();
    using Eval = AdFull::ValueType;
    std::array<Eval, 6> x{};
    const std::array<double, 6> xd{
        0.3246914, 0.0128351, 0.2278401, 0.2606985, 0.1134144, 0.0605206};
    for (int i = 0; i < 6; ++i) x[static_cast<std::size_t>(i)] = Eval(xd[static_cast<std::size_t>(i)]);
    Eval p = Eval::createVariable(15.0e6, AdFull::Primary::pressure);
    const auto phase = eos.phaseResult(p, Eval(387.45), x, true);
    require(std::isfinite(phase.compressibility.value()), "AD EOS value must be finite");
    require(std::isfinite(phase.compressibility.derivative(AdFull::Primary::pressure)), "AD EOS pressure derivative must be finite");

    MPMC::AqueousCO2Model aqueous(387.45, 0.04401, 0.01801528);
    Eval xw = Eval::createVariable(0.01, AdFull::Primary::aqueousCO2MoleFraction);
    const Eval mass = aqueous.massFraction(xw);
    require(mass.derivative(AdFull::Primary::aqueousCO2MoleFraction) > 0.0,
            "aqueous mass fraction derivative must be positive");

    MPMC::LandTrappingModel land(2.0);
    Eval sg = Eval::createVariable(0.3, AdFull::Primary::vaporSaturation);
    const Eval trapped = land.trappedGasSaturation(sg, 0.6);
    require(std::isfinite(trapped.derivative(AdFull::Primary::vaporSaturation)),
            "Land AD derivative must be finite");
}

} // namespace

int main()
{
    testEOSAndProperties();
    testAqueousCO2AdsorptionAndLand();
    testAccumulationAndFlux();
    testWellAndLimiter();
    testNaturalLayoutRegistration();
    testAD();
    std::cout << "Natural core regression: ALL PASS\n";
    return 0;
}
