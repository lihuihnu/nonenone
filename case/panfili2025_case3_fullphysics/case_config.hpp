#pragma once

#include <indices/model_config.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <common/units.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>

/**
 * @file case_config.hpp
 * @brief Panfili 2025 Case-3 全物理算例的网格、流体、时间步和求解参数。
 *
 * Literature target:
 *   - depleted gas-condensate CO2-storage problem;
 *   - 30% connate-water variant;
 *   - pure-CO2 injection;
 *   - full live-water physics (dissolution + water vaporization);
 *   - Soreide-Whitson EOS and all-component O/G/W equilibrium.
 *
 * The paper reports the Hamilton CCS grid dimensions and average rock data but
 * does not publish the complete Eclipse deck, exact completion cell indices,
 * hydrocarbon BIP table, case-#3 brine salinity, or raw relative-permeability
 * tables.  Those gaps are kept explicit below.  This production case is a
 * runnable, paper-derived average-property reproduction, not a claim that the
 * source/archived field model has been reproduced cell-for-cell.
 */
namespace CaseConfig
{

inline constexpr char name[] = "panfili2025_case3_fullphysics";
inline constexpr double bar = 1.0e5;
inline constexpr double psi = 6894.757293168;
inline constexpr double foot = 0.3048;
inline constexpr double mD = 9.869232667160130e-16;
inline constexpr double gasConstant = MPMC::units::gasConstant;
inline constexpr double secondsPerDay = 86400.0;
inline constexpr double daysPerYear = 365.25;
inline constexpr double standardCubicFoot = 0.028316846592;

// ============================================================================
// 1. Model
// ============================================================================
struct Model
{
    // H2O + the nine non-water components reported in Panfili Table 3.
    static constexpr int numberOfComponents = 10;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false; // H2O/CO2 are ordinary components.
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

// ============================================================================
// 2. Paper-reported Hamilton geometry/average rock properties.
// ============================================================================
struct Grid
{
    inline static constexpr const char *meshDirectory = "";

    // Panfili et al. state 47 x 99 x 15 cells.  Note that this product is
    // 69,795, although the same paragraph also says "approximately 400k active
    // cells".  Those two statements are internally inconsistent; the explicit
    // grid dimensions are retained here because they are directly reproducible.
    static constexpr int nx = 47;
    static constexpr int ny = 99;
    static constexpr int nz = 15;

    static constexpr double dx = 328.0 * foot;
    static constexpr double dy = 328.0 * foot;
    static constexpr double dz = 8.0 * foot;
    static constexpr double lx = nx * dx;
    static constexpr double ly = ny * dy;
    static constexpr double lz = nz * dz;
};

struct Rock
{
    // Paper: average horizontal permeability 815 mD and porosity 15%.
    static constexpr double kx = 815.0 * mD;
    static constexpr double ky = 815.0 * mD;

    // Vertical permeability is not reported in the paper text.  The runnable
    // structured surrogate therefore uses isotropic permeability rather than
    // inventing an anisotropy ratio.  Replace this with Hamilton deck values
    // when the official Eclipse data are imported.
    static constexpr double kz = 815.0 * mD;
    static constexpr double porosity = 0.15;
};

// ============================================================================
// 3. Soreide-Whitson fluid: Panfili Table 3 + documented closures.
// ============================================================================
struct Fluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::SoreideWhitson;

    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int co2Component = 1;
    static constexpr int nitrogenComponent = 2;
    static constexpr int methaneComponent = 3;

    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "CO2", "N2", "C1", "C2", "C3", "C4-6", "C7+1", "C7+2", "C7+3"};

    // Panfili Table 3 gives the nine non-water components in degR / psia.
    // H2O is the standard SW water component.
    inline static constexpr std::array<double, N> criticalTemperature{
        647.30,
        548.46 * 5.0 / 9.0,
        227.16 * 5.0 / 9.0,
        343.08 * 5.0 / 9.0,
        549.77 * 5.0 / 9.0,
        665.64 * 5.0 / 9.0,
        806.54 * 5.0 / 9.0,
        838.11 * 5.0 / 9.0,
        1058.04 * 5.0 / 9.0,
        1291.89 * 5.0 / 9.0};

    inline static constexpr std::array<double, N> criticalPressure{
        220.48e5,
        1071.33 * psi,
        492.31 * psi,
        667.78 * psi,
        708.34 * psi,
        618.70 * psi,
        514.93 * psi,
        410.75 * psi,
        247.56 * psi,
        160.42 * psi};

    inline static constexpr std::array<double, N> acentricFactor{
        0.344, 0.22500, 0.04000, 0.01300, 0.09860,
        0.15240, 0.21575, 0.31230, 0.55670, 0.91692};

    inline static constexpr std::array<double, N> molarMass{
        0.01801528,
        0.04401,
        0.028013,
        0.016043,
        0.03007,
        0.044097,
        0.066869,
        0.107779,
        0.198562,
        0.335198};

    // Table 3 does not report critical volumes.  Current Natural LBC-style
    // viscosity needs Vc although PR/SW fugacity does not.  Use a transparent
    // Zc=0.27 estimate for pseudo-components; keep the standard water Vc.
    inline static constexpr std::array<double, N> criticalVolume{
        5.6e-5,
        0.27 * gasConstant * criticalTemperature[1] / criticalPressure[1],
        0.27 * gasConstant * criticalTemperature[2] / criticalPressure[2],
        0.27 * gasConstant * criticalTemperature[3] / criticalPressure[3],
        0.27 * gasConstant * criticalTemperature[4] / criticalPressure[4],
        0.27 * gasConstant * criticalTemperature[5] / criticalPressure[5],
        0.27 * gasConstant * criticalTemperature[6] / criticalPressure[6],
        0.27 * gasConstant * criticalTemperature[7] / criticalPressure[7],
        0.27 * gasConstant * criticalTemperature[8] / criticalPressure[8],
        0.27 * gasConstant * criticalTemperature[9] / criticalPressure[9]};

    // The paper says the case-#3 EOS is "inspired from" SPE comparative
    // solution project #3, but does not reproduce its BIP table.  We therefore
    // use the published SPE3 nine-component interaction pattern for the
    // non-water block as an explicit MPMC closure, while SW-recommended
    // non-aqueous H2O interactions are used in row/column 0.
    //
    // IMPORTANT: this matrix is not claimed to be the undisclosed exact BIP
    // matrix used by Panfili et al.
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction = [] {
        std::array<std::array<double, N>, N> k{};
        auto set = [&](int i, int j, double value) {
            k[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
            k[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = value;
        };

        // SW-recommended non-aqueous water interactions.
        set(0, 1, 0.1896); // H2O-CO2
        set(0, 2, 0.5000); // H2O-N2
        set(0, 3, 0.4850); // H2O-C1
        for (int c = 4; c < static_cast<int>(N); ++c)
            set(0, c, 0.5000);

        // SPE3-inspired non-water block, same component grouping as Table 3.
        set(1, 2, -0.0200);
        set(1, 3,  0.1000); set(2, 3, 0.0360);
        set(1, 4,  0.1300); set(2, 4, 0.0500);
        set(1, 5,  0.1350); set(2, 5, 0.0800);
        set(1, 6,  0.1277); set(2, 6, 0.1002); set(3, 6, 0.09281);
        set(1, 7,  0.1000); set(2, 7, 0.1000); set(4, 7, 0.00385); set(5, 7, 0.00385);
        set(1, 8,  0.1000); set(2, 8, 0.1000); set(4, 8, 0.00630); set(5, 8, 0.00630);
        set(1, 9,  0.1000); set(2, 9, 0.1000); set(3, 9, 0.13920);
        set(4, 9,  0.00600); set(5, 9, 0.00600);
        return k;
    }();

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    // Panfili explicitly transitions to the original Peng-Robinson EOS on
    // which SW is based, rather than a corrected-PR legacy variant.
    static constexpr int eosModelFlag = 1;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;

    static constexpr double temperature = (200.0 - 32.0) * 5.0 / 9.0 + 273.15;

    // Section 7 does not state case-#3 salinity.  The runnable wrapper therefore
    // explicitly assumes fresh water (zero molality) until the source value is recovered.  Keep this as one explicit parameter so an exact
    // Hamilton brine value can be inserted when recovered from the source deck.
    static constexpr double soreideWhitsonSalinityMolality = 0.0;

    static double soreideWhitsonAqueousWaterBip(
        int component,
        double temperatureK,
        double salinityMolality)
    {
        if (component == waterComponent)
            return 0.0;
        if (component == co2Component)
            return MPMC::SoreideWhitsonCorrelations::co2AqueousBip(
                temperatureK, criticalTemperature[co2Component], salinityMolality);
        if (component == nitrogenComponent)
            return MPMC::SoreideWhitsonCorrelations::nitrogenAqueousBip(
                temperatureK, criticalTemperature[nitrogenComponent], salinityMolality);

        // Panfili uses SW for the complete pseudo-component fluid but does not
        // report case-specific aqueous BIPs.  Use the original SW hydrocarbon
        // correlation for the pseudo-components.  C7+ entries are an explicit
        // extrapolation beyond the original light-hydrocarbon fit range.
        if (component >= methaneComponent && component < static_cast<int>(N))
            return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
                temperatureK,
                criticalTemperature[static_cast<std::size_t>(component)],
                acentricFactor[static_cast<std::size_t>(component)],
                salinityMolality);

        throw std::out_of_range("Unexpected Panfili-2025 SW component index.");
    }

    // Current common well API uses one surface density per phase rather than a
    // composition-dependent separator model.  These are transparent engineering
    // closures; they do not enter SW fugacity calculations.
    inline static constexpr std::array<double, 3> surfaceDensity{
        750.0, 1.0, 1000.0};
    inline static constexpr std::array<double, 3> viscosity{
        5.0e-4, 2.0e-5, 5.0e-4};
    static constexpr double waterViscosity = 5.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

// ============================================================================
// 4. Literature reference values and thermodynamically consistent initialization.
// ============================================================================
struct LiteratureReference
{
    // Panfili Table 3: NON-WATER feed composition, normalized to one.
    inline static constexpr std::array<double, 9> hydrocarbonComposition{
        0.01210, 0.01940, 0.65990, 0.08690, 0.05910,
        0.12970, 0.02745, 0.00515, 0.00030};

    static constexpr double initialPressure = 3000.0 * psi;
    static constexpr double referenceDepth = 2910.0 * foot;
    static constexpr double temperature = Fluid::temperature;
    static constexpr double connateWaterSaturation = 0.30;

    static constexpr double depletionYears = 32.0;
    static constexpr double idleYears = 2.0;
    static constexpr double injectionStartYear = depletionYears + idleYears;

    // Section 7.3/Fig.25 shows cumulative injection becoming flat at about
    // year 50, but the text does not give the exact stop time.  This is a
    // figure-derived comparison value, not a verbatim table entry.
    static constexpr double figure25InjectionEndYear = 50.0;

    static constexpr double producerGasRateMscfPerDay = 50000.0;
    static constexpr double producerMinimumBhp = 435.0 * psi;
    static constexpr double injectorGasRateMscfPerDay = 90000.0;
    static constexpr double injectorMaximumBhp = 2850.0 * psi;
};

struct InitialState
{
    static constexpr double pressure = LiteratureReference::initialPressure;
    static constexpr double temperature = Fluid::temperature;

    // The paper specifies 30% CONNATE-WATER SATURATION, not an overall H2O
    // mole fraction.  The value below was obtained by solving the production
    // SW P-T-z flash at 3000 psi / 200 F while preserving the nine-component
    // Table-3 non-water ratios, until Sw=0.30.  This distinction is essential:
    // z_H2O != Sw because phase molar volumes differ strongly.
    //
    // With the v25 SW configuration this gives a G+W initial state with
    // Sg=0.70 and Sw=0.30; a pressure sweep produces liquid condensate below
    // roughly 1.8 ksi, exercising the intended G+W -> O+G+W transition.
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        0.68344517564934826,
        0.0038303133746428858,
        0.0061411635924026441,
        0.20889452858899510,
        0.027508614236071639,
        0.018708390119123518,
        0.041057160718279535,
        0.0086894299284253906,
        0.0016302573454058565,
        0.000094966447305195513};
};

// ============================================================================
// 5. Optional legacy blocks required by common case support.
// ============================================================================
struct Dissolution
{
    static constexpr int component = Fluid::co2Component;
    static constexpr double waterMolarMass = 0.01801528;
    static constexpr double salinityMolality = Fluid::soreideWhitsonSalinityMolality;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};

struct Land { static constexpr double constant = 0.4; };

struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
    inline static constexpr std::array<double, Fluid::N> thetaMax{};
    inline static constexpr std::array<double, Fluid::N> coefficient{};
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
};

struct Numerics
{
    static constexpr double fugacityScalingFactor = 1.0;
    static constexpr bool useVariableBounds = false;
};

// ============================================================================
// 6. Section-7.3 comparison window: 0-100 years.
// ============================================================================
struct Time
{
    // Fig.25/Fig.26 are shown over ~100 years.  One output target per year also
    // aligns the 32-y depletion end, 34-y injection start and ~50-y stop.
    static constexpr int numberOfSteps = 100;
    static constexpr double dtDays = daysPerYear;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-3;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.5;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 14;
    static constexpr int maximumWellControlIterations = 8;
    static constexpr bool printAdaptiveSteps = true;
};

struct Output
{
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 1;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = false;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = true;
    static constexpr bool printComponentMassBalance = true;
    static constexpr bool writeComponentMassBalance = true;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = true;
    static constexpr bool writePhaseStateSnapshots = true;
    static constexpr bool printNewtonIterations = true;
    static constexpr bool writeSolverHistory = true;
    static constexpr bool printWellControlSwitches = true;
    static constexpr bool writeWellControlSwitchHistory = true;
    static constexpr bool printReservoirDiagnostics = true;
    static constexpr bool writeReservoirDiagnostics = true;
    static constexpr bool printFinalSummary = true;
    static constexpr bool writeFinalSummary = true;
};

struct Config final
{
    inline static constexpr const char *name = CaseConfig::name;
    using Model = CaseConfig::Model;
    using Grid = CaseConfig::Grid;
    using Fluid = CaseConfig::Fluid;
    using InitialState = CaseConfig::InitialState;
    using Dissolution = CaseConfig::Dissolution;
    using Land = CaseConfig::Land;
    using Adsorption = CaseConfig::Adsorption;
    using Numerics = CaseConfig::Numerics;
    using Time = CaseConfig::Time;
    using Output = CaseConfig::Output;
    using Rock = CaseConfig::Rock;
};

} // namespace CaseConfig
