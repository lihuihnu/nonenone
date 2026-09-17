#pragma once

#include <array>

/**
 * @file bsb_reference_properties.hpp
 * @brief BSB pseudo-component reference properties used by the uploaded LMH cases.
 *
 * Values are the same preregistered H2O/BSB properties already present in the
 * H2O-CO2-BSB laboratory benchmark.  The LMH case removes XH(C28+) and
 * renormalizes only the L/M/H hydrocarbon mole fractions.
 */
namespace BsbReference
{

struct ComponentProperties
{
    double criticalTemperatureK;
    double criticalPressurePa;
    double criticalVolumeM3PerMol;
    double acentricFactor;
    double molarMassKgPerMol;
};

inline constexpr ComponentProperties water{
    647.30, 22.048e6, 5.594803743e-5, 0.344, 0.018015};
inline constexpr ComponentProperties light{
    354.1916431226766, 4.065799256505576e6, 1.95564637471291e-4,
    0.1498936802973978, 0.04606110037174722};
inline constexpr ComponentProperties middle{
    605.78, 2.175e6, 6.252498825299838e-4, 0.618, 0.140960};
inline constexpr ComponentProperties heavy{
    751.00, 1.654e6, 1.0193008374141067e-3, 0.957, 0.280990};
inline constexpr ComponentProperties extraHeavy{
    942.50, 1.642e6, 1.2885644791440596e-3, 1.268, 0.519620};

inline constexpr double waterHydrocarbonKij = 0.5;
inline constexpr double initialWaterMoleFraction = 0.20;

// BSB L/M/H mole ratios after excluding XH(C28+) and renormalizing L+M+H.
inline constexpr std::array<double, 3> lmhMoleFraction{
    0.45083798882681564,
    0.36916201117318430,
    0.18000000000000000};

} // namespace BsbReference
