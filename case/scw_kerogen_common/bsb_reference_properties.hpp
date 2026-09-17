#pragma once

#include <array>
#include <cstddef>

/**
 * @file bsb_reference_properties.hpp
 * @brief H2O-CO2-BSB 算例的水和烃拟组分参考物性。
 *
 * 数值逐项来自 h2o_co2_bsb_lumped_3d_lab/component_lumping.csv。当前
 * SCW 算例直接使用表中的 L、M、H 物性，不用全馏分平均值冒充重组分。
 */
namespace BsbReference
{

struct Component
{
    double criticalTemperatureK;
    double criticalPressurePa;
    double criticalVolumeM3PerMol;
    double acentricFactor;
    double molarMassKgPerMol;
};

inline constexpr Component water{
    647.30, 22.048e6, 5.594803743e-5, 0.344, 0.018015};
inline constexpr Component light{
    354.1916431226766, 4.065799256505576e6,
    1.95564637471291e-4, 0.1498936802973978,
    0.04606110037174722};
inline constexpr Component middle{
    605.78, 2.175e6, 6.252498825299838e-4, 0.618, 0.140960};
inline constexpr Component heavy{
    751.00, 1.654e6, 1.0193008374141067e-3, 0.957, 0.280990};
inline constexpr Component extraHeavy{
    942.50, 1.642e6, 1.2885644791440596e-3, 1.268, 0.519620};

// Hydrocarbon-normalized mole fractions from the BSB reference case.
inline constexpr std::array<double, 4> hydrocarbonMoleFraction{
    0.41757218255200246,
    0.34192279830280453,
    0.16671841043154300,
    0.07378660871365000};

// Renormalized L/M/H distribution after deliberately excluding the separate
// XH(C28+) class from the requested three-component L+M+H experiment.
inline constexpr std::array<double, 3> lmhMoleFraction{
    0.45083798882681564,
    0.36916201117318440,
    0.18000000000000002};

inline constexpr double waterHydrocarbonKij = 0.5000;
inline constexpr double initialWaterMoleFraction = 0.20;

} // namespace BsbReference
