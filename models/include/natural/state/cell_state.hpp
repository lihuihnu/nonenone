/**
 * @file cell_state.hpp
 * @brief 一个网格单元的主变量、相态与历史状态聚合。
 */
#pragma once

#include <natural/phase_state.hpp>

#include <array>

namespace MPMC
{

/**
 * @brief 一个单元的主状态快照，不含任何 PETSc/Grid 所有权。
 */
template <class Indices, class Scalar>
struct CellState
{
    static constexpr int numComponents = Indices::numComponents;
    static constexpr int numPhases = Indices::numPhases;

    Scalar pressure{0.0};
    Scalar liquidSaturation{0.0};
    Scalar vaporSaturation{0.0};
    Scalar waterSaturation{0.0};
    Scalar wellPressure{0.0};
    Scalar aqueousCO2MoleFraction{0.0};

    std::array<Scalar, numComponents> liquidMoleFraction{};
    std::array<Scalar, numComponents> vaporMoleFraction{};
    std::array<Scalar, numComponents> aqueousMoleFraction{};

    HydrocarbonPhaseState hydrocarbonPhaseState{HydrocarbonPhaseState::TwoPhase};
    PhasePresence phasePresence = PhasePresence::all();
    PhasePresence phaseSuppression{std::uint8_t{0}};
};

/**
 * @brief 由 CellState 计算得到的局部物性。
 */
template <class Indices, class Scalar>
struct CellProperties
{
    static constexpr int numComponents = Indices::numComponents;
    static constexpr int numPhases = Indices::numPhases;

    std::array<Scalar, numPhases> density{};
    std::array<Scalar, numPhases> viscosity{};
    std::array<Scalar, numPhases> mobility{};
    std::array<Scalar, numPhases> saturation{};

    std::array<std::array<Scalar, numComponents>, numPhases> massFraction{};
    std::array<std::array<Scalar, numComponents>, Indices::numThermodynamicPhases> fugacity{};

    Scalar porosity{0.0};
    Scalar aqueousCO2MassFraction{0.0};
    Scalar aqueousCO2Fugacity{0.0};
    Scalar trappedGasSaturation{0.0};
    Scalar freeGasSaturation{0.0};
    std::array<Scalar, numComponents> adsorbedVolume{};
};

} // namespace MPMC
