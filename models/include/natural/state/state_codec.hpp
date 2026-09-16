/**
 * @file state_codec.hpp
 * @brief Natural 主变量向量与 CellState 之间的编码/解码。
 */
#pragma once

#include <common/math.hpp>
#include <natural/numerics.hpp>
#include <natural/primary_variables.hpp>
#include <natural/state/cell_state.hpp>
#include <natural/state/phase_equilibrium.hpp>

#include <array>
#include <cstddef>

namespace MPMC
{

/**
 * @brief 主未知量数组与强类型 CellState 之间的唯一转换入口。
 *
 * 这里集中维护 dependent-last-component 约束，避免 EOS、通量、井等模块各自
 * 重复 `x_N = 1 - sum(x_i)` 并产生不一致。
 */
template <class Indices>
class CellStateCodec final
{
public:
    using Scalar = typename Indices::ValueType;
    using PrimaryArray =
        std::array<double, Indices::numPrimaryVariables>;
    using State = CellState<Indices, Scalar>;

    [[nodiscard]] static State decode(
        const PrimaryArray &primary,
        HydrocarbonPhaseState phaseState)
    {
        State state;

        state.hydrocarbonPhaseState =
            phaseState;

        state.pressure =
            PrimaryVariables<Indices>::make(
                primary[Indices::Primary::pressure],
                Indices::Primary::pressure);

        state.liquidSaturation =
            PrimaryVariables<Indices>::make(
                primary[Indices::Primary::liquidSaturation],
                Indices::Primary::liquidSaturation);

        state.vaporSaturation =
            PrimaryVariables<Indices>::make(
                primary[Indices::Primary::vaporSaturation],
                Indices::Primary::vaporSaturation);

        if constexpr (Indices::hasWater)
        {
            state.waterSaturation =
                PrimaryVariables<Indices>::make(
                    primary[Indices::Primary::waterSaturation],
                    Indices::Primary::waterSaturation);
        }

        if constexpr (Indices::hasWellUnknown)
        {
            state.wellPressure =
                PrimaryVariables<Indices>::make(
                    primary[Indices::Primary::wellPressure],
                    Indices::Primary::wellPressure);
        }

        if constexpr (Indices::hasAqueousCO2Dissolution)
        {
            state.aqueousCO2MoleFraction =
                PrimaryVariables<Indices>::make(
                    primary[
                        Indices::Primary::
                            aqueousCO2MoleFraction],
                    Indices::Primary::
                        aqueousCO2MoleFraction);
        }

        decodeComposition_(
            primary,
            Indices::Primary::liquidComposition,
            state.liquidMoleFraction);

        decodeComposition_(
            primary,
            Indices::Primary::vaporComposition,
            state.vaporMoleFraction);

        return state;
    }

    [[nodiscard]] static State decode(
        const PrimaryArray &primary,
        const PhaseStateData<Indices> &phaseState)
    {
        if constexpr (!Indices::fullyCompositionalThreePhase)
        {
            return decode(primary, phaseState.phase);
        }
        else
        {
            State state;
            state.phasePresence = phaseState.phasePresence;
            state.phaseSuppression = phaseState.phaseSuppression;
            state.pressure = PrimaryVariables<Indices>::make(
                primary[Indices::Primary::pressure], Indices::Primary::pressure);
            state.liquidSaturation = PrimaryVariables<Indices>::make(
                primary[Indices::Primary::liquidSaturation], Indices::Primary::liquidSaturation);
            state.vaporSaturation = PrimaryVariables<Indices>::make(
                primary[Indices::Primary::vaporSaturation], Indices::Primary::vaporSaturation);
            state.waterSaturation = PrimaryVariables<Indices>::make(
                primary[Indices::Primary::waterSaturation], Indices::Primary::waterSaturation);
            if constexpr (Indices::hasWellUnknown)
            {
                state.wellPressure = PrimaryVariables<Indices>::make(
                    primary[Indices::Primary::wellPressure], Indices::Primary::wellPressure);
            }
            decodeComposition_(primary, Indices::Primary::liquidComposition, state.liquidMoleFraction);
            decodeComposition_(primary, Indices::Primary::vaporComposition, state.vaporMoleFraction);
            decodeComposition_(primary, Indices::Primary::waterComposition, state.aqueousMoleFraction);
            return state;
        }
    }

private:
    template <std::size_t Size>
    static void decodeComposition_(
        const PrimaryArray &primary,
        const std::array<int, Size> &indices,
        std::array<Scalar, Indices::numComponents> &composition)
    {
        Scalar dependent = 1.0;

        for (std::size_t component = 0;
             component < Size;
             ++component)
        {
            const int primaryIndex =
                indices[component];

            composition[component] =
                PrimaryVariables<Indices>::make(
                    primary[
                        static_cast<std::size_t>(
                            primaryIndex)],
                    primaryIndex);

            dependent -=
                composition[component];
        }

        // 数值：当真实 dependent component 位于 trace 边界（例如 PTz flash
        // 给出的 1e-30）时，浮点消去可能产生极小负值。这里保留 AD 导数，仅把
        // 标量值平移回物理边界 x=0，避免把舍入伪影送入 cubic EOS；明显的负值
        // 则保留给上层 nonlinear sanitizer 识别和修正。
        if constexpr (Indices::fullyCompositionalThreePhase)
        {
            const double dependentValue = scalarValue(dependent);
            if (dependentValue < 0.0 &&
                dependentValue >= -NaturalNumerics::phaseEquilibriumTraceComposition)
            {
                dependent += -dependentValue;
            }
        }

        composition.back() = dependent;
    }
};

} // namespace MPMC
