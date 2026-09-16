/**
 * @file phase_behavior.hpp
 * @brief 统一描述单相、两相和三相存在性的相行为类型。
 */
#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>

namespace MPMC
{

/** @brief 全组分路径使用的 O/G/W 公共相槽位标识。 */
enum class CompositionalPhase : int
{
    Oil = 0,
    Gas = 1,
    Water = 2
};

/**
 * @brief 油富集、气相和水富集相的运行时存在性掩码。
 *
 * 该状态与 PETSc/EOS 解耦，是 flash、稳定性分析和 Natural 相态转换之间交换的公共协议。
 */
class PhasePresence final
{
public:
    static constexpr std::uint8_t oilBit = 1u << 0;
    static constexpr std::uint8_t gasBit = 1u << 1;
    static constexpr std::uint8_t waterBit = 1u << 2;
    static constexpr std::uint8_t allBits = oilBit | gasBit | waterBit;

    constexpr PhasePresence() noexcept : bits_(allBits) {}
    explicit constexpr PhasePresence(std::uint8_t bits) noexcept : bits_(bits & allBits) {}

    [[nodiscard]] static constexpr PhasePresence all() noexcept { return PhasePresence(allBits); }
    [[nodiscard]] static constexpr PhasePresence oilOnly() noexcept { return PhasePresence(oilBit); }
    [[nodiscard]] static constexpr PhasePresence gasOnly() noexcept { return PhasePresence(gasBit); }
    [[nodiscard]] static constexpr PhasePresence waterOnly() noexcept { return PhasePresence(waterBit); }

    [[nodiscard]] constexpr bool contains(CompositionalPhase phase) const noexcept
    {
        return (bits_ & bit(phase)) != 0;
    }

    constexpr void add(CompositionalPhase phase) noexcept { bits_ |= bit(phase); }
    constexpr void remove(CompositionalPhase phase) noexcept { bits_ &= static_cast<std::uint8_t>(~bit(phase)); }

    [[nodiscard]] constexpr int count() const noexcept
    {
        return (contains(CompositionalPhase::Oil) ? 1 : 0) +
               (contains(CompositionalPhase::Gas) ? 1 : 0) +
               (contains(CompositionalPhase::Water) ? 1 : 0);
    }

    [[nodiscard]] constexpr std::uint8_t bits() const noexcept { return bits_; }
    [[nodiscard]] constexpr bool empty() const noexcept { return bits_ == 0; }

    [[nodiscard]] static constexpr std::uint8_t bit(CompositionalPhase phase) noexcept
    {
        switch (phase)
        {
        case CompositionalPhase::Oil: return oilBit;
        case CompositionalPhase::Gas: return gasBit;
        case CompositionalPhase::Water: return waterBit;
        }
        return 0;
    }

private:
    std::uint8_t bits_;
};

[[nodiscard]] constexpr int phaseIndex(CompositionalPhase phase) noexcept
{
    return static_cast<int>(phase);
}

[[nodiscard]] inline PhasePresence decodePhasePresence(double flag)
{
    const int raw = static_cast<int>(flag + (flag >= 0.0 ? 0.5 : -0.5));
    if (raw <= 0 || raw > static_cast<int>(PhasePresence::allBits))
        throw std::runtime_error("Invalid fully-compositional phase-presence flag.");
    return PhasePresence(static_cast<std::uint8_t>(raw));
}

[[nodiscard]] constexpr double encodePhasePresence(PhasePresence presence) noexcept
{
    return static_cast<double>(presence.bits());
}

} // namespace MPMC
