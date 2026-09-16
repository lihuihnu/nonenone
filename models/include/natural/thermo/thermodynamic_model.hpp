/**
 * @file thermodynamic_model.hpp
 * @brief PR、SW、CPA 热力学后端的公共模型枚举与选择。
 */
#pragma once

namespace MPMC
{

/** @brief 全组分相可选择的热力学后端；上层 flash/状态接口保持一致。 */
enum class CubicThermodynamicModel
{
    PengRobinson = 0,
    SoreideWhitson = 1,
    CubicPlusAssociation = 2
};

/** @brief CPA 缔合项下方采用的立方 EOS 物理贡献类型。 */
enum class CpaCubicPhysicalTerm
{
    SoaveRedlichKwong = 0,
    PengRobinson = 1
};

[[nodiscard]] inline const char *thermodynamicModelName(
    CubicThermodynamicModel model) noexcept
{
    switch (model)
    {
    case CubicThermodynamicModel::PengRobinson: return "Peng-Robinson (PR)";
    case CubicThermodynamicModel::SoreideWhitson: return "Soreide-Whitson (SW/PR)";
    case CubicThermodynamicModel::CubicPlusAssociation: return "Cubic-Plus-Association (CPA)";
    }
    return "Unknown";
}

} // namespace MPMC
