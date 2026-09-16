/**
 * @file compositional_mixture.hpp
 * @brief 组分临界性质、分子量和 EOS 二元作用参数容器。
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <utility>

namespace MPMC
{

/**
 * @brief 固定组分数的混合物基础参数。
 *
 * 所有内部存储均为 std::array，避免原始 C 数组的尺寸/拷贝隐患。
 * 单位约定：Tc[K]、Pc[Pa]、Vc[m^3/mol]、M[kg/mol]。
 */
template <class Indices>
class CompositionalMixture final
{
public:
    static constexpr int numComponents = Indices::numComponents;
    using ScalarArray = std::array<double, numComponents>;
    using BinaryInteractionMatrix =
        std::array<std::array<double, numComponents>, numComponents>;

    CompositionalMixture() = default;

    CompositionalMixture(
        ScalarArray criticalTemperature,
        ScalarArray criticalPressure,
        ScalarArray criticalVolume,
        ScalarArray acentricFactor,
        ScalarArray molarMass,
        BinaryInteractionMatrix binaryInteraction)
        : criticalTemperature_(std::move(criticalTemperature)),
          criticalPressure_(std::move(criticalPressure)),
          criticalVolume_(std::move(criticalVolume)),
          acentricFactor_(std::move(acentricFactor)),
          molarMass_(std::move(molarMass)),
          binaryInteraction_(std::move(binaryInteraction))
    {
        validate();
    }

    /** @brief 便于在算例参数表中直接使用 initializer_list 构造固定组分数据。 */
    CompositionalMixture(
        std::initializer_list<double> criticalTemperature,
        std::initializer_list<double> criticalPressure,
        std::initializer_list<double> criticalVolume,
        std::initializer_list<double> acentricFactor,
        std::initializer_list<double> molarMass,
        std::initializer_list<std::initializer_list<double>> binaryInteraction)
    {
        assign(criticalTemperature, criticalTemperature_);
        assign(criticalPressure, criticalPressure_);
        assign(criticalVolume, criticalVolume_);
        assign(acentricFactor, acentricFactor_);
        assign(molarMass, molarMass_);

        if (binaryInteraction.size() != static_cast<std::size_t>(numComponents))
            throw std::invalid_argument("Binary-interaction matrix row count is invalid.");
        std::size_t i = 0;
        for (const auto &row : binaryInteraction)
        {
            if (row.size() != static_cast<std::size_t>(numComponents))
                throw std::invalid_argument("Binary-interaction matrix must be square.");
            std::size_t j = 0;
            for (double value : row)
                binaryInteraction_[i][j++] = value;
            ++i;
        }
        validate();
    }

    /** @brief 返回组分 i 的临界温度 T_c [K]。 */
    [[nodiscard]] double criticalTemperature(int i) const
    {
        return criticalTemperature_.at(static_cast<std::size_t>(i));
    }

    /** @brief 返回组分 i 的临界压力 P_c [Pa]。 */
    [[nodiscard]] double criticalPressure(int i) const
    {
        return criticalPressure_.at(static_cast<std::size_t>(i));
    }

    /** @brief 返回组分 i 的临界摩尔体积 V_c [m3/mol]。 */
    [[nodiscard]] double criticalVolume(int i) const
    {
        return criticalVolume_.at(static_cast<std::size_t>(i));
    }

    /** @brief 返回组分 i 的 Pitzer 偏心因子 [-]。 */
    [[nodiscard]] double acentricFactor(int i) const
    {
        return acentricFactor_.at(static_cast<std::size_t>(i));
    }

    /** @brief 返回组分 i 的摩尔质量 [kg/mol]。 */
    [[nodiscard]] double molecularWeight(int i) const
    {
        return molarMass_.at(static_cast<std::size_t>(i));
    }

    /** @brief 返回 EOS 混合规则使用的二元交互系数 k_ij。 */
    [[nodiscard]] double binaryInteractionCoefficient(int i, int j) const
    {
        return binaryInteraction_.at(static_cast<std::size_t>(i))
            .at(static_cast<std::size_t>(j));
    }

private:
    static void assign(
        std::initializer_list<double> source,
        ScalarArray &destination)
    {
        if (source.size() != destination.size())
            throw std::invalid_argument("Mixture property count does not match Indices::numComponents.");
        std::copy(source.begin(), source.end(), destination.begin());
    }

    void validate() const
    {
        for (int i = 0; i < numComponents; ++i)
        {
            if (!(criticalTemperature(i) > 0.0) || !std::isfinite(criticalTemperature(i)))
                throw std::invalid_argument("Critical temperatures must be finite and positive.");
            if (!(criticalPressure(i) > 0.0) || !std::isfinite(criticalPressure(i)))
                throw std::invalid_argument("Critical pressures must be finite and positive.");
            if (!(criticalVolume(i) > 0.0) || !std::isfinite(criticalVolume(i)))
                throw std::invalid_argument("Critical volumes must be finite and positive.");
            if (!std::isfinite(acentricFactor(i)))
                throw std::invalid_argument("Acentric factors must be finite.");
            if (!(molecularWeight(i) > 0.0) || !std::isfinite(molecularWeight(i)))
                throw std::invalid_argument("Molar masses must be finite and positive.");
            for (int j = 0; j < numComponents; ++j)
            {
                const double kij = binaryInteractionCoefficient(i, j);
                if (!std::isfinite(kij))
                    throw std::invalid_argument("Binary interaction coefficients must be finite.");
                if (j > i)
                {
                    const double kji = binaryInteractionCoefficient(j, i);
                    const double scale = std::max({1.0, std::abs(kij), std::abs(kji)});
                    if (std::abs(kij - kji) > 1.0e-12 * scale)
                        throw std::invalid_argument(
                            "Binary interaction matrix must satisfy k_ij = k_ji.");
                }
            }
        }
    }

    ScalarArray criticalTemperature_{};
    ScalarArray criticalPressure_{};
    ScalarArray criticalVolume_{};
    ScalarArray acentricFactor_{};
    ScalarArray molarMass_{};
    BinaryInteractionMatrix binaryInteraction_{};
};

} // namespace MPMC
