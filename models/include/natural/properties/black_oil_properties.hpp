/**
 * @file black_oil_properties.hpp
 * @brief 表格式/黑油型相密度和黏度物性。
 */
#pragma once

#include <common/math.hpp>

#include <array>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace MPMC
{

/**
 * @brief K-value 路径使用的表驱动液/气密度与黏度模型。
 *
 * Each callback receives `(pressure, compositionSelector)`. The selector is
 * component 0 for the liquid table and component 1 for the vapor table,
 * matching the existing MPMC black-oil/K-value parameterization. The first
 * array slot is liquid, the second vapor.
 */
template <class Indices>
class BlackOilPropertyModel
{
public:
    using ValueType = typename Indices::ValueType;
    using TableFunction =
        std::function<ValueType(const ValueType &, const ValueType &)>;
    using PhaseFunctions = std::array<TableFunction, 2>;

    /** @brief 一次配置完整的液/气密度和黏度表。 */
    void configure(PhaseFunctions density, PhaseFunctions viscosity)
    {
        densityFunctions_ = std::move(density);
        viscosityFunctions_ = std::move(viscosity);
        validate_();
    }

    /** @brief 在储层压力下计算已配置的液/气密度量。 */
    template <class Scalar>
    [[nodiscard]] Scalar computeMolarDensity(
        const Scalar &pressure,
        const std::array<Scalar, Indices::numComponents> &composition,
        bool liquid) const
    {
        validate_();
        const std::size_t phase = liquid ? 0u : 1u;
        const std::size_t selector =
            (liquid || Indices::numComponents == 1) ? 0u : 1u;
        const ValueType result = densityFunctions_[phase](
            ValueType(pressure), ValueType(composition[selector]));
        if constexpr (std::is_same_v<Scalar, double>)
        {
            return scalarValue(result);
        }
        else
        {
            return Scalar(result);
        }
    }

    /** @brief 在储层压力下计算液/气黏度 [Pa s]。 */
    template <class Scalar>
    [[nodiscard]] Scalar computeViscosity(
        const Scalar &pressure,
        const std::array<Scalar, Indices::numComponents> &composition,
        bool liquid) const
    {
        validate_();
        const std::size_t phase = liquid ? 0u : 1u;
        const std::size_t selector =
            (liquid || Indices::numComponents == 1) ? 0u : 1u;
        const ValueType result = viscosityFunctions_[phase](
            ValueType(pressure), ValueType(composition[selector]));
        if constexpr (std::is_same_v<Scalar, double>)
        {
            return scalarValue(result);
        }
        else
        {
            return Scalar(result);
        }
    }

private:
    void validate_() const
    {
        if (!densityFunctions_[0] || !densityFunctions_[1] ||
            !viscosityFunctions_[0] || !viscosityFunctions_[1])
        {
            throw std::logic_error(
                "BlackOilPropertyModel requires liquid and vapor density/viscosity callbacks.");
        }
    }

    PhaseFunctions densityFunctions_{};
    PhaseFunctions viscosityFunctions_{};
};

} // namespace MPMC
