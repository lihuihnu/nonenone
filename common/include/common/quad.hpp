/**
 * @file quad.hpp
 * @brief 小型定长四元数据结构及基础运算。
 */
#pragma once

#include <common/math.hpp>

#if defined(HAVE_QUAD) && HAVE_QUAD

extern "C"
{
#include <quadmath.h>
}

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace MPMC
{

/**
 * @brief GCC `__float128` 的显式项目别名。
 *
 * 本文件不修改 `std::is_floating_point`、
 * `std::is_arithmetic` 等标准 type traits。
 */
using Quad = __float128;

/**
 * @brief `MathToolbox<Quad>` 的显式数值扩展。
 */
template <>
struct MathToolbox<Quad>
{
    using Scalar = Quad;
    using ValueType = Quad;
    using InnerToolbox = MathToolbox<Quad>;

    [[nodiscard]] static constexpr Quad value(Quad input) noexcept
    {
        return input;
    }

    [[nodiscard]] static constexpr Quad scalarValue(Quad input) noexcept
    {
        return input;
    }

    [[nodiscard]] static constexpr Quad createBlank(Quad) noexcept
    {
        return static_cast<Quad>(0);
    }

    [[nodiscard]] static constexpr Quad createConstant(Quad input) noexcept
    {
        return input;
    }

    [[nodiscard]] static Quad createConstant(
        unsigned numDerivatives,
        Quad input)
    {
        if (numDerivatives != 0)
        {
            throw std::logic_error(
                "Quad scalar cannot store derivatives.");
        }

        return input;
    }

    [[nodiscard]] static constexpr Quad createConstant(
        Quad,
        Quad input) noexcept
    {
        return input;
    }

    [[nodiscard]] static Quad createVariable(
        Quad,
        unsigned)
    {
        throw std::logic_error(
            "Quad scalar cannot represent an AD variable.");
    }

    [[nodiscard]] static Quad createVariable(
        unsigned,
        Quad,
        unsigned)
    {
        throw std::logic_error(
            "Quad scalar cannot represent an AD variable.");
    }

    [[nodiscard]] static Quad createVariable(
        Quad,
        Quad,
        unsigned)
    {
        throw std::logic_error(
            "Quad scalar cannot represent an AD variable.");
    }

    template <class LhsEval>
    [[nodiscard]] static LhsEval decay(Quad input)
    {
        return static_cast<LhsEval>(input);
    }

    [[nodiscard]] static bool isSame(
        Quad a,
        Quad b,
        Quad tolerance) noexcept
    {
        const Quad difference =
            fabsq(a - b);

        Quad scale =
            fabsq(a);

        if (fabsq(b) > scale)
        {
            scale = fabsq(b);
        }

        if (scale < static_cast<Quad>(1))
        {
            scale = static_cast<Quad>(1);
        }

        return difference <= tolerance ||
               difference <= tolerance * scale;
    }

    [[nodiscard]] static Quad max(Quad a, Quad b) noexcept
    {
        return a > b ? a : b;
    }

    [[nodiscard]] static Quad min(Quad a, Quad b) noexcept
    {
        return a < b ? a : b;
    }

    [[nodiscard]] static Quad abs(Quad input) noexcept
    {
        return fabsq(input);
    }

    [[nodiscard]] static Quad tan(Quad input) { return tanq(input); }
    [[nodiscard]] static Quad atan(Quad input) { return atanq(input); }
    [[nodiscard]] static Quad atan2(Quad y, Quad x) { return atan2q(y, x); }
    [[nodiscard]] static Quad sin(Quad input) { return sinq(input); }
    [[nodiscard]] static Quad asin(Quad input) { return asinq(input); }
    [[nodiscard]] static Quad sinh(Quad input) { return sinhq(input); }
    [[nodiscard]] static Quad asinh(Quad input) { return asinhq(input); }
    [[nodiscard]] static Quad cos(Quad input) { return cosq(input); }
    [[nodiscard]] static Quad acos(Quad input) { return acosq(input); }
    [[nodiscard]] static Quad cosh(Quad input) { return coshq(input); }
    [[nodiscard]] static Quad acosh(Quad input) { return acoshq(input); }
    [[nodiscard]] static Quad sqrt(Quad input) { return sqrtq(input); }
    [[nodiscard]] static Quad exp(Quad input) { return expq(input); }
    [[nodiscard]] static Quad log(Quad input) { return logq(input); }
    [[nodiscard]] static Quad log10(Quad input) { return log10q(input); }
    [[nodiscard]] static Quad pow(Quad base, Quad exponent)
    {
        return powq(base, exponent);
    }

    [[nodiscard]] static bool isfinite(Quad input) noexcept
    {
        return finiteq(input) != 0;
    }

    [[nodiscard]] static bool isnan(Quad input) noexcept
    {
        return isnanq(input) != 0;
    }
};

} // namespace MPMC

#endif
