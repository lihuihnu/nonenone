/**
 * @file math.hpp
 * @brief 项目通用数值数学辅助函数。
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace MPMC
{

/**
 * @brief 数值类型统一适配器。
 *
 * MathToolbox 是 common 与 AD 模块之间唯一保留的数值扩展点。
 *
 * 对普通浮点类型，本模板直接转发到 `<cmath>`；自动微分类型通过在 AD 模块中
 * 对 `MathToolbox<Evaluation<...>>` 做特化，使 EOS、插值和物性代码能够使用
 * 同一组数学接口。
 *
 * @tparam ScalarT 普通浮点标量类型。
 */
template <class ScalarT>
struct MathToolbox
{
    static_assert(
        std::is_floating_point_v<ScalarT>,
        "MathToolbox primary template requires a built-in floating type. "
        "Provide an explicit specialization for custom numeric types.");

    using Scalar = ScalarT;
    using ValueType = ScalarT;
    using InnerToolbox = MathToolbox<Scalar>;

    [[nodiscard]] static constexpr Scalar value(Scalar input) noexcept
    {
        return input;
    }

    [[nodiscard]] static constexpr Scalar scalarValue(Scalar input) noexcept
    {
        return input;
    }

    [[nodiscard]] static constexpr Scalar createBlank(Scalar) noexcept
    {
        return Scalar{};
    }

    [[nodiscard]] static constexpr Scalar createConstant(Scalar input) noexcept
    {
        return input;
    }

    [[nodiscard]] static Scalar createConstant(
        unsigned numDerivatives,
        Scalar input)
    {
        if (numDerivatives != 0)
        {
            throw std::logic_error(
                "A primitive floating scalar cannot store derivatives.");
        }

        return input;
    }

    [[nodiscard]] static constexpr Scalar createConstant(
        Scalar,
        Scalar input) noexcept
    {
        return input;
    }

    [[nodiscard]] static Scalar createVariable(
        Scalar,
        unsigned)
    {
        throw std::logic_error(
            "A primitive floating scalar cannot represent an AD variable.");
    }

    [[nodiscard]] static Scalar createVariable(
        unsigned,
        Scalar,
        unsigned)
    {
        throw std::logic_error(
            "A primitive floating scalar cannot represent an AD variable.");
    }

    [[nodiscard]] static Scalar createVariable(
        Scalar,
        Scalar,
        unsigned)
    {
        throw std::logic_error(
            "A primitive floating scalar cannot represent an AD variable.");
    }

    template <class LhsEval>
    [[nodiscard]] static LhsEval decay(Scalar input)
    {
        static_assert(
            std::is_floating_point_v<LhsEval>,
            "Primitive MathToolbox can only decay to another floating type.");

        return static_cast<LhsEval>(input);
    }

    /**
     * @brief 同时使用绝对和相对容差比较两个数。
     *
     * 相对尺度使用 `max(1, |a|, |b|)`，因此在 `a≈-b` 时仍保持稳定尺度。
     */
    [[nodiscard]] static bool isSame(
        Scalar a,
        Scalar b,
        Scalar tolerance) noexcept
    {
        const Scalar difference =
            std::abs(a - b);

        const Scalar scale =
            std::max(
                Scalar{1},
                std::max(
                    std::abs(a),
                    std::abs(b)));

        return difference <= tolerance ||
               difference <= tolerance * scale;
    }

    [[nodiscard]] static Scalar max(Scalar a, Scalar b) noexcept
    {
        return std::max(a, b);
    }

    [[nodiscard]] static Scalar min(Scalar a, Scalar b) noexcept
    {
        return std::min(a, b);
    }

    [[nodiscard]] static Scalar abs(Scalar input) noexcept
    {
        using std::abs;
        return abs(input);
    }

    [[nodiscard]] static Scalar tan(Scalar input)
    {
        using std::tan;
        return tan(input);
    }

    [[nodiscard]] static Scalar atan(Scalar input)
    {
        using std::atan;
        return atan(input);
    }

    [[nodiscard]] static Scalar atan2(Scalar y, Scalar x)
    {
        using std::atan2;
        return atan2(y, x);
    }

    [[nodiscard]] static Scalar sin(Scalar input)
    {
        using std::sin;
        return sin(input);
    }

    [[nodiscard]] static Scalar asin(Scalar input)
    {
        using std::asin;
        return asin(input);
    }

    [[nodiscard]] static Scalar sinh(Scalar input)
    {
        using std::sinh;
        return sinh(input);
    }

    [[nodiscard]] static Scalar asinh(Scalar input)
    {
        using std::asinh;
        return asinh(input);
    }

    [[nodiscard]] static Scalar cos(Scalar input)
    {
        using std::cos;
        return cos(input);
    }

    [[nodiscard]] static Scalar acos(Scalar input)
    {
        using std::acos;
        return acos(input);
    }

    [[nodiscard]] static Scalar cosh(Scalar input)
    {
        using std::cosh;
        return cosh(input);
    }

    [[nodiscard]] static Scalar acosh(Scalar input)
    {
        using std::acosh;
        return acosh(input);
    }

    [[nodiscard]] static Scalar sqrt(Scalar input)
    {
        using std::sqrt;
        return sqrt(input);
    }

    [[nodiscard]] static Scalar exp(Scalar input)
    {
        using std::exp;
        return exp(input);
    }

    [[nodiscard]] static Scalar log(Scalar input)
    {
        using std::log;
        return log(input);
    }

    [[nodiscard]] static Scalar log10(Scalar input)
    {
        using std::log10;
        return log10(input);
    }

    [[nodiscard]] static Scalar pow(Scalar base, Scalar exponent)
    {
        using std::pow;
        return pow(base, exponent);
    }

    [[nodiscard]] static bool isfinite(Scalar input) noexcept
    {
        return std::isfinite(input);
    }

    [[nodiscard]] static bool isnan(Scalar input) noexcept
    {
        return std::isnan(input);
    }
};

/**
 * @brief 推导两个数值表达式混合运算后的 Evaluation 类型。
 *
 * 该类型萃取用于 DenseAd 混合表达式：先去掉 cv/ref 修饰，
 * 再选择可由另一侧构造的 Evaluation 类型。
 */
template <class Eval1, class Eval2>
struct ReturnEval_
{
    using T = std::remove_cv_t<std::remove_reference_t<Eval1>>;
    using U = std::remove_cv_t<std::remove_reference_t<Eval2>>;

    using type =
        std::conditional_t<
            std::is_constructible_v<T, U>,
            T,
            U>;
};

template <class Evaluation>
[[nodiscard]] Evaluation blank(const Evaluation &x)
{
    return MathToolbox<Evaluation>::createBlank(x);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation constant(const Scalar &input)
{
    return MathToolbox<Evaluation>::createConstant(input);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation constant(
    unsigned numDerivatives,
    const Scalar &input)
{
    return MathToolbox<Evaluation>::createConstant(
        numDerivatives,
        input);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation constant(
    const Evaluation &reference,
    const Scalar &input)
{
    return MathToolbox<Evaluation>::createConstant(
        reference,
        input);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation variable(
    unsigned numDerivatives,
    const Scalar &input,
    unsigned derivativeIndex)
{
    return MathToolbox<Evaluation>::createVariable(
        numDerivatives,
        input,
        derivativeIndex);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation variable(
    const Evaluation &reference,
    const Scalar &input,
    unsigned derivativeIndex)
{
    return MathToolbox<Evaluation>::createVariable(
        reference,
        input,
        derivativeIndex);
}

template <class Evaluation, class Scalar>
[[nodiscard]] Evaluation variable(
    const Scalar &input,
    unsigned derivativeIndex)
{
    return MathToolbox<Evaluation>::createVariable(
        input,
        derivativeIndex);
}

template <class ResultEval, class Evaluation>
[[nodiscard]] auto decay(const Evaluation &input)
    -> decltype(
        MathToolbox<Evaluation>::template decay<ResultEval>(
            input))
{
    return MathToolbox<Evaluation>::template decay<ResultEval>(
        input);
}

template <class Evaluation>
[[nodiscard]] auto getValue(const Evaluation &input)
    -> decltype(MathToolbox<Evaluation>::value(input))
{
    return MathToolbox<Evaluation>::value(input);
}

template <class Evaluation>
[[nodiscard]] auto scalarValue(const Evaluation &input)
    -> decltype(MathToolbox<Evaluation>::scalarValue(input))
{
    return MathToolbox<Evaluation>::scalarValue(input);
}

template <class Evaluation1, class Evaluation2>
[[nodiscard]] typename ReturnEval_<Evaluation1, Evaluation2>::type
max(const Evaluation1 &a, const Evaluation2 &b)
{
    using Result =
        typename ReturnEval_<Evaluation1, Evaluation2>::type;

    return MathToolbox<Result>::max(a, b);
}

template <class Evaluation1, class Evaluation2>
[[nodiscard]] typename ReturnEval_<Evaluation1, Evaluation2>::type
min(const Evaluation1 &a, const Evaluation2 &b)
{
    using Result =
        typename ReturnEval_<Evaluation1, Evaluation2>::type;

    return MathToolbox<Result>::min(a, b);
}

template <class Evaluation>
[[nodiscard]] Evaluation abs(const Evaluation &input)
{
    return MathToolbox<Evaluation>::abs(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation tan(const Evaluation &input)
{
    return MathToolbox<Evaluation>::tan(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation atan(const Evaluation &input)
{
    return MathToolbox<Evaluation>::atan(input);
}

template <class Evaluation1, class Evaluation2>
[[nodiscard]] typename ReturnEval_<Evaluation1, Evaluation2>::type
atan2(const Evaluation1 &y, const Evaluation2 &x)
{
    using Result =
        typename ReturnEval_<Evaluation1, Evaluation2>::type;

    return MathToolbox<Result>::atan2(y, x);
}

template <class Evaluation>
[[nodiscard]] Evaluation sin(const Evaluation &input)
{
    return MathToolbox<Evaluation>::sin(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation asin(const Evaluation &input)
{
    return MathToolbox<Evaluation>::asin(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation sinh(const Evaluation &input)
{
    return MathToolbox<Evaluation>::sinh(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation asinh(const Evaluation &input)
{
    return MathToolbox<Evaluation>::asinh(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation cos(const Evaluation &input)
{
    return MathToolbox<Evaluation>::cos(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation acos(const Evaluation &input)
{
    return MathToolbox<Evaluation>::acos(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation cosh(const Evaluation &input)
{
    return MathToolbox<Evaluation>::cosh(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation acosh(const Evaluation &input)
{
    return MathToolbox<Evaluation>::acosh(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation sqrt(const Evaluation &input)
{
    return MathToolbox<Evaluation>::sqrt(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation exp(const Evaluation &input)
{
    return MathToolbox<Evaluation>::exp(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation log(const Evaluation &input)
{
    return MathToolbox<Evaluation>::log(input);
}

template <class Evaluation>
[[nodiscard]] Evaluation log10(const Evaluation &input)
{
    return MathToolbox<Evaluation>::log10(input);
}

template <class Evaluation1, class Evaluation2>
[[nodiscard]] typename ReturnEval_<Evaluation1, Evaluation2>::type
pow(const Evaluation1 &base, const Evaluation2 &exponent)
{
    using Result =
        typename ReturnEval_<Evaluation1, Evaluation2>::type;

    return MathToolbox<Result>::pow(base, exponent);
}

template <class Evaluation>
[[nodiscard]] bool isfinite(const Evaluation &input)
{
    return MathToolbox<Evaluation>::isfinite(input);
}

template <class Evaluation>
[[nodiscard]] bool isnan(const Evaluation &input)
{
    return MathToolbox<Evaluation>::isnan(input);
}

} // namespace MPMC
