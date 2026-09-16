/**
 * @file common_test.cpp
 * @brief 单元测试：验证 `common` 的核心语义、边界条件和回归行为。
 */
#include <common/bilinear_interpolation.hpp>
#include <common/console.hpp>
#include <common/interval_search.hpp>
#include <common/math.hpp>
#include <common/piecewise_linear_interpolation.hpp>
#include <common/small_vector.hpp>
#include <common/type_traits.hpp>
#include <common/units.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Dual
{
    double value{0.0};
    double derivative{0.0};

    Dual() = default;
    Dual(double v) : value(v) {}
    Dual(double v, double d)
        : value(v), derivative(d)
    {
    }
};

Dual operator+(const Dual &a, const Dual &b)
{
    return {a.value + b.value,
            a.derivative + b.derivative};
}

Dual operator-(const Dual &a, const Dual &b)
{
    return {a.value - b.value,
            a.derivative - b.derivative};
}

Dual operator*(const Dual &a, const Dual &b)
{
    return {
        a.value * b.value,
        a.derivative * b.value +
            a.value * b.derivative};
}

Dual operator/(const Dual &a, const Dual &b)
{
    return {
        a.value / b.value,
        (a.derivative * b.value -
         a.value * b.derivative) /
            (b.value * b.value)};
}

Dual operator+(double a, const Dual &b)
{
    return Dual(a) + b;
}

Dual operator-(const Dual &a, double b)
{
    return a - Dual(b);
}

Dual operator-(double a, const Dual &b)
{
    return Dual(a) - b;
}

Dual operator*(const Dual &a, double b)
{
    return a * Dual(b);
}

Dual operator/(const Dual &a, double b)
{
    return a / Dual(b);
}

bool near(
    double actual,
    double expected,
    double tolerance = 1.0e-12)
{
    return std::abs(actual - expected) <=
           tolerance;
}

void require(
    bool condition,
    const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

} // namespace

namespace MPMC
{

template <>
struct MathToolbox<Dual>
{
    using Scalar = double;
    using ValueType = Dual;

    static double scalarValue(
        const Dual &input) noexcept
    {
        return input.value;
    }
};

} // namespace MPMC

int main()
{
    using MPMC::BoundaryPolicy;
    using MPMC::BilinearInterpolation;
    using MPMC::PiecewiseLinearInterpolation;

    // --------------------------------------------------------------
    // interval search: ascending + descending + boundaries
    // --------------------------------------------------------------
    {
        const std::vector<double> ascending{
            0.0, 1.0, 2.0};

        require(
            MPMC::findInterval(
                ascending,
                -1.0) == 0,
            "ascending left extrapolation interval failed");

        require(
            MPMC::findInterval(
                ascending,
                2.0) == 1,
            "ascending right boundary interval failed");

        const std::vector<double> descending{
            2.0, 1.0, 0.0};

        require(
            MPMC::findInterval(
                descending,
                0.0) == 1,
            "descending minimum boundary must return last valid segment");

        require(
            MPMC::findInterval(
                descending,
                -1.0) == 1,
            "descending right extrapolation interval failed");
    }

    // --------------------------------------------------------------
    // piecewise linear interpolation + AD derivative
    // --------------------------------------------------------------
    {
        const std::vector<double> x{
            0.0, 1.0, 2.0};

        const std::vector<double> y{
            0.0, 2.0, 4.0};

        const Dual query{
            0.25,
            1.0};

        const Dual result =
            PiecewiseLinearInterpolation::evaluate(
                x,
                y,
                query);

        require(
            near(result.value, 0.5),
            "piecewise interpolation value failed");

        require(
            near(result.derivative, 2.0),
            "piecewise interpolation derivative failed");

        const Dual clamped =
            PiecewiseLinearInterpolation::evaluate(
                x,
                y,
                Dual(-1.0, 1.0),
                BoundaryPolicy::Clamp);

        require(
            near(clamped.value, 0.0) &&
                near(clamped.derivative, 0.0),
            "piecewise clamp failed");

        const Dual extrapolated =
            PiecewiseLinearInterpolation::evaluate(
                x,
                y,
                Dual(-1.0, 1.0),
                BoundaryPolicy::LinearExtrapolation);

        require(
            near(extrapolated.value, -2.0) &&
                near(extrapolated.derivative, 2.0),
            "piecewise extrapolation failed");

        const std::vector<double> xd{
            2.0, 1.0, 0.0};

        const std::vector<double> yd{
            4.0, 2.0, 0.0};

        const Dual descending =
            PiecewiseLinearInterpolation::evaluate(
                xd,
                yd,
                query);

        require(
            near(descending.value, 0.5) &&
                near(descending.derivative, 2.0),
            "descending piecewise interpolation failed");
    }

    // --------------------------------------------------------------
    // bilinear interpolation on plane f = 1 + 2x + 3y
    // --------------------------------------------------------------
    {
        const std::vector<double> x{
            0.0, 1.0};

        const std::vector<double> y{
            0.0, 1.0};

        const std::vector<std::vector<double>> values{
            {1.0, 3.0},
            {4.0, 6.0}};

        const Dual resultX =
            BilinearInterpolation::evaluate(
                x,
                y,
                values,
                Dual(0.2, 1.0),
                Dual(0.4, 0.0));

        require(
            near(resultX.value, 2.6) &&
                near(resultX.derivative, 2.0),
            "bilinear d/dx failed");

        const Dual resultY =
            BilinearInterpolation::evaluate(
                x,
                y,
                values,
                Dual(0.2, 0.0),
                Dual(0.4, 1.0));

        require(
            near(resultY.value, 2.6) &&
                near(resultY.derivative, 3.0),
            "bilinear d/dy failed");

        const Dual extrapolated =
            BilinearInterpolation::evaluate(
                x,
                y,
                values,
                Dual(-1.0, 1.0),
                Dual(0.0, 0.0));

        require(
            near(extrapolated.value, -1.0) &&
                near(extrapolated.derivative, 2.0),
            "bilinear default linear extrapolation failed");
    }

    // --------------------------------------------------------------
    // SmallVector inline / heap / copy / move / self-move
    // --------------------------------------------------------------
    {
        MPMC::SmallVector<double, 4> small(
            3,
            2.0);

        require(
            small.size() == 3 &&
                near(small[2], 2.0),
            "SmallVector inline storage failed");

        MPMC::SmallVector<double, 4> large(
            8,
            3.0);

        large[7] = 9.0;

        auto copied = large;

        require(
            copied.size() == 8 &&
                near(copied[7], 9.0),
            "SmallVector copy failed");

        auto moved =
            std::move(large);

        require(
            moved.size() == 8 &&
                near(moved[7], 9.0) &&
                large.size() == 0,
            "SmallVector move failed");

        auto *self = &moved;
        moved = std::move(*self);

        require(
            moved.size() == 8 &&
                near(moved[7], 9.0),
            "SmallVector self-move must preserve the object");
    }

    // --------------------------------------------------------------
    // traits + units + numerical comparison
    // --------------------------------------------------------------
    {
        static_assert(
            MPMC::is_specialization_of_v<
                std::vector<int>,
                std::vector>);

        static_assert(
            !MPMC::is_specialization_of_v<
                int,
                std::vector>);

        require(
            near(
                MPMC::units::psi,
                6894.757293168361,
                1.0e-9),
            "psi conversion failed");

        require(
            MPMC::MathToolbox<double>::isSame(
                1.0e8,
                1.0e8 + 0.01,
                1.0e-9),
            "relative isSame failed");

        require(
            !MPMC::MathToolbox<double>::isSame(
                1.0,
                -1.0,
                1.0e-6),
            "opposite-sign isSame scaling failed");
    }

    // --------------------------------------------------------------
    // common console/report formatting
    // --------------------------------------------------------------
    {
        require(
            MPMC::consoleRule().size() == MPMC::consoleWidth,
            "console rule width failed");
        require(
            MPMC::consoleCentered("GRID").size() == MPMC::consoleWidth,
            "centered console title width failed");
        const std::string row = MPMC::consoleKeyValue("Porosity", "0.25");
        require(
            row.find("Porosity") != std::string::npos &&
                row.find("0.25") != std::string::npos,
            "console key/value formatting failed");
        MPMC::ConsoleSection section("TEST SECTION");
        section.row("Status", "SUCCESS").separator().line("Readable output");
        const std::string report = section.str();
        require(
            report.find("TEST SECTION") != std::string::npos &&
                report.find("Status") != std::string::npos &&
                report.find("SUCCESS") != std::string::npos &&
                report.find("Readable output") != std::string::npos,
            "console section formatting failed");
    }

    // --------------------------------------------------------------
    // invalid non-monotonic table must fail deterministically
    // --------------------------------------------------------------
    {
        bool threw = false;

        try
        {
            const std::vector<double> x{
                0.0, 2.0, 1.0};

            const std::vector<double> y{
                0.0, 2.0, 1.0};

            (void)PiecewiseLinearInterpolation::evaluate(
                x,
                y,
                1.0);
        }
        catch (const std::invalid_argument &)
        {
            threw = true;
        }

        require(
            threw,
            "non-monotonic interpolation table must throw");
    }

    std::cout
        << "Common core validation: ALL PASS\n";

    return 0;
}
