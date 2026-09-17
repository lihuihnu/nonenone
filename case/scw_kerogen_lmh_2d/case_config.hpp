#pragma once

#include "../scw_kerogen_common/benchmark_2d_common.hpp"
#include "../scw_kerogen_lmh_1d/case_config.hpp"

#include <natural/thermo/soreide_whitson.hpp>

#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 60x20x1 水–BSB轻中重拟组分 PR/SW/CPA 二维对比配置。
 */
namespace ScwKerogenLmh2D
{

inline constexpr char name[] = "scw_kerogen_lmh_2d";
using Model = CaseConfig::Model;
using Grid = ScwKerogen2D::Grid;
using Rock = ScwKerogen2D::Rock;
using Numerics = ScwKerogen2D::Numerics;
using Time = ScwKerogen2D::Time;
using Output = ScwKerogen2D::Output;
using InitialState = CaseConfig::InitialState;
using Dissolution = CaseConfig::Dissolution;
using Land = CaseConfig::Land;
using Adsorption = CaseConfig::Adsorption;

struct CommonFluid : CaseConfig::Fluid
{
    static double swTemperatureShape(int component, double temperatureK)
    {
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperatureK,
            criticalTemperature[static_cast<std::size_t>(component)],
            acentricFactor[static_cast<std::size_t>(component)], 0.0);
    }

    static double cubicBinaryInteractionCoefficient(
        int i, int j, double temperatureK)
    {
        if (i == j)
            return 0.0;
        const int other = i == waterComponent ? j
            : (j == waterComponent ? i : -1);
        if (other < 0)
            return 0.0;

        (void)temperatureK;
        return BsbReference::waterHydrocarbonKij;
    }
};

struct PrFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
};

struct SwFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::SoreideWhitson;
    static constexpr double soreideWhitsonSalinityMolality = 0.0;
    static double soreideWhitsonAqueousWaterBip(
        int component, double temperatureK, double salinityMolality)
    {
        if (component == waterComponent)
            return 0.0;
        return MPMC::SoreideWhitsonCorrelations::hydrocarbonAqueousBip(
            temperatureK,
            criticalTemperature[static_cast<std::size_t>(component)],
            acentricFactor[static_cast<std::size_t>(component)],
            salinityMolality);
    }
};

struct CpaFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::CubicPlusAssociation;
    static constexpr auto cpaCubicPhysicalTerm =
        MPMC::CpaCubicPhysicalTerm::PengRobinson;
};

struct PrFactoryConfig { using Fluid = PrFluid; };
struct SwFactoryConfig { using Fluid = SwFluid; };
struct CpaFactoryConfig { using Fluid = CpaFluid; };

struct Config final
{
    inline static constexpr const char *name = ScwKerogenLmh2D::name;
    using Model = ScwKerogenLmh2D::Model;
    using Grid = ScwKerogenLmh2D::Grid;
    using Fluid = ScwKerogenLmh2D::PrFluid;
    using InitialState = ScwKerogenLmh2D::InitialState;
    using Dissolution = ScwKerogenLmh2D::Dissolution;
    using Land = ScwKerogenLmh2D::Land;
    using Adsorption = ScwKerogenLmh2D::Adsorption;
    using Numerics = ScwKerogenLmh2D::Numerics;
    using Time = ScwKerogenLmh2D::Time;
    using Output = ScwKerogenLmh2D::Output;
    using Rock = ScwKerogenLmh2D::Rock;
};

} // namespace ScwKerogenLmh2D
