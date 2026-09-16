#pragma once

#include "../scw_kerogen_common/benchmark_2d_common.hpp"
#include "../scw_kerogen_squalane_1d/case_config.hpp"

#include <natural/thermo/soreide_whitson.hpp>

#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 60x20x1 水–角鲨烷 PR/SW/CPA 二维对比配置。
 */
namespace ScwKerogenSqualane2D
{

inline constexpr char name[] = "scw_kerogen_squalane_2d";
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
    static double cubicBinaryInteractionCoefficient(
        int i, int j, double temperatureK)
    {
        if (i == j)
            return 0.0;
        return ScwKerogenCalibration::prSqualaneKij(temperatureK);
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
        // Explicitly an extrapolation: the original hydrocarbon correlation
        // was fitted only through nC4, not C30 squalane.
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
    inline static constexpr const char *name = ScwKerogenSqualane2D::name;
    using Model = ScwKerogenSqualane2D::Model;
    using Grid = ScwKerogenSqualane2D::Grid;
    using Fluid = ScwKerogenSqualane2D::PrFluid;
    using InitialState = ScwKerogenSqualane2D::InitialState;
    using Dissolution = ScwKerogenSqualane2D::Dissolution;
    using Land = ScwKerogenSqualane2D::Land;
    using Adsorption = ScwKerogenSqualane2D::Adsorption;
    using Numerics = ScwKerogenSqualane2D::Numerics;
    using Time = ScwKerogenSqualane2D::Time;
    using Output = ScwKerogenSqualane2D::Output;
    using Rock = ScwKerogenSqualane2D::Rock;
};

} // namespace ScwKerogenSqualane2D
