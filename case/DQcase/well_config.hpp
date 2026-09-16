#pragma once

#include "case_config.hpp"

#include <case/well_config_schema.hpp>

#include <initializer_list>
#include <vector>

/**
 * @file well_config.hpp
 * @brief DQ 网格三相多组分超算算例的井位置、完井和控制配置。
 *
 * 完井 cell 统一使用原始网格数据文件的 0 基行号（input cell index）。
 * 该序号与 `G/cells/indexMap/data.csv`、`rock/perm/data.csv` 和
 * `rock/poro/data.csv` 的行顺序一致，且不会因 METIS 分区改变。
 * DQcase 适配层会在求解前自动转换为 CpGrid current id。
 */
namespace WellConfig
{

using Type = MPMC::cases::well_config::Type;
using Control = MPMC::cases::well_config::Control;
using InjectionPhase = MPMC::cases::well_config::InjectionPhase;
using Perforation = MPMC::cases::well_config::Perforation;

struct WellDefaults
{
    static constexpr double radius = 0.0762;
};

using WellDefinition =
    MPMC::cases::well_config::PerforatedWellDefinition<CaseConfig::Config, WellDefaults>;

inline double perDay(double value)
{
    return value / CaseConfig::secondsPerDay;
}

inline WellDefinition injector(
    int id,
    const char *name,
    std::initializer_list<Perforation> perforations)
{
    WellDefinition w;
    w.id = id;
    w.name = name;
    w.type = Type::Injector;
    w.control = Control::TotalRate;
    w.target = perDay(600.0);
    w.initialBhp = CaseConfig::InitialState::pressure;
    w.perforations.assign(perforations.begin(), perforations.end());
    w.injectionPhase = InjectionPhase::Gas;
    w.injectedComponent = 0; // pure CO2
    w.maximumBhp = 1000.0 * CaseConfig::bar;
    w.openTime = 0.0;
    w.closeTime = 2555.0 * CaseConfig::secondsPerDay;
    return w;
}

inline WellDefinition producer(
    int id,
    const char *name,
    Control control,
    double target,
    std::initializer_list<Perforation> perforations)
{
    WellDefinition w;
    w.id = id;
    w.name = name;
    w.type = Type::Producer;
    w.control = control;
    w.target = target;
    w.initialBhp = CaseConfig::InitialState::pressure;
    w.perforations.assign(perforations.begin(), perforations.end());
    return w;
}

/**
 * @brief 当前 DQ 对齐基线井表。
 *
 * 所有完井单元号均为原始输入序号，不是分区后的 current id。
 * 生产率控制目标统一转为 m3/s；BHP 使用 Pa。
 */
inline const std::vector<WellDefinition> wells = [] {
    std::vector<WellDefinition> v{
        injector(0, "OY64-511", {{269, 4.9351847789466236e-15}}),
        injector(1, "OY66-50",  {{1035, 7.4640453743079468e-15},
                                  {2771, 1.9903872481047633e-15},
                                  {4057, 7.4749444587259123e-15}}),
        injector(2, "OY66-51",  {{914, 8.7528790009105677e-15},
                                  {2640, 2.1475630087345985e-15}}),
        injector(3, "OY66-52",  {{752, 1.0916826092691624e-14},
                                  {2469, 4.0005537084195953e-15}}),
        injector(4, "OY66-53",  {{674, 7.6298967604894955e-15},
                                  {2388, 4.9211250910250061e-15}}),
        injector(5, "OY66-54",  {{519, 1.9296985553810100e-15},
                                  {2214, 7.3641149972464782e-15}}),
        injector(6, "OY68-51",  {{1754, 6.1199939022610706e-15},
                                  {3500, 5.7866344861430104e-15},
                                  {4498, 1.4294172219340555e-14}}),
        injector(7, "OY68-52",  {{1642, 6.2959825026003330e-15},
                                  {3388, 5.0544736441596166e-15}}),
        injector(8, "OY68-53",  {{1467, 1.0667917598898876e-15},
                                  {3215, 9.8349560471826005e-15}}),

        producer(9,  "11Y65-50", Control::OilRate, perDay(0.8206666708),
                 {{602, 8.0968805928442583e-15}, {2313, 2.5783136236893869e-15}}),
        producer(10, "11Y65-51", Control::OilRate, perDay(0.2343333364),
                 {{526, 9.9453229902339043e-15}, {2231, 2.4842012086524136e-16}}),
        producer(11, "11Y65-52", Control::OilRate, perDay(0.1953333318),
                 {{408, 1.3424000781333505e-14}, {2066, 1.5614355719304994e-16}}),
        producer(12, "11Y65-53", Control::OilRate, perDay(0.1953333318),
                 {{335, 8.5929424684441359e-15}, {1978, 5.6226763139644488e-15}}),
        producer(13, "11Y67-50", Control::OilRate, perDay(0.1953333318),
                 {{1485, 4.9408498224653527e-15}, {3232, 7.5694107161995126e-15},
                  {4282, 1.0868158687694259e-14}}),
        producer(14, "11Y67-51", Control::OilRate, perDay(0.2343333364),
                 {{1358, 2.9726903332036368e-15}, {3105, 2.3225320733842404e-15},
                  {4192, 5.1799184515639604e-15}}),
        producer(15, "11Y67-52", Control::OilRate, perDay(0.2343333364),
                 {{1140, 9.7188272297666911e-15}}),
        producer(16, "11Y67-53", Control::OilRate, perDay(0.2343333364),
                 {{1065, 8.1886368250394252e-15}, {2801, 1.1334098830169261e-14}}),
        producer(17, "11Y67-54", Control::Bhp, 2.0 * CaseConfig::bar,
                 {{943, 3.0912955853569950e-15}})
    };

    v[9].minimumBhp = 1.01 * CaseConfig::bar;
    v[15].maximumWaterRate = 0.0;
    v[17].initialBhp = 2.0 * CaseConfig::bar;
    return v;
}();

} // namespace WellConfig
