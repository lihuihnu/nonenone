/**
 * @file h2o_co2_binary_2d_conventional_new.cpp
 * @brief 常规温压统一组成 H2O-CO2 新型状态方程算例入口。
 */
#define MPMC_H2O_CO2_NC10_CASE_NAME "h2o_co2_binary_2d_conventional_new"
#define MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA 5.16e6
#define MPMC_H2O_CO2_NC10_TEMPERATURE_K 333.15
#define MPMC_H2O_CO2_NC10_FIXED_INITIAL_H2O_MOLE_FRACTION 0.80
#define MPMC_H2O_CO2_NC10_FIXED_INITIAL_CO2_MOLE_FRACTION 0.20
#include "../h2o_co2_nc10_2d_benchmark/h2o_co2_nc10_2d_benchmark.cpp"

static_assert(BenchmarkCommon::initialPressure == 5.16e6);
static_assert(BenchmarkCommon::temperature == 333.15);
