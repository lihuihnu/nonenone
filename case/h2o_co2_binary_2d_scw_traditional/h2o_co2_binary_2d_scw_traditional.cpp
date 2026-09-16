/**
 * @file h2o_co2_binary_2d_scw_traditional.cpp
 * @brief 超临界水温压统一组成 H2O-CO2 传统模型算例入口。
 */
#define MPMC_H2O_CO2_NC10_TRADITIONAL_CASE_NAME "h2o_co2_binary_2d_scw_traditional"
#define MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA 28.0e6
#define MPMC_H2O_CO2_NC10_TEMPERATURE_K 653.15
#define MPMC_H2O_CO2_BINARY_INITIAL_H2O_MOLE_FRACTION 0.80
#include "../h2o_co2_nc10_2d_traditional/h2o_co2_nc10_2d_traditional.cpp"

static_assert(BenchmarkCommon::initialPressure > 22.064e6);
static_assert(BenchmarkCommon::temperature > 647.096);
