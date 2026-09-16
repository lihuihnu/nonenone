/**
 * @file h2o_co2_nc10_2d_scw_benchmark.cpp
 * @brief 653.15 K、28 MPa 下复用 New-PR、New-SW 与 New-CPA 的二维算例。
 */
#define MPMC_H2O_CO2_NC10_CASE_NAME "h2o_co2_nc10_2d_scw_benchmark"
#define MPMC_H2O_CO2_NC10_INITIAL_PRESSURE_PA 28.0e6
#define MPMC_H2O_CO2_NC10_TEMPERATURE_K 653.15
#define MPMC_H2O_CO2_NC10_FIXED_INITIAL_H2O_MOLE_FRACTION 0.20

#include "../h2o_co2_nc10_2d_benchmark/h2o_co2_nc10_2d_benchmark.cpp"

static_assert(BenchmarkCommon::initialPressure > 22.064e6);
static_assert(BenchmarkCommon::temperature > 647.096);
