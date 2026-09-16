#pragma once

/**
 * @file calibrated_binary_properties.hpp
 * @brief 超临界水算例可审计的 H2O–烃相平衡参数锚点。
 *
 * Values in this header deliberately retain their evidence status.  The
 * squalane parameters are a direct fit to the local Stevenson tie-line table.
 * The nC10 relation is a temperature regression of published PR kij values
 * below the target window and is therefore an extrapolation at 653.2 K.
 * The nC4 value is the published Soreide-Whitson non-aqueous prior, not a new
 * regression to raw target-window tie-line data.
 */
namespace ScwKerogenCalibration
{

inline constexpr double referenceTemperatureK = 653.2;

struct SqualaneTargetState
{
    static constexpr double temperatureK = 653.2;
    static constexpr double pressurePa = 27.74e6;
    static constexpr double waterRichH2OMoleFraction = 0.9966;
    static constexpr double hydrocarbonRichH2OMoleFraction = 0.9110;
    static constexpr double hydrocarbonRichEndpointUncertainty = 0.0020;
    // The common PR/SW/CPA comparison must start from the same phase count and
    // role.  CPA predicts two phases closer to the experimental PR-calibrated
    // endpoint, so the model-comparison baseline uses a conservative point on
    // the hydrocarbon-side one-phase branch.  The experimental endpoint remains
    // the explicit upper bound, rather than an artificial initial water phase.
    static constexpr double initialOverallH2OMoleFraction = 0.6000;
    static_assert(
        initialOverallH2OMoleFraction < hydrocarbonRichH2OMoleFraction,
        "The experimental-endpoint baseline must start outside the two-phase interval.");
};

/** Direct classical-PR composition fit to the selected Stevenson tie lines. */
[[nodiscard]] inline double prSqualaneKij(double temperatureK)
{
    constexpr double kReference = 0.053233659543451335;
    constexpr double inverseTemperatureSlopeK = 1005.045450882816;
    return kReference + inverseTemperatureSlopeK
        * (1.0 / temperatureK - 1.0 / referenceTemperatureK);
}

/** Published SW non-aqueous nC4 prior; direct target-window regression pending. */
[[nodiscard]] inline double nonAqueousNc4Kij(double /*temperatureK*/)
{
    return 0.5091;
}

/**
 * PR H2O-nC10 regression of Teratani et al. kij(lij=0) at 573.2, 593.2,
 * and 613.2 K: kij = -0.164670020 + 278.594624/T.
 * Evaluation above 613.2 K is explicitly an extrapolation.
 */
[[nodiscard]] inline double extrapolatedPrNc10Kij(double temperatureK)
{
    return -0.1646700200657025 + 278.594624154873 / temperatureK;
}

} // namespace ScwKerogenCalibration
