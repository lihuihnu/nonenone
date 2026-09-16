/**
 * @file writer.hpp
 * @brief 网格无关的 CSV/文本行写入接口。
 */
#pragma once

#include <output/core/format.hpp>
#include <output/core/text_file.hpp>
#include <output/core/types.hpp>
#include <output/metrics/component_totals.hpp>

#include <array>
#include <filesystem>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 标准 MPMC CSV 文件的 rank-local 写入器。
 *
 * MPI 运行时仅指定 writer rank 真正打开文件，其余 rank 可调用同一 API，
 * 从而避免上层 collective 输出代码散布 rank 分支。文件名和数值格式集中在
 * output 模块，网格/PETSc 适配层只负责采集数据。
 */
template <class Indices>
class OutputWriter final
{
public:
    explicit OutputWriter(
        std::filesystem::path outputDirectory,
        bool writerEnabled = true)
        : outputDirectory_(std::move(outputDirectory)),
          writerEnabled_(writerEnabled)
    {
        solution_.open(outputDirectory_ / "solution_samples.csv", writerEnabled_);
        phaseState_.open(outputDirectory_ / "phase_state_samples.csv", writerEnabled_);
        well_.open(outputDirectory_ / "well_phase_rates.csv", writerEnabled_);
        if (writerEnabled_)
        {
            for (int phase = 0; phase < Indices::numPhases; ++phase)
            {
                if (phase != 0)
                    well_.stream() << ',';
                well_.stream() << "phase_" << phase << "_surface_rate_m3_s";
            }
            well_.stream() << '\n';
        }
        openMassFiles();
    }

    OutputWriter(const OutputWriter &) = delete;
    OutputWriter &operator=(const OutputWriter &) = delete;
    OutputWriter(OutputWriter &&) = delete;
    OutputWriter &operator=(OutputWriter &&) = delete;

    /** @brief 追加一个采样时刻的主变量 CSV 记录。 */
    void writeSolutionSamples(const std::vector<double> &values)
    {
        if (!writerEnabled_)
            return;
        ensureSampleHeader_(solution_.stream(), solutionHeaderWritten_, values.size());
        writeCsvValues(solution_.stream(), values);
    }

    /** @brief 追加一个采样时刻的相态变量记录。 */
    void writePhaseStateSamples(const std::vector<double> &values)
    {
        if (!writerEnabled_)
            return;
        ensureSampleHeader_(phaseState_.stream(), phaseStateHeaderWritten_, values.size());
        writeCsvValues(phaseState_.stream(), values);
    }

    /** @brief 按配置的相顺序追加一口井的各相流量。 */
    void writeWellPhaseRates(const std::array<double, Indices::numPhases> &rates)
    {
        if (writerEnabled_)
            MPMC::writeWellPhaseRatesCsv(well_.stream(), rates);
    }

    /** @brief 为一个启用的库存序列追加各组分质量。 */
    void writeMassSamples(MassSeries series, const std::vector<double> &values)
    {
        requireMassSeries(series);
        if (!writerEnabled_)
            return;
        const auto index = massSeriesIndex(series);
        ensureSampleHeader_(massFields_[index].stream(), massFieldHeaderWritten_[index], values.size());
        writeCsvValues(massFields_[index].stream(), values);
    }

    /** @brief 为一个质量库存序列追加全局组分总量。 */
    void writeMassTotals(
        MassSeries series,
        double time,
        const ComponentTotals<static_cast<std::size_t>(Indices::numComponents)> &totals)
    {
        requireMassSeries(series);
        if (writerEnabled_)
            writeMassTotalLine(
                massTotals_[massSeriesIndex(series)].stream(), series, time, totals);
    }

private:
    static constexpr std::size_t componentCount =
        static_cast<std::size_t>(Indices::numComponents);
    static constexpr std::size_t massSeriesCount =
        static_cast<std::size_t>(MassSeries::Count);

    static void requireMassSeries(MassSeries series)
    {
        if (!massSeriesEnabled<Indices>(series))
            throw std::logic_error(
                "Requested mass output is disabled by the model configuration.");
    }

    static void ensureSampleHeader_(
        std::ostream &stream,
        bool &written,
        std::size_t valueCount)
    {
        if (written)
            return;
        writeIndexedCsvHeader(stream, valueCount);
        written = true;
    }

    void openMassFiles()
    {
        for (std::size_t index = 0; index < massSeriesCount; ++index)
        {
            const auto series = static_cast<MassSeries>(index);
            const auto &info = massSeriesInfo(series);
            const bool enabled = writerEnabled_ && massSeriesEnabled<Indices>(series);

            massFields_[index].open(outputDirectory_ / info.fieldFilename, enabled);
            massTotals_[index].open(outputDirectory_ / info.totalFilename, enabled);

            if (enabled)
            {
                writeMassHeader<componentCount>(massTotals_[index].stream(), series);
                massTotals_[index].stream().flush();
            }
        }
    }

    std::filesystem::path outputDirectory_;
    bool writerEnabled_{true};
    OutputTextFile solution_;
    OutputTextFile phaseState_;
    OutputTextFile well_;
    bool solutionHeaderWritten_{false};
    bool phaseStateHeaderWritten_{false};
    std::array<bool, massSeriesCount> massFieldHeaderWritten_{};
    std::array<OutputTextFile, massSeriesCount> massFields_{};
    std::array<OutputTextFile, massSeriesCount> massTotals_{};
};

} // namespace MPMC
