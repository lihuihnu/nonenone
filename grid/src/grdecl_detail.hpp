/**
 * @file grdecl_detail.hpp
 * @brief GRDECL 读取实现的内部共享声明，不属于公共 API。
 */
#pragma once

#include <cpgrid/grdecl.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace MPMC::grdecl_detail
{

/** 递归展开 INCLUDE 后保存的原始 keyword token。 */
struct RawDeck final
{
    std::unordered_map<std::string, std::vector<std::string>> keywordValues;
    std::optional<GrdeclUnitSystem> units;
};

[[nodiscard]] std::vector<std::string> readTokens(const std::filesystem::path &file);

[[nodiscard]] RawDeck parseDeck(const std::filesystem::path &filename);

[[nodiscard]] std::vector<double> expandReal(
    const std::vector<std::string> &tokens,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize = std::nullopt);

[[nodiscard]] std::vector<int> expandInt(
    const std::vector<std::string> &tokens,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize = std::nullopt);

[[nodiscard]] const std::vector<std::string> &required(
    const RawDeck &deck,
    const std::string &keyword);

[[nodiscard]] std::optional<std::vector<double>> optionalReal(
    const RawDeck &deck,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize = std::nullopt);

[[nodiscard]] std::optional<GrdeclUnitSystem> keywordLengthUnit(
    const RawDeck &deck,
    const std::string &keyword);

[[nodiscard]] std::array<int, 3> dimensions(const RawDeck &deck);

[[nodiscard]] double geometryLengthFactor(GrdeclUnitSystem units) noexcept;

void convertGeometryToSi(
    std::vector<double> &coord,
    std::vector<double> &zcorn,
    GrdeclUnitSystem units);

void applyMapAxes(
    std::vector<double> &coord,
    const std::vector<double> &mapAxes,
    double originLengthFactor);

void convertPermeabilityToSi(std::vector<double> &values) noexcept;

/** active corner-point topology 构造结果。 */
struct ActiveGridBuild final
{
    std::vector<GrdeclCell> activeCells;
    std::vector<int> storageByCartesian;
};

[[nodiscard]] ActiveGridBuild buildActiveGrid(
    const std::vector<double> &coord,
    const std::vector<double> &zcorn,
    const std::vector<int> &actnum,
    int nx,
    int ny,
    int nz);

[[nodiscard]] std::size_t retainLargestConnectedComponent(
    std::vector<int> &actnum,
    int nx,
    int ny,
    int nz);

[[nodiscard]] std::vector<double> activeProperty(
    const std::optional<std::vector<double>> &raw,
    const std::vector<int> &storageByCartesian,
    std::size_t cartesianCount,
    const std::string &keyword);

} // namespace MPMC::grdecl_detail
