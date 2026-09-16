/**
 * @file grdecl_parser.cpp
 * @brief GRDECL token、INCLUDE、keyword 与数值展开实现。
 */
#include "grdecl_detail.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace MPMC::grdecl_detail
{
namespace
{

std::string upper(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

bool isStandaloneKeyword(const std::string &keyword)
{
    static const std::unordered_set<std::string> standalone{
        "RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION",
        "SUMMARY", "SCHEDULE", "ECHO", "NOECHO"};
    return standalone.count(keyword) != 0;
}

void parseFileRecursive(
    const std::filesystem::path &file,
    RawDeck &deck,
    std::vector<std::filesystem::path> &includeStack)
{
    if (includeStack.size() >= 17)
        throw std::runtime_error("GRDECL INCLUDE nesting exceeds 16 levels.");

    const auto normalized = std::filesystem::absolute(file).lexically_normal();
    if (std::find(includeStack.begin(), includeStack.end(), normalized) !=
        includeStack.end())
        throw std::runtime_error(
            "GRDECL INCLUDE cycle detected at " + normalized.string());

    includeStack.push_back(normalized);
    struct StackPop final
    {
        std::vector<std::filesystem::path> &stack;
        ~StackPop() { stack.pop_back(); }
    } stackPop{includeStack};

    const auto tokens = readTokens(normalized);
    std::size_t pos = 0;
    while (pos < tokens.size())
    {
        const std::string keyword = upper(tokens[pos++]);
        if (keyword == "/")
            continue;

        if (keyword == "METRIC" || keyword == "FIELD" || keyword == "LAB")
        {
            deck.units = keyword == "FIELD" ? GrdeclUnitSystem::Field
                       : keyword == "LAB" ? GrdeclUnitSystem::Lab
                                           : GrdeclUnitSystem::Metric;
            if (pos < tokens.size() && tokens[pos] == "/")
                ++pos;
            continue;
        }

        if (isStandaloneKeyword(keyword))
        {
            if (pos < tokens.size() && tokens[pos] == "/")
                ++pos;
            continue;
        }

        if (keyword == "INCLUDE")
        {
            if (pos >= tokens.size())
                throw std::runtime_error("INCLUDE is missing a filename in " + normalized.string());
            const std::filesystem::path included = normalized.parent_path() / tokens[pos++];
            while (pos < tokens.size() && tokens[pos] != "/")
                ++pos;
            if (pos == tokens.size())
                throw std::runtime_error(
                    "INCLUDE is missing its terminating '/' in " + normalized.string());
            ++pos;
            parseFileRecursive(included.lexically_normal(), deck, includeStack);
            continue;
        }

        std::vector<std::string> values;
        while (pos < tokens.size() && tokens[pos] != "/")
            values.push_back(tokens[pos++]);
        if (pos == tokens.size())
            throw std::runtime_error(
                keyword + " is missing its terminating '/' in " + normalized.string());
        ++pos;

        // 数值：后出现的定义覆盖前值，保持既有 deck override 语义。
        deck.keywordValues[keyword] = std::move(values);
    }
}

GrdeclUnitSystem namedLengthUnit(
    const std::string &value,
    const std::string &keyword)
{
    const std::string name = upper(value);
    if (name == "METRES" || name == "METERS" || name == "METRE" ||
        name == "METER" || name == "METRIC" || name == "M")
        return GrdeclUnitSystem::Metric;
    if (name == "FEET" || name == "FOOT" || name == "FT" || name == "FIELD")
        return GrdeclUnitSystem::Field;
    if (name == "CM" || name == "CENTIMETRES" || name == "CENTIMETERS" ||
        name == "CENTIMETRE" || name == "CENTIMETER" || name == "LAB")
        return GrdeclUnitSystem::Lab;
    throw std::runtime_error(keyword + " contains unsupported length unit '" + value + "'.");
}

} // namespace

RawDeck parseDeck(const std::filesystem::path &filename)
{
    RawDeck deck;
    std::vector<std::filesystem::path> includeStack;
    parseFileRecursive(filename, deck, includeStack);
    return deck;
}

std::vector<double> expandReal(
    const std::vector<std::string> &tokens,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize)
{
    std::vector<double> out;
    if (expectedSize)
        out.reserve(*expectedSize);

    for (const auto &token : tokens)
    {
        const auto star = token.find('*');
        std::size_t count = 1;
        std::string valueToken = token;
        if (star != std::string::npos)
        {
            const std::string countToken = token.substr(0, star);
            if (countToken.empty())
                throw std::runtime_error(keyword + " contains malformed repetition token: " + token);
            std::size_t parsed = 0;
            const auto parsedCount = std::stoull(countToken, &parsed);
            if (parsed != countToken.size() ||
                parsedCount > static_cast<unsigned long long>(
                    std::numeric_limits<std::size_t>::max()))
                throw std::runtime_error(
                    keyword + " contains invalid repetition count: " + token);
            count = static_cast<std::size_t>(parsedCount);
            valueToken = token.substr(star + 1);
            if (valueToken.empty())
                throw std::runtime_error(
                    keyword + " uses defaulted repetition 'n*', which is not supported by the direct GRDECL reader: " + token);
        }
        std::replace(valueToken.begin(), valueToken.end(), 'D', 'E');
        std::replace(valueToken.begin(), valueToken.end(), 'd', 'e');
        std::size_t parsed = 0;
        const double value = std::stod(valueToken, &parsed);
        if (parsed != valueToken.size() || !std::isfinite(value))
            throw std::runtime_error(
                keyword + " contains an invalid finite real value: " + token);

        if (expectedSize)
        {
            if (out.size() > *expectedSize || count > *expectedSize - out.size())
                throw std::runtime_error(
                    keyword + " size exceeds its expected value count.");
        }
        if (count > out.max_size() - out.size())
            throw std::length_error(
                keyword + " repetition count exceeds vector capacity.");
        out.insert(out.end(), count, value);
    }

    if (expectedSize && out.size() != *expectedSize)
        throw std::runtime_error(
            keyword + " size does not match its expected value count.");
    return out;
}

std::vector<int> expandInt(
    const std::vector<std::string> &tokens,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize)
{
    const auto real = expandReal(tokens, keyword, expectedSize);
    std::vector<int> out;
    out.reserve(real.size());
    for (double value : real)
    {
        const double rounded = std::round(value);
        if (std::abs(value - rounded) > 1.0e-9 ||
            rounded < static_cast<double>(std::numeric_limits<int>::min()) ||
            rounded > static_cast<double>(std::numeric_limits<int>::max()))
            throw std::runtime_error(keyword + " requires integer values.");
        out.push_back(static_cast<int>(rounded));
    }
    return out;
}

const std::vector<std::string> &required(
    const RawDeck &deck,
    const std::string &keyword)
{
    const auto it = deck.keywordValues.find(keyword);
    if (it == deck.keywordValues.end())
        throw std::runtime_error("GRDECL is missing required keyword " + keyword + '.');
    return it->second;
}

std::optional<std::vector<double>> optionalReal(
    const RawDeck &deck,
    const std::string &keyword,
    std::optional<std::size_t> expectedSize)
{
    const auto it = deck.keywordValues.find(keyword);
    if (it == deck.keywordValues.end())
        return std::nullopt;
    return expandReal(it->second, keyword, expectedSize);
}

std::optional<GrdeclUnitSystem> keywordLengthUnit(
    const RawDeck &deck,
    const std::string &keyword)
{
    const auto it = deck.keywordValues.find(keyword);
    if (it == deck.keywordValues.end())
        return std::nullopt;
    if (it->second.empty())
        throw std::runtime_error(keyword + " requires a length unit.");
    return namedLengthUnit(it->second.front(), keyword);
}

std::array<int, 3> dimensions(const RawDeck &deck)
{
    auto it = deck.keywordValues.find("SPECGRID");
    if (it == deck.keywordValues.end())
        it = deck.keywordValues.find("DIMENS");
    if (it == deck.keywordValues.end())
        throw std::runtime_error("GRDECL requires SPECGRID or DIMENS.");
    if (it->second.size() < 3)
        throw std::runtime_error("SPECGRID/DIMENS must contain nx ny nz.");

    const std::vector<std::string> firstThree{it->second[0], it->second[1], it->second[2]};
    const auto vals = expandInt(firstThree, "SPECGRID/DIMENS");
    if (vals.size() != 3 || vals[0] <= 0 || vals[1] <= 0 || vals[2] <= 0)
        throw std::runtime_error("Invalid GRDECL grid dimensions.");
    return {vals[0], vals[1], vals[2]};
}

} // namespace MPMC::grdecl_detail
