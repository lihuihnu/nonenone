/**
 * @file csv_reader.hpp
 * @brief 网格属性 CSV 文件的轻量读取工具。
 */
#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace MPMC
{

namespace detail
{

inline std::string trimCsvField(const std::string &value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

template <class Number, class Parser>
Number parseCsvNumber(const std::string &value, Parser parser)
{
    const std::string field = trimCsvField(value);
    if (field.empty())
        throw std::invalid_argument("empty numeric CSV field");

    std::size_t parsed = 0;
    const Number result = parser(field, &parsed);
    if (parsed != field.size())
        throw std::invalid_argument("numeric CSV field contains trailing characters");
    return result;
}

template <class T>
T parseCsvValue(const std::string &value)
{
    if constexpr (std::is_same_v<T, int>)
    {
        return parseCsvNumber<int>(
            value,
            [](const std::string &field, std::size_t *parsed) {
                return std::stoi(field, parsed);
            });
    }
    else if constexpr (std::is_same_v<T, long long>)
    {
        return parseCsvNumber<long long>(
            value,
            [](const std::string &field, std::size_t *parsed) {
                return std::stoll(field, parsed);
            });
    }
    else if constexpr (std::is_same_v<T, double>)
    {
        const double result = parseCsvNumber<double>(
            value,
            [](const std::string &field, std::size_t *parsed) {
                return std::stod(field, parsed);
            });
        if (!std::isfinite(result))
            throw std::invalid_argument("numeric CSV field must be finite");
        return result;
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
        return value;
    }
    else
    {
        static_assert(
            !sizeof(T),
            "Unsupported CSV value type.");
    }
}

} // namespace detail

/**
 * @brief 读取普通 CSV，不自动修改数值。
 */
template <class T>
std::vector<std::vector<T>>
readCsv(const std::string &filename)
{
    std::ifstream file(filename);

    if (!file)
    {
        throw std::runtime_error(
            "Cannot open CSV file: " + filename);
    }

    std::vector<std::vector<T>> data;
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(file, line))
    {
        ++lineNumber;

        if (line.empty())
        {
            continue;
        }

        std::vector<T> row;
        std::stringstream stream(line);
        std::string value;

        if (!line.empty() && line.back() == ',')
        {
            throw std::runtime_error(
                "Empty trailing CSV field at " + filename + ":" +
                std::to_string(lineNumber));
        }

        while (std::getline(stream, value, ','))
        {
            try
            {
                row.push_back(
                    detail::parseCsvValue<T>(value));
            }
            catch (const std::exception &)
            {
                throw std::runtime_error(
                    "Invalid CSV value at " +
                    filename +
                    ":" +
                    std::to_string(lineNumber) +
                    " -> " +
                    value);
            }
        }

        data.push_back(std::move(row));
    }

    return data;
}

/**
 * @brief 读取 MATLAB/MRST 1 基索引 CSV，并转为 C++ 0 基索引。
 *
 * MATLAB 中的 0（常用于 boundary neighbor）会被转换成 -1。
 */
template <class T = int>
std::vector<std::vector<T>>
readMatlabIndexCsv(
    const std::string &filename)
{
    static_assert(
        std::is_integral_v<T> && std::is_signed_v<T>,
        "MATLAB index CSV values require a signed integral type.");

    auto data = readCsv<T>(filename);

    for (auto &row : data)
    {
        for (T &value : row)
        {
            --value;
        }
    }

    return data;
}

} // namespace MPMC
