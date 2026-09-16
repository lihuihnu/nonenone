/**
 * @file console.hpp
 * @brief 并行运行时的统一控制台输出辅助。
 */
#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

namespace MPMC
{

inline constexpr std::size_t consoleWidth = 78;
inline constexpr std::size_t consoleLabelWidth = 28;

inline std::string consoleRule(char fill = '=')
{
    return std::string(consoleWidth, fill);
}

inline std::string consoleCentered(std::string_view text, char fill = '=')
{
    const std::string title = " " + std::string(text) + " ";
    if (title.size() >= consoleWidth)
        return title;
    const std::size_t remaining = consoleWidth - title.size();
    const std::size_t left = remaining / 2;
    return std::string(left, fill) + title + std::string(remaining - left, fill);
}

inline std::string consoleKeyValue(
    std::string_view label,
    std::string_view value,
    std::string_view unit = {})
{
    std::ostringstream out;
    out << "  " << std::left << std::setw(static_cast<int>(consoleLabelWidth)) << label
        << " : " << value;
    if (!unit.empty())
        out << ' ' << unit;
    return out.str();
}

inline std::string consoleNumber(double value, int precision = 6)
{
    std::ostringstream out;
    const double magnitude = value < 0.0 ? -value : value;
    if ((magnitude != 0.0 && magnitude < 1.0e-3) || magnitude >= 1.0e6)
        out << std::scientific;
    else
        out << std::fixed;
    out << std::setprecision(precision) << value;
    return out.str();
}

inline std::string consoleScientific(double value, int precision = 6)
{
    std::ostringstream out;
    out << std::scientific << std::setprecision(precision) << value;
    return out.str();
}

inline const char *consoleOnOff(bool value) noexcept
{
    return value ? "ON" : "OFF";
}

/** Build one visually consistent console section. */
class ConsoleSection final
{
public:
    explicit ConsoleSection(std::string title)
    {
        stream_ << '\n' << consoleCentered(title) << '\n';
    }

    ConsoleSection &row(
        std::string_view label,
        std::string_view value,
        std::string_view unit = {})
    {
        stream_ << consoleKeyValue(label, value, unit) << '\n';
        return *this;
    }

    template <class T>
    ConsoleSection &rowValue(
        std::string_view label,
        const T &value,
        std::string_view unit = {})
    {
        std::ostringstream text;
        text << value;
        return row(label, text.str(), unit);
    }

    ConsoleSection &line(std::string_view value)
    {
        stream_ << "  " << value << '\n';
        return *this;
    }

    ConsoleSection &separator(char fill = '-')
    {
        stream_ << consoleRule(fill) << '\n';
        return *this;
    }

    [[nodiscard]] std::string str(bool close = true) const
    {
        std::string value = stream_.str();
        if (close)
            value += consoleRule() + "\n";
        return value;
    }

private:
    std::ostringstream stream_;
};

} // namespace MPMC
