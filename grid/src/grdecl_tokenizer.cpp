/**
 * @file grdecl_tokenizer.cpp
 * @brief GRDECL 行注释、引号与分隔符 token 化实现。
 */
#include "grdecl_detail.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace MPMC::grdecl_detail
{
namespace
{

std::string stripLineComment(const std::string &line)
{
    bool inQuote = false;
    char quote = '\0';
    for (std::size_t i = 0; i + 1 < line.size(); ++i)
    {
        const char c = line[i];
        if (c == '\'' || c == '"')
        {
            if (!inQuote)
            {
                inQuote = true;
                quote = c;
            }
            else if (quote == c)
            {
                inQuote = false;
            }
        }
        if (!inQuote && c == '-' && line[i + 1] == '-')
            return line.substr(0, i);
    }
    return line;
}

std::vector<std::string> tokenize(const std::string &text)
{
    std::vector<std::string> tokens;
    std::string current;
    bool inQuote = false;
    char quote = '\0';

    auto flush = [&]() {
        if (!current.empty())
        {
            tokens.push_back(current);
            current.clear();
        }
    };

    for (char c : text)
    {
        if (inQuote)
        {
            if (c == quote)
            {
                flush();
                inQuote = false;
            }
            else
            {
                current.push_back(c);
            }
            continue;
        }

        if (c == '\'' || c == '"')
        {
            flush();
            inQuote = true;
            quote = c;
        }
        else if (c == '/')
        {
            flush();
            tokens.emplace_back("/");
        }
        else if (std::isspace(static_cast<unsigned char>(c)) || c == ',')
        {
            flush();
        }
        else
        {
            current.push_back(c);
        }
    }
    if (inQuote)
        throw std::runtime_error("GRDECL contains an unterminated quoted string.");
    flush();
    return tokens;
}

} // namespace

std::vector<std::string> readTokens(const std::filesystem::path &file)
{
    std::ifstream stream(file);
    if (!stream)
        throw std::runtime_error("Cannot open GRDECL file: " + file.string());

    std::ostringstream cleaned;
    std::string line;
    while (std::getline(stream, line))
        cleaned << stripLineComment(line) << '\n';
    return tokenize(cleaned.str());
}

} // namespace MPMC::grdecl_detail
