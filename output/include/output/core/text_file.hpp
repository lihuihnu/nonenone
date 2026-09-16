/**
 * @file text_file.hpp
 * @brief 文本输出文件的创建、覆盖和目录管理。
 */
#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace MPMC
{

/**
 * @brief 可选 rank-local 文本文件的轻量 RAII 封装。
 *
 * `open(path, false)` keeps the object disabled; this lets all MPI ranks share
 * the same output code while only the writer rank owns an actual file stream.
 */
class OutputTextFile final
{
public:
    OutputTextFile() = default;

    /** @brief 启用时创建父目录并打开/截断 `path`。 */
    void open(const std::filesystem::path &path, bool enabled)
    {
        enabled_ = enabled;
        if (!enabled_)
            return;

        const auto parent = path.parent_path();
        if (!parent.empty())
            std::filesystem::create_directories(parent);

        stream_.open(path, std::ios::out | std::ios::trunc);
        if (!stream_.is_open())
            throw std::runtime_error("Failed to open output file: " + path.string());
    }

    /** @brief 在启用的 writer 进程上返回可写流。 */
    std::ostream &stream()
    {
        if (!enabled_ || !stream_.is_open())
            throw std::logic_error("Output file is not open on this process.");
        return stream_;
    }

private:
    bool enabled_{false};
    std::ofstream stream_;
};

} // namespace MPMC
