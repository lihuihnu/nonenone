/**
 * @file exceptions.hpp
 * @brief 自适应时间步失败时使用的异常类型。
 */
#pragma once

#include <stdexcept>
#include <string>

namespace MPMC
{

class AdaptiveTimeStepFailure final : public std::runtime_error
{
public:
    AdaptiveTimeStepFailure(
        std::string message,
        double targetTime,
        double attemptedDt,
        int reasonCode,
        int retryCount)
        : std::runtime_error(std::move(message)),
          targetTime_(targetTime),
          attemptedDt_(attemptedDt),
          reasonCode_(reasonCode),
          retryCount_(retryCount)
    {
    }

    [[nodiscard]] double targetTime() const noexcept { return targetTime_; }
    [[nodiscard]] double attemptedDt() const noexcept { return attemptedDt_; }
    [[nodiscard]] int reasonCode() const noexcept { return reasonCode_; }
    [[nodiscard]] int retryCount() const noexcept { return retryCount_; }

private:
    double targetTime_{0.0};
    double attemptedDt_{0.0};
    int reasonCode_{0};
    int retryCount_{0};
};

} // namespace MPMC
