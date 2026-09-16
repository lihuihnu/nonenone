/**
 * @file specification.hpp
 * @brief 与网格无关的完整井定义和用户配置接口。
 */
#pragma once

#include <well/control.hpp>
#include <well/perforation.hpp>
#include <well/schedule.hpp>
#include <well/state.hpp>
#include <well/types.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 与网格实现无关的完整井定义。
 *
 * `CellId` 是稳定的整数单元标识。井层只保存物理参数和控制数据；网格归属、
 * PETSc DOF 和 MPI 通信由 Natural runtime 在装配阶段提供。
 *
 * 流量目标内部保存为非负工程量；装配井方程时由 `signedTarget()` 转换为
 * 模拟器符号约定（注入为正、采出为负）。
 */
template <class Indices, class CellId>
struct WellSpecification final
{
    static_assert(std::is_integral_v<CellId>,
                  "WellSpecification requires an integral stable cell id.");

    int id{0};
    std::string name;
    WellType type{WellType::Producer};
    WellControl control{WellControl::Bhp};
    double target{0.0};

    // 主控制是用户指定模式；活动控制可因安全/生产限制临时切换，条件解除后
    // 允许返回主控制。
    WellControl primaryControl{WellControl::Bhp};
    double primaryTarget{0.0};
    bool primaryControlInitialized{false};

    CellId bhpCellId{-1};
    std::vector<WellPerforation<CellId>> perforations;
    std::array<double, Indices::numPhases> injectionPhaseFraction{};
    std::array<double, Indices::numComponents> injectionComponentMassFraction{};

    double radius{0.1};
    double skin{0.0};
    double referenceDepth{0.0};
    WellSchedule schedule{};
    WellControlLimits limits{};

    WellSpecification() = default;

    /** @brief 由主控制方式和穿孔列表构造井。 */
    WellSpecification(
        int wellId,
        std::string wellName,
        WellType wellType,
        WellControl wellControl,
        double controlTarget,
        CellId representativeCell,
        std::vector<WellPerforation<CellId>> wellPerforations)
        : id(wellId),
          name(std::move(wellName)),
          type(wellType),
          control(wellControl),
          target(isRateControl(wellControl) ? std::abs(controlTarget) : controlTarget),
          primaryControl(wellControl),
          primaryTarget(isRateControl(wellControl) ? std::abs(controlTarget) : controlTarget),
          primaryControlInitialized(true),
          bhpCellId(representativeCell),
          perforations(std::move(wellPerforations))
    {
    }

    /** @brief 当前活动井方程为地面流量控制而非 BHP 控制时返回 true。 */
    [[nodiscard]] bool isRateControlled() const noexcept
    {
        return isRateControl(control);
    }

    /** @brief 返回采用模拟器符号约定的活动控制目标。 */
    [[nodiscard]] double signedTarget() const
    {
        return isRateControlled() ? signedRateTarget(type, target) : target;
    }

    /** @brief 给定物理时间位于井的启用 schedule 内时返回 true。 */
    [[nodiscard]] bool isActive(double time) const noexcept
    {
        return schedule.isActive(time);
    }

    /** @brief 设置注入井最大 BHP 安全限制 [Pa]。 */
    void setMaximumBhp(double value)
    {
        auto updated = limits;
        updated.maximumBhp = value;
        updated.validate();
        limits = std::move(updated);
    }

    /** @brief 设置生产井最小 BHP 安全限制 [Pa]。 */
    void setMinimumBhp(double value)
    {
        auto updated = limits;
        updated.minimumBhp = value;
        updated.validate();
        limits = std::move(updated);
    }

    /** @brief 设置最大产水地面流量工程量。 */
    void setMaximumWaterRate(double value)
    {
        auto updated = limits;
        updated.maximumWaterRate = value;
        updated.validate();
        limits = std::move(updated);
    }

    /** @brief 设置自动控制切换使用的压力和相对流量滞回。 */
    void setControlSwitchTolerances(double pressureTolerance, double rateRelativeTolerance)
    {
        auto updated = limits;
        updated.pressureTolerance = pressureTolerance;
        updated.rateRelativeTolerance = rateRelativeTolerance;
        updated.validate();
        limits = std::move(updated);
    }

    /**
     * @brief 在储层/井方程收敛后更新活动控制方式。
     *
     * The returned record identifies whether a switch occurred, old/new control,
     * trigger reason, measured value and violated limit. No reservoir state is
     * changed here; a changed control simply requests another nonlinear solve.
     */
    [[nodiscard]] WellControlUpdate updateControl(const WellState<Indices> &state)
    {
        if (!primaryControlInitialized)
            capturePrimaryControl_();

        WellControlSelection selection{
            primaryControl, primaryTarget, control, target};
        const WellControlUpdate update =
            updateWellControl(type, limits, state, selection);
        control = selection.activeControl;
        target = selection.activeTarget;
        return update;
    }

    /**
     * @brief 校验井 id、网格无关参数、schedule 和注入相/组分分数。
     *
     * The representative BHP cell must be one of the perforations. Fraction arrays
     * are intentionally not normalized silently: invalid case input is rejected.
     */
    void validate(CellId globalCellCount) const
    {
        if (globalCellCount <= 0)
            throw std::invalid_argument("Global cell count must be positive.");

        schedule.validate();
        limits.validate();
        if (bhpCellId < 0 || bhpCellId >= globalCellCount)
            throw std::invalid_argument("Well BHP representative cell id is outside the grid.");
        if (perforations.empty())
            throw std::invalid_argument("A well requires at least one perforation.");

        bool representativeFound = false;
        std::unordered_set<CellId> uniqueCells;
        uniqueCells.reserve(perforations.size());
        for (const auto &perforation : perforations)
        {
            perforation.validate();
            if (perforation.currentCellId >= globalCellCount)
                throw std::invalid_argument("Well perforation cell id is outside the grid.");
            if (!uniqueCells.insert(perforation.currentCellId).second)
                throw std::invalid_argument("A well cannot contain duplicate perforation cells.");
            representativeFound |= perforation.currentCellId == bhpCellId;
        }
        if (!representativeFound)
            throw std::invalid_argument("Well BHP representative cell must be one of its perforations.");

        if (!std::isfinite(target))
            throw std::invalid_argument("Well target must be finite.");
        if (!std::isfinite(radius) || radius <= 0.0)
            throw std::invalid_argument("Well radius must be positive and finite.");
        if (!std::isfinite(skin))
            throw std::invalid_argument("Well skin must be finite.");
        if (!std::isfinite(referenceDepth))
            throw std::invalid_argument("Well reference depth must be finite.");
        if (control == WellControl::WaterRate && !Indices::hasWater)
            throw std::invalid_argument("Water-rate control requires a water phase.");

        if (type == WellType::Injector)
        {
            validatePositiveFractions_(
                injectionPhaseFraction,
                "Injector phase fractions must be finite and non-negative.",
                "Injector requires a positive phase-fraction sum.");
            validatePositiveFractions_(
                injectionComponentMassFraction,
                "Injector component fractions must be finite and non-negative.",
                "Injector requires a positive component-fraction sum.");
        }
    }

private:
    void capturePrimaryControl_()
    {
        primaryControl = control;
        primaryTarget = isRateControl(control) ? std::abs(target) : target;
        target = isRateControl(control) ? std::abs(target) : target;
        primaryControlInitialized = true;
    }

    template <class Array>
    static void validatePositiveFractions_(
        const Array &fractions,
        const char *invalidMessage,
        const char *zeroMessage)
    {
        double sum = 0.0;
        for (double value : fractions)
        {
            if (value < 0.0 || !std::isfinite(value))
                throw std::invalid_argument(invalidMessage);
            sum += value;
        }
        if (!(sum > 0.0))
            throw std::invalid_argument(zeroMessage);
    }
};

} // namespace MPMC
