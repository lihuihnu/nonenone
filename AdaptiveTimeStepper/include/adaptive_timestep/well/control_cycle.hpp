/**
 * @file control_cycle.hpp
 * @brief 收敛时间步内的井控制切换与重复求解流程。
 */
#pragma once

#include <well/control.hpp>
#include <well/specification.hpp>
#include <well/state.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace MPMC
{

struct NoWellControlCycle final
{
    void backup() noexcept {}
    void restore() noexcept {}
    [[nodiscard]] bool update() noexcept { return false; }
};

struct NoWellControlSwitchCallback final
{
    template <class Well, class State>
    void operator()(const Well &, const State &, const WellControlUpdate &) const noexcept
    {
    }
};

template <class Well>
struct WellControlSnapshot final
{
    decltype(Well::control) control{};
    double target{0.0};
    decltype(Well::primaryControl) primaryControl{};
    double primaryTarget{0.0};
    bool primaryControlInitialized{false};
};

/**
 * @brief 已收敛时间步内的井控制切换与重复求解循环。
 *
 * StateProvider signature:
 *     void(std::vector<WellState<Indices>>& states)
 *
 * SyncCallback signature:
 *     void(const std::vector<Well>& wells)
 *
 * SwitchCallback signature:
 *     void(const Well&, const WellState<Indices>&, const WellControlUpdate&)
 *
 * The switch callback is invoked immediately after each actual control change,
 * making control transitions easy to audit without coupling the well module to
 * any particular logging/output implementation.
 */
template <
    class Indices,
    class Well,
    class StateProvider,
    class SyncCallback,
    class SwitchCallback = NoWellControlSwitchCallback>
class WellControlCycle final
{
public:
    WellControlCycle(
        std::vector<Well> &wells,
        StateProvider stateProvider,
        SyncCallback syncCallback,
        SwitchCallback switchCallback = {})
        : wells_(wells),
          stateProvider_(std::move(stateProvider)),
          syncCallback_(std::move(syncCallback)),
          switchCallback_(std::move(switchCallback)),
          states_(wells.size())
    {
    }

    void backup()
    {
        snapshots_.clear();
        snapshots_.reserve(wells_.size());

        for (const auto &well : wells_)
        {
            snapshots_.push_back(
                Snapshot{
                    well.control,
                    well.target,
                    well.primaryControl,
                    well.primaryTarget,
                    well.primaryControlInitialized});
        }
    }

    void restore()
    {
        const std::size_t count =
            snapshots_.size() < wells_.size()
                ? snapshots_.size()
                : wells_.size();

        for (std::size_t i = 0; i < count; ++i)
        {
            wells_[i].control = snapshots_[i].control;
            wells_[i].target = snapshots_[i].target;
            wells_[i].primaryControl = snapshots_[i].primaryControl;
            wells_[i].primaryTarget = snapshots_[i].primaryTarget;
            wells_[i].primaryControlInitialized = snapshots_[i].primaryControlInitialized;
        }

        syncCallback_(wells_);
    }

    [[nodiscard]] bool update()
    {
        if (states_.size() != wells_.size())
            states_.resize(wells_.size());

        stateProvider_(states_);

        bool switched = false;
        for (std::size_t i = 0; i < wells_.size(); ++i)
        {
            const WellControlUpdate update = wells_[i].updateControl(states_[i]);
            if (!update.changed)
                continue;

            switched = true;
            switchCallback_(wells_[i], states_[i], update);
        }

        if (switched)
            syncCallback_(wells_);

        return switched;
    }

    [[nodiscard]] const std::vector<WellState<Indices>> &states() const noexcept
    {
        return states_;
    }

private:
    using Snapshot = WellControlSnapshot<Well>;

    std::vector<Well> &wells_;
    StateProvider stateProvider_;
    SyncCallback syncCallback_;
    SwitchCallback switchCallback_;
    std::vector<WellState<Indices>> states_;
    std::vector<Snapshot> snapshots_;
};

} // namespace MPMC
