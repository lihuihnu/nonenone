/**
 * @file manager.hpp
 * @brief 多井集合、状态和控制更新的统一管理器。
 */
#pragma once

#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace MPMC
{

/**
 * @brief 具有值语义的井注册表。
 *
 * Manager 直接拥有井对象并按稳定整数 id 查询；刻意不包含网格、PETSc、
 * schedule 或求解器逻辑，以保持井数据层独立。
 */
template <class Well>
class WellManager final
{
public:
    /** @brief 添加一口井；重复井 id 直接拒绝。 */
    void addWell(Well well)
    {
        const int id = well.id;
        if (idToIndex_.find(id) != idToIndex_.end())
            throw std::invalid_argument("Duplicate well id.");

        idToIndex_.emplace(id, wells_.size());
        wells_.push_back(std::move(well));
    }

    /** @brief 当前注册井数量。 */
    [[nodiscard]] std::size_t size() const noexcept
    {
        return wells_.size();
    }

    /** @brief 按稳定井 id 获取可写井对象。 */
    [[nodiscard]] Well &atId(int id)
    {
        return wells_.at(idToIndex_.at(id));
    }

    /** @brief 按稳定井 id 获取只读井对象。 */
    [[nodiscard]] const Well &atId(int id) const
    {
        return wells_.at(idToIndex_.at(id));
    }

    /** @brief 以插入顺序访问连续的可写井存储。 */
    [[nodiscard]] std::vector<Well> &wells() noexcept
    {
        return wells_;
    }

    /** @brief 以插入顺序访问连续的只读井存储。 */
    [[nodiscard]] const std::vector<Well> &wells() const noexcept
    {
        return wells_;
    }


private:
    std::vector<Well> wells_;
    std::unordered_map<int, std::size_t> idToIndex_;
};

} // namespace MPMC
