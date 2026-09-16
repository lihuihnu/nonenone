/**
 * @file primary_variables.hpp
 * @brief Natural 主变量布局、读取和受限更新辅助。
 */
#pragma once

#include <type_traits>

namespace MPMC
{

/**
 * @brief 将解向量中的 double 主变量构造成当前 Indices::ValueType。
 *
 * 对 ScalarIndices 直接返回 double；对 ADIndices 使用统一 Evaluation 主模板，
 * 不再依赖任何 EvaluationXX.hpp。
 */
template <class Indices>
class PrimaryVariables final
{
public:
    using ValueType = typename Indices::ValueType;

    [[nodiscard]] static ValueType make(double value, [[maybe_unused]] int primaryIndex)
    {
        if constexpr (std::is_same_v<ValueType, double>)
        {
            return value;
        }
        else
        {
            return ValueType::createVariable(value, primaryIndex);
        }
    }
};

} // namespace MPMC
