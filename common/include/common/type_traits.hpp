/**
 * @file type_traits.hpp
 * @brief 项目通用 C++ 类型萃取辅助。
 */
#pragma once

#include <type_traits>

namespace MPMC
{

/**
 * @brief 判断 T 是否是给定类模板 Template 的特化。
 */
template <
    class T,
    template <class...> class Template>
struct is_specialization_of
    : std::false_type
{
};

template <
    template <class...> class Template,
    class... Args>
struct is_specialization_of<
    Template<Args...>,
    Template>
    : std::true_type
{
};

template <
    class T,
    template <class...> class Template>
inline constexpr bool
    is_specialization_of_v =
        is_specialization_of<
            std::remove_cv_t<
                std::remove_reference_t<T>>,
            Template>::value;

} // namespace MPMC
