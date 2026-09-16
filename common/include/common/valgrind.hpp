/**
 * @file valgrind.hpp
 * @brief 调试构建下的 Valgrind 数据定义状态辅助。
 */
#pragma once

#include <cstddef>

#if defined(HAVE_VALGRIND) && HAVE_VALGRIND
#include <valgrind/memcheck.h>
#endif

namespace MPMC::valgrind
{

/**
 * @brief 当前进程是否正在 Valgrind 下运行。
 */
[[nodiscard]] inline bool isRunning() noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    return RUNNING_ON_VALGRIND != 0;
#else
    return false;
#endif
}

/**
 * @brief 检查单个对象是否全部初始化。
 */
template <class T>
inline bool checkDefined(
    [[maybe_unused]] const T &value) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    return VALGRIND_CHECK_MEM_IS_DEFINED(
               &value,
               sizeof(T)) == 0;
#else
    return true;
#endif
}

/**
 * @brief 检查连续数组是否全部初始化。
 */
template <class T>
inline bool checkDefined(
    [[maybe_unused]] const T *data,
    [[maybe_unused]] std::size_t count) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    return VALGRIND_CHECK_MEM_IS_DEFINED(
               data,
               count * sizeof(T)) == 0;
#else
    return true;
#endif
}

/**
 * @brief 检查对象内存是否可访问。
 */
template <class T>
inline bool checkAddressable(
    [[maybe_unused]] const T &value) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    return VALGRIND_CHECK_MEM_IS_ADDRESSABLE(
               &value,
               sizeof(T)) == 0;
#else
    return true;
#endif
}

template <class T>
inline void markUndefined(
    [[maybe_unused]] T &value) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(
        &value,
        sizeof(T));
#endif
}

template <class T>
inline void markUndefined(
    [[maybe_unused]] T *data,
    [[maybe_unused]] std::size_t count) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(
        data,
        count * sizeof(T));
#endif
}

template <class T>
inline void markDefined(
    [[maybe_unused]] T &value) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_DEFINED(
        &value,
        sizeof(T));
#endif
}

template <class T>
inline void markDefined(
    [[maybe_unused]] T *data,
    [[maybe_unused]] std::size_t count) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_DEFINED(
        data,
        count * sizeof(T));
#endif
}

template <class T>
inline void markNoAccess(
    [[maybe_unused]] T &value) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_NOACCESS(
        &value,
        sizeof(T));
#endif
}

template <class T>
inline void markNoAccess(
    [[maybe_unused]] T *data,
    [[maybe_unused]] std::size_t count) noexcept
{
#if !defined(NDEBUG) && defined(HAVE_VALGRIND) && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_NOACCESS(
        data,
        count * sizeof(T));
#endif
}

} // namespace MPMC::valgrind
