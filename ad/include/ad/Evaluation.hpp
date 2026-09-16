/**
 * @file Evaluation.hpp
 * @brief 固定导数维数的前向自动微分变量及基础运算。
 */
#pragma once

#ifndef NDEBUG
#include <common/valgrind.hpp>
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace MPMC
{
namespace DenseAd
{

/**
 * @brief 表示运行时动态导数维数的特殊模板参数。
 *
 * 当 Evaluation 的第二个模板参数取该值时，实际类型由 DynamicEvaluation.hpp 中的偏特化实现。
 */
inline constexpr int DynamicSize = -1;

/**
 * @brief 固定导数维数的前向自动微分变量。
 *
 * Evaluation 同时存储函数值和相对于若干独立变量的一阶导数。对于固定维数版本，导数数量
 * 在编译期由 numDerivs 决定，因此内部使用 std::array 连续存储，不发生动态内存分配。
 *
 * 内部数据布局固定为：
 * - data_[0]：函数值；
 * - data_[1 + i]：函数值对第 i 个独立变量的一阶偏导数。
 *
 * @tparam ValueT 函数值及导数使用的标量类型，通常为 double，也允许嵌套 Evaluation。
 * @tparam numDerivs 一阶导数数量，必须为非负整数；DynamicSize 由动态维数偏特化处理。
 * @tparam staticSize 为与动态 Evaluation 保持统一模板签名而保留的参数；固定维数版本不使用该值。
 *
 * @note 本类只负责自动微分数据及基本代数运算，不包含流动、相态或热力学模型逻辑。
 * @note 二元 Evaluation 运算要求两侧类型具有相同的固定导数维数。
 */
template <class ValueT, int numDerivs, unsigned staticSize = 0>
class Evaluation
{
  static_assert(numDerivs >= 0,
                "Fixed-size Evaluation requires a non-negative derivative count.");

public:
  /** 编译期导数维数。 */
  static constexpr int numVars = numDerivs;

  /** 函数值与导数使用的标量类型。 */
  using ValueType = ValueT;

  /** 函数值与全部导数的连续存储类型。 */
  using StorageType =
      std::array<ValueType, static_cast<std::size_t>(numDerivs) + 1U>;

  /** 仅包含一阶导数的固定长度存储类型。 */
  using DerivativeStorageType =
      std::array<ValueType, static_cast<std::size_t>(numDerivs)>;

  /**
   * @brief 默认构造零值自动微分对象。
   *
   * 函数值和全部导数均进行值初始化，因此对于常用数值标量类型其初始值均为 0。
   */
  Evaluation() = default;

  /**
   * @brief 拷贝构造自动微分对象。
   *
   * @param[in] other 待拷贝对象。
   */
  Evaluation(const Evaluation &other) = default;

  /**
   * @brief 移动构造自动微分对象。
   *
   * 固定维数版本内部使用 std::array，移动语义等价于逐元素移动。
   *
   * @param[in] other 待移动对象。
   */
  Evaluation(Evaluation &&other) = default;

  /**
   * @brief 由标量值构造常量自动微分对象。
   *
   * 函数值设置为 c，全部一阶导数设置为 0。构造函数保持隐式，使 Evaluation 可直接参与
   * 与普通标量的算术表达式。
   *
   * @tparam RhsValueType 输入标量类型。
   * @param[in] c 常量函数值。
   */
  template <class RhsValueType>
  Evaluation(const RhsValueType &c)
  {
    setValue(c);
    clearDerivatives();
    checkDefined_();
  }

  /**
   * @brief 构造一个固定维数的独立自动微分变量。
   *
   * 函数值设置为 c，第 varPos 个导数设置为 1，其余导数设置为 0。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] c 独立变量函数值。
   * @param[in] varPos 独立变量在导数向量中的位置。
   * @pre varPos 必须满足 0 <= varPos < size()。
   */
  template <class RhsValueType>
  Evaluation(const RhsValueType &c, int varPos)
  {
    assertValidDerivativeIndex_(varPos);

    setValue(c);
    clearDerivatives();
    data_[derivativeStart_() + varPos] = ValueType(1.0);
    checkDefined_();
  }

  /**
   * @brief 拷贝赋值。
   *
   * @param[in] other 待拷贝对象。
   * @return 当前对象引用。
   */
  Evaluation &operator=(const Evaluation &other) = default;

  /**
   * @brief 移动赋值。
   *
   * @param[in] other 待移动对象。
   * @return 当前对象引用。
   */
  Evaluation &operator=(Evaluation &&other) = default;

  /**
   * @brief 返回一阶导数数量。
   *
   * @return 编译期固定导数数量 numDerivs。
   */
  [[nodiscard]] static constexpr int size() noexcept
  {
    return numDerivs;
  }

  /**
   * @brief 清零全部一阶导数并保留函数值。
   */
  void clearDerivatives()
  {
    std::fill(data_.begin() + derivativeStart_(), data_.end(), ValueType{});
  }

  /**
   * @brief 创建与参考对象同类型的零值对象。
   *
   * 对固定维数 Evaluation 而言，参考对象仅用于统一静态和动态 Evaluation 的工厂接口。
   *
   * @param[in] reference 参考对象。
   * @return 函数值和全部导数均为 0 的新对象。
   */
  [[nodiscard]] static Evaluation createBlank(const Evaluation &reference)
  {
    static_cast<void>(reference);
    return Evaluation();
  }

  /**
   * @brief 创建与参考对象同类型的常量 0。
   *
   * @param[in] reference 参考对象。
   * @return 函数值为 0、全部导数为 0 的新对象。
   */
  [[nodiscard]] static Evaluation createConstantZero(
      const Evaluation &reference)
  {
    static_cast<void>(reference);
    return Evaluation(ValueType{});
  }

  /**
   * @brief 创建与参考对象同类型的常量 1。
   *
   * @param[in] reference 参考对象。
   * @return 函数值为 1、全部导数为 0 的新对象。
   */
  [[nodiscard]] static Evaluation createConstantOne(
      const Evaluation &reference)
  {
    static_cast<void>(reference);
    return Evaluation(ValueType(1.0));
  }

  /**
   * @brief 创建固定维数独立变量。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] value 独立变量函数值。
   * @param[in] varPos 独立变量在导数向量中的位置。
   * @return 新创建的自动微分独立变量。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createVariable(const RhsValueType &value,
                                                 int varPos)
  {
    return Evaluation(value, varPos);
  }

  /**
   * @brief 创建具有显式导数总数校验的固定维数独立变量。
   *
   * 该接口用于与动态维数 Evaluation 保持统一。对固定维数对象，nVars 必须与 numDerivs
   * 完全一致。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] nVars 调用方期望的导数总数。
   * @param[in] value 独立变量函数值。
   * @param[in] varPos 独立变量在导数向量中的位置。
   * @return 新创建的自动微分独立变量。
   * @throws std::logic_error 当 nVars 与当前固定导数维数不一致时抛出异常。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createVariable(int nVars,
                                                 const RhsValueType &value,
                                                 int varPos)
  {
    validateDerivativeCount_(nVars);
    return Evaluation(value, varPos);
  }

  /**
   * @brief 创建与参考对象具有相同类型的独立变量。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] reference 参考对象，仅用于统一静态和动态 Evaluation 接口。
   * @param[in] value 独立变量函数值。
   * @param[in] varPos 独立变量在导数向量中的位置。
   * @return 新创建的自动微分独立变量。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createVariable(const Evaluation &reference,
                                                 const RhsValueType &value,
                                                 int varPos)
  {
    static_cast<void>(reference);
    return Evaluation(value, varPos);
  }

  /**
   * @brief 创建具有显式导数总数校验的固定维数常量。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] nVars 调用方期望的导数总数。
   * @param[in] value 常量函数值。
   * @return 新创建的自动微分常量。
   * @throws std::logic_error 当 nVars 与当前固定导数维数不一致时抛出异常。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createConstant(int nVars,
                                                 const RhsValueType &value)
  {
    validateDerivativeCount_(nVars);
    return Evaluation(value);
  }

  /**
   * @brief 创建固定维数常量。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] value 常量函数值。
   * @return 新创建的自动微分常量。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createConstant(const RhsValueType &value)
  {
    return Evaluation(value);
  }

  /**
   * @brief 创建与参考对象具有相同类型的常量。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] reference 参考对象，仅用于统一静态和动态 Evaluation 接口。
   * @param[in] value 常量函数值。
   * @return 新创建的自动微分常量。
   */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation createConstant(const Evaluation &reference,
                                                 const RhsValueType &value)
  {
    static_cast<void>(reference);
    return Evaluation(value);
  }

  /**
   * @brief 复制另一个对象的全部导数并保持当前函数值不变。
   *
   * @param[in] other 导数来源对象。
   */
  void copyDerivatives(const Evaluation &other)
  {
    std::copy(other.data_.begin() + derivativeStart_(), other.data_.end(),
              data_.begin() + derivativeStart_());
  }

  /**
   * @brief 与另一个自动微分对象原位相加。
   *
   * @param[in] other 右操作数。
   * @return 当前对象引用。
   */
  Evaluation &operator+=(const Evaluation &other)
  {
    for (int i = 0; i < length_(); ++i)
    {
      data_[static_cast<std::size_t>(i)] +=
          other.data_[static_cast<std::size_t>(i)];
    }
    return *this;
  }

  /**
   * @brief 与普通标量原位相加。
   *
   * 标量只改变函数值，不改变任何导数。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator+=(const RhsValueType &other)
  {
    data_[valuePosition_()] += other;
    return *this;
  }

  /**
   * @brief 与另一个自动微分对象原位相减。
   *
   * @param[in] other 右操作数。
   * @return 当前对象引用。
   */
  Evaluation &operator-=(const Evaluation &other)
  {
    for (int i = 0; i < length_(); ++i)
    {
      data_[static_cast<std::size_t>(i)] -=
          other.data_[static_cast<std::size_t>(i)];
    }
    return *this;
  }

  /**
   * @brief 与普通标量原位相减。
   *
   * 标量只改变函数值，不改变任何导数。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator-=(const RhsValueType &other)
  {
    data_[valuePosition_()] -= other;
    return *this;
  }

  /**
   * @brief 与另一个自动微分对象原位相乘。
   *
   * 导数按照乘积法则 d(uv) = v du + u dv 更新。函数值和导数计算均使用运算前的
   * u、v 值，因此即使执行 x *= x 也能够得到正确结果。
   *
   * @param[in] other 右操作数。
   * @return 当前对象引用。
   */
  Evaluation &operator*=(const Evaluation &other)
  {
    const ValueType u = value();
    const ValueType v = other.value();

    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      const ValueType uPrime = data_[static_cast<std::size_t>(i)];
      const ValueType vPrime = other.data_[static_cast<std::size_t>(i)];
      data_[static_cast<std::size_t>(i)] = uPrime * v + vPrime * u;
    }

    data_[valuePosition_()] = u * v;
    return *this;
  }

  /**
   * @brief 与普通标量原位相乘。
   *
   * 函数值和全部导数同时乘以该标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator*=(const RhsValueType &other)
  {
    for (auto &entry : data_)
    {
      entry *= other;
    }
    return *this;
  }

  /**
   * @brief 除以另一个自动微分对象。
   *
   * 导数按照商法则 d(u/v) = (v du - u dv) / v^2 更新。函数值与导数均使用运算前的
   * u、v 值，因此 x /= x 能稳定得到函数值 1 和零导数。
   *
   * @param[in] other 右操作数。
   * @return 当前对象引用。
   * @pre other.value() 对所使用的 ValueType 必须可作为合法除数。
   */
  Evaluation &operator/=(const Evaluation &other)
  {
    const ValueType u = value();
    const ValueType v = other.value();
    const ValueType denominator = v * v;

    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      const ValueType uPrime = data_[static_cast<std::size_t>(i)];
      const ValueType vPrime = other.data_[static_cast<std::size_t>(i)];
      data_[static_cast<std::size_t>(i)] =
          (v * uPrime - u * vPrime) / denominator;
    }

    data_[valuePosition_()] = u / v;
    return *this;
  }

  /**
   * @brief 除以普通标量。
   *
   * 函数值和全部导数同时乘以该标量的倒数。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator/=(const RhsValueType &other)
  {
    const ValueType reciprocal = ValueType(1.0) / other;
    for (auto &entry : data_)
    {
      entry *= reciprocal;
    }
    return *this;
  }

  /**
   * @brief 计算两个自动微分对象之和。
   *
   * @param[in] other 右操作数。
   * @return 加法结果。
   */
  [[nodiscard]] Evaluation operator+(const Evaluation &other) const
  {
    Evaluation result(*this);
    result += other;
    return result;
  }

  /**
   * @brief 计算自动微分对象与普通标量之和。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 加法结果。
   */
  template <class RhsValueType>
  [[nodiscard]] Evaluation operator+(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result += other;
    return result;
  }

  /**
   * @brief 计算两个自动微分对象之差。
   *
   * @param[in] other 右操作数。
   * @return 减法结果。
   */
  [[nodiscard]] Evaluation operator-(const Evaluation &other) const
  {
    Evaluation result(*this);
    result -= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象减去普通标量的结果。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 减法结果。
   */
  template <class RhsValueType>
  [[nodiscard]] Evaluation operator-(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result -= other;
    return result;
  }

  /**
   * @brief 计算当前自动微分对象的相反数。
   *
   * @return 函数值和全部导数均取反后的对象。
   */
  [[nodiscard]] Evaluation operator-() const
  {
    Evaluation result(*this);
    for (auto &entry : result.data_)
    {
      entry = -entry;
    }
    return result;
  }

  /**
   * @brief 计算两个自动微分对象的乘积。
   *
   * @param[in] other 右操作数。
   * @return 乘法结果。
   */
  [[nodiscard]] Evaluation operator*(const Evaluation &other) const
  {
    Evaluation result(*this);
    result *= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象与普通标量的乘积。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 乘法结果。
   */
  template <class RhsValueType>
  [[nodiscard]] Evaluation operator*(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result *= other;
    return result;
  }

  /**
   * @brief 计算两个自动微分对象的商。
   *
   * @param[in] other 右操作数。
   * @return 除法结果。
   */
  [[nodiscard]] Evaluation operator/(const Evaluation &other) const
  {
    Evaluation result(*this);
    result /= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象除以普通标量的结果。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 右侧标量。
   * @return 除法结果。
   */
  template <class RhsValueType>
  [[nodiscard]] Evaluation operator/(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result /= other;
    return result;
  }

  /**
   * @brief 将普通标量赋给当前对象。
   *
   * 赋值后函数值等于 other，全部一阶导数清零。
   *
   * @tparam RhsValueType 输入标量类型。
   * @param[in] other 新函数值。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator=(const RhsValueType &other)
  {
    setValue(other);
    clearDerivatives();
    return *this;
  }

  /**
   * @brief 判断当前对象函数值是否等于普通标量。
   *
   * 该比较仅比较函数值，不比较导数。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 函数值相等时返回 true。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator==(const RhsValueType &other) const
  {
    return value() == other;
  }

  /**
   * @brief 判断两个自动微分对象是否完全相等。
   *
   * 函数值和全部一阶导数均相等时才返回 true。
   *
   * @param[in] other 待比较对象。
   * @return 两个对象完全相等时返回 true。
   */
  [[nodiscard]] bool operator==(const Evaluation &other) const
  {
    return data_ == other.data_;
  }

  /**
   * @brief 判断两个自动微分对象是否不完全相等。
   *
   * @param[in] other 待比较对象。
   * @return 任一函数值或导数不相等时返回 true。
   */
  [[nodiscard]] bool operator!=(const Evaluation &other) const
  {
    return !(*this == other);
  }

  /**
   * @brief 判断当前对象函数值是否不等于普通标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 函数值不相等时返回 true。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator!=(const RhsValueType &other) const
  {
    return !(*this == other);
  }

  /**
   * @brief 比较当前对象函数值是否大于普通标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator>(const RhsValueType &other) const
  {
    return value() > other;
  }

  /**
   * @brief 比较两个自动微分对象的函数值大小。
   *
   * 导数不参与大小比较。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值更大时返回 true。
   */
  [[nodiscard]] bool operator>(const Evaluation &other) const
  {
    return value() > other.value();
  }

  /**
   * @brief 比较当前对象函数值是否小于普通标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator<(const RhsValueType &other) const
  {
    return value() < other;
  }

  /**
   * @brief 比较两个自动微分对象的函数值大小。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值更小时返回 true。
   */
  [[nodiscard]] bool operator<(const Evaluation &other) const
  {
    return value() < other.value();
  }

  /**
   * @brief 比较当前对象函数值是否大于等于普通标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator>=(const RhsValueType &other) const
  {
    return value() >= other;
  }

  /**
   * @brief 比较两个自动微分对象的函数值大小。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值大于等于 other 时返回 true。
   */
  [[nodiscard]] bool operator>=(const Evaluation &other) const
  {
    return value() >= other.value();
  }

  /**
   * @brief 比较当前对象函数值是否小于等于普通标量。
   *
   * @tparam RhsValueType 右操作数类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  [[nodiscard]] bool operator<=(const RhsValueType &other) const
  {
    return value() <= other;
  }

  /**
   * @brief 比较两个自动微分对象的函数值大小。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值小于等于 other 时返回 true。
   */
  [[nodiscard]] bool operator<=(const Evaluation &other) const
  {
    return value() <= other.value();
  }

  /**
   * @brief 返回函数值只读引用。
   *
   * @return 函数值。
   */
  [[nodiscard]] const ValueType &value() const noexcept
  {
    return data_[valuePosition_()];
  }

  /**
   * @brief 设置函数值并保留现有导数。
   *
   * @tparam RhsValueType 输入值类型。
   * @param[in] value 新函数值。
   */
  template <class RhsValueType>
  void setValue(const RhsValueType &value)
  {
    data_[valuePosition_()] = value;
  }

  /**
   * @brief 返回指定独立变量对应的一阶导数。
   *
   * @param[in] varIdx 导数索引。
   * @return 指定一阶导数的只读引用。
   * @pre varIdx 必须满足 0 <= varIdx < size()。
   */
  [[nodiscard]] const ValueType &derivative(int varIdx) const
  {
    assertValidDerivativeIndex_(varIdx);
    return data_[static_cast<std::size_t>(derivativeStart_() + varIdx)];
  }

  /**
   * @brief 设置指定独立变量对应的一阶导数。
   *
   * @param[in] varIdx 导数索引。
   * @param[in] derivativeValue 新导数值。
   * @pre varIdx 必须满足 0 <= varIdx < size()。
   */
  void setDerivative(int varIdx, const ValueType &derivativeValue)
  {
    assertValidDerivativeIndex_(varIdx);
    data_[static_cast<std::size_t>(derivativeStart_() + varIdx)] =
        derivativeValue;
  }

  /**
   * @brief 将内部连续存储交给序列化器处理。
   *
   * @tparam Serializer 序列化器类型。
   * @param[in,out] serializer 序列化器对象。
   */
  template <class Serializer>
  void serializeOp(Serializer &serializer)
  {
    serializer(data_);
  }

  /**
   * @brief 返回函数值和全部导数的副本。
   *
   * 返回顺序为 [value, d0, d1, ..., dn]。采用值返回与当前固定维数特化及动态维数接口的
   * 语义保持一致；性能敏感路径应使用 rawData()。
   *
   * @return 内部数据副本。
   */
  [[nodiscard]] StorageType data() const
  {
    return data_;
  }

  /**
   * @brief 用完整连续数据覆盖当前对象。
   *
   * 输入顺序必须为 [value, d0, d1, ..., dn]。
   *
   * @param[in] data 新的函数值和导数数据。
   */
  void setData(const StorageType &data)
  {
    data_ = data;
    checkDefined_();
  }

  /**
   * @brief 返回全部一阶导数的副本。
   *
   * @return 长度为 numDerivs 的导数数组。
   */
  [[nodiscard]] DerivativeStorageType derivatives() const
  {
    DerivativeStorageType result{};
    std::copy(data_.begin() + derivativeStart_(), data_.end(), result.begin());
    return result;
  }

  /**
   * @brief 返回内部连续数据首地址。
   *
   * 数据布局为 [value, d0, d1, ..., dn]。该接口主要供 Jacobian 装配等性能敏感代码使用。
   *
   * @return 可写连续数据指针。
   */
  [[nodiscard]] ValueType *rawData() noexcept
  {
    return data_.data();
  }

  /**
   * @brief 返回内部连续数据首地址的只读版本。
   *
   * @return 只读连续数据指针。
   */
  [[nodiscard]] const ValueType *rawData() const noexcept
  {
    return data_.data();
  }

  /**
   * @brief 返回第一个一阶导数的地址。
   *
   * 当 numDerivs 为 0 时返回 nullptr；否则返回连续导数区域首地址。
   *
   * @return 可写导数数据指针，或 nullptr。
   */
  [[nodiscard]] ValueType *derivativeData() noexcept
  {
    if constexpr (numDerivs == 0)
    {
      return nullptr;
    }
    else
    {
      return data_.data() + derivativeStart_();
    }
  }

  /**
   * @brief 返回第一个一阶导数地址的只读版本。
   *
   * @return 只读导数数据指针，或 nullptr。
   */
  [[nodiscard]] const ValueType *derivativeData() const noexcept
  {
    if constexpr (numDerivs == 0)
    {
      return nullptr;
    }
    else
    {
      return data_.data() + derivativeStart_();
    }
  }

private:
  /**
   * @brief 返回函数值与全部导数的总存储长度。
   *
   * @return 1 + numDerivs。
   */
  [[nodiscard]] static constexpr int length_() noexcept
  {
    return numDerivs + 1;
  }

  /**
   * @brief 返回函数值在内部数组中的位置。
   *
   * @return 固定返回 0。
   */
  [[nodiscard]] static constexpr std::size_t valuePosition_() noexcept
  {
    return 0U;
  }

  /**
   * @brief 返回第一个导数在内部数组中的位置。
   *
   * @return 固定返回 1。
   */
  [[nodiscard]] static constexpr std::size_t derivativeStart_() noexcept
  {
    return 1U;
  }

  /**
   * @brief 返回导数区间的尾后位置。
   *
   * @return 等于内部总存储长度。
   */
  [[nodiscard]] static constexpr int derivativeEnd_() noexcept
  {
    return length_();
  }

  /**
   * @brief 调试模式下检查内部数据是否均已定义。
   *
   * Release 模式下该函数为空操作。
   */
  void checkDefined_() const
  {
#ifndef NDEBUG
    for (const auto &entry : data_)
    {
      valgrind::checkDefined(entry);
    }
#endif
  }

  /**
   * @brief 检查调用方给出的导数总数是否匹配当前固定维数。
   *
   * @param[in] nVars 调用方给出的导数数量。
   * @throws std::logic_error 当 nVars != numDerivs 时抛出异常。
   */
  static void validateDerivativeCount_(int nVars)
  {
    if (nVars != numDerivs)
    {
      throw std::logic_error(
          "Fixed-size Evaluation derivative count does not match numDerivs.");
    }
  }

  /**
   * @brief 断言导数索引处于合法范围。
   *
   * @param[in] varIdx 待检查索引。
   */
  static void assertValidDerivativeIndex_(int varIdx)
  {
    assert(0 <= varIdx && varIdx < numDerivs);
    static_cast<void>(varIdx);
  }

  /** 函数值和全部一阶导数的连续固定长度存储区。 */
  StorageType data_{};
};

/**
 * @brief 比较普通标量是否小于 Evaluation 的函数值。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 比较结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] bool operator<(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return rhs > lhs;
}

/**
 * @brief 比较普通标量是否大于 Evaluation 的函数值。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 比较结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] bool operator>(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return rhs < lhs;
}

/**
 * @brief 比较普通标量是否小于等于 Evaluation 的函数值。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 比较结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] bool operator<=(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return rhs >= lhs;
}

/**
 * @brief 比较普通标量是否大于等于 Evaluation 的函数值。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 比较结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] bool operator>=(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return rhs <= lhs;
}

/**
 * @brief 判断普通标量是否不等于 Evaluation 的函数值。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 函数值不相等时返回 true。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] bool operator!=(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return lhs != rhs.value();
}

/**
 * @brief 计算普通标量与 Evaluation 之和。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 加法结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize> operator+(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  Evaluation<ValueType, numVars, staticSize> result(rhs);
  result += lhs;
  return result;
}

/**
 * @brief 计算普通标量减去 Evaluation 的结果。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 减法结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize> operator-(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  return -(rhs - lhs);
}

/**
 * @brief 计算普通标量除以 Evaluation 的结果。
 *
 * 通过 createConstant(rhs, lhs) 按右操作数的导数维数创建常量对象，从而同时适用于
 * 固定维数和运行时动态维数 Evaluation。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 除法结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize> operator/(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  using EvaluationType = Evaluation<ValueType, numVars, staticSize>;
  EvaluationType result = EvaluationType::createConstant(rhs, lhs);
  result /= rhs;
  return result;
}

/**
 * @brief 计算普通标量与 Evaluation 的乘积。
 *
 * @tparam RhsValueType 左侧标量类型。
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in] lhs 左侧标量。
 * @param[in] rhs 右侧自动微分对象。
 * @return 乘法结果。
 */
template <class RhsValueType, class ValueType, int numVars,
          unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize> operator*(
    const RhsValueType &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  Evaluation<ValueType, numVars, staticSize> result(rhs);
  result *= lhs;
  return result;
}

/**
 * @brief 判断一个类型是否为 DenseAd::Evaluation。
 *
 * @tparam T 待检查类型。
 */
template <class T>
struct is_evaluation : std::false_type
{
};

/**
 * @brief Evaluation 类型的 is_evaluation 特化。
 *
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 */
template <class ValueType, int numVars, unsigned staticSize>
struct is_evaluation<Evaluation<ValueType, numVars, staticSize>>
    : std::true_type
{
};

/**
 * @brief 去除引用和 cv 限定后判断类型是否为 Evaluation 的便捷常量。
 *
 * @tparam T 待检查类型。
 */
template <class T>
inline constexpr bool is_evaluation_v =
    is_evaluation<std::remove_cv_t<std::remove_reference_t<T>>>::value;

/**
 * @brief 将 Evaluation 输出到指定流。
 *
 * 输出格式固定为 `v: <value>`；当 withDer 为 true 时，继续输出
 * ` / d: <d0> <d1> ...`。函数模板定义保留在头文件中，从而允许任意固定
 * 导数维数和动态导数维数的 Evaluation 直接使用流式输出，而无需在
 * Evaluation.cpp 中维护显式模板实例化列表。
 *
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in,out] os 输出流。
 * @param[in] eval 待输出对象。
 * @param[in] withDer 是否输出全部导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
void printEvaluation(std::ostream &os,
                     const Evaluation<ValueType, numVars, staticSize> &eval,
                     bool withDer = false)
{
  os << "v: " << eval.value();

  if (!withDer)
  {
    return;
  }

  os << " / d:";

  for (int variableIndex = 0; variableIndex < eval.size(); ++variableIndex)
  {
    os << ' ' << eval.derivative(variableIndex);
  }
}

/**
 * @brief 使用流运算符输出 Evaluation。
 *
 * 对普通标量 Evaluation 输出函数值和导数；对于嵌套 Evaluation，仅递归输出其内部函数值，
 * 避免输出结构无限扩张。
 *
 * @tparam ValueType Evaluation 标量类型。
 * @tparam numVars Evaluation 导数维数。
 * @tparam staticSize 动态 Evaluation 的栈上缓冲容量；固定维数 Evaluation 忽略该值。
 * @param[in,out] os 输出流。
 * @param[in] eval 待输出对象。
 * @return 输出流引用。
 */
template <class ValueType, int numVars, unsigned staticSize>
std::ostream &operator<<(
    std::ostream &os,
    const Evaluation<ValueType, numVars, staticSize> &eval)
{
  if constexpr (is_evaluation_v<ValueType>)
  {
    printEvaluation(os, eval.value(), false);
  }
  else
  {
    printEvaluation(os, eval, true);
  }

  return os;
}

} // namespace DenseAd
} // namespace MPMC

/*
 * 固定维数统一使用本文件中的通用 Evaluation 主模板。
 * 这里只引入运行时动态维数的偏特化。
 *
 * DynamicEvaluation.hpp 会再次 include Evaluation.hpp，但由于本文件使用
 * #pragma once，不会形成递归定义；此时固定维数主模板已经完整定义。
 */
#include <ad/DynamicEvaluation.hpp>
