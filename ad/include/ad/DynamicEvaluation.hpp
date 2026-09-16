/**
 * @file DynamicEvaluation.hpp
 * @brief 运行时导数维数的前向自动微分变量。
 */
#pragma once

#include <ad/Evaluation.hpp>
#include <common/small_vector.hpp>

#ifndef NDEBUG
#include <common/valgrind.hpp>
#endif

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace MPMC
{
namespace DenseAd
{

/**
 * @brief 运行时动态导数维数的前向自动微分变量。
 *
 * 该类是 Evaluation 的动态维数偏特化。与固定维数 Evaluation 不同，导数数量在运行时确定，
 * 因而适用于未知量数量可能随模型配置变化的场景，例如后续增加新的相、组分或热力学约束方程。
 *
 * 内部数据按照以下顺序连续存储：
 * - data_[0]：变量函数值；
 * - data_[1 + i]：变量对第 i 个独立未知量的偏导数。
 *
 * @tparam ValueT 函数值及导数的标量类型，通常为 double，也允许嵌套 Evaluation。
 * @tparam staticSize SmallVector 的栈上静态缓冲区容量。
 *
 * @note 本类只负责自动微分数据及基本代数运算，不包含任何具体流动物理或热力学模型。
 * @note 动态维数对象之间进行二元运算时，要求导数维数完全一致。
 */
template <class ValueT, unsigned staticSize>
class Evaluation<ValueT, DynamicSize, staticSize>
{
public:
  /** 动态维数标记，与固定维数 Evaluation 的 numVars 接口保持一致。 */
  static constexpr int numVars = DynamicSize;

  /** 函数值和导数所使用的标量类型。 */
  using ValueType = ValueT;

  /** 内部连续存储类型。 */
  using StorageType = SmallVector<ValueType, staticSize>;

  /**
   * @brief 构造一个不包含导数的零值对象。
   *
   * 默认对象始终包含一个函数值槽位，导数维数为 0，因此默认构造后即可安全赋标量值。
   */
  Evaluation() : data_(1, ValueType{})
  {
  }

  /**
   * @brief 拷贝构造对象。
   *
   * @param[in] other 待拷贝的自动微分对象。
   */
  Evaluation(const Evaluation &other) = default;

  /**
   * @brief 移动构造对象。
   *
   * 移动完成后，将源对象恢复为“0 个导数、函数值为 0”的有效状态，避免后续误用已移动对象时
   * 破坏内部存储不变量。
   *
   * @param[in,out] other 待移动的自动微分对象。
   */
  Evaluation(Evaluation &&other) : data_(std::move(other.data_))
  {
    other.resetToScalarZero_();
  }

  /**
   * @brief 拷贝赋值。
   *
   * @param[in] other 待拷贝的自动微分对象。
   * @return 当前对象引用。
   */
  Evaluation &operator=(const Evaluation &other) = default;

  /**
   * @brief 移动赋值。
   *
   * @param[in,out] other 待移动的自动微分对象。
   * @return 当前对象引用。
   */
  Evaluation &operator=(Evaluation &&other)
  {
    if (this != &other)
    {
      data_ = std::move(other.data_);
      other.resetToScalarZero_();
    }
    return *this;
  }

  /**
   * @brief 构造指定导数维数的零值对象。
   *
   * @param[in] numDerivatives 导数数量，必须大于或等于 0。
   * @throws std::invalid_argument 当 numDerivatives 小于 0 时抛出异常。
   */
  explicit Evaluation(int numDerivatives)
      : data_(storageSize_(numDerivatives), ValueType{})
  {
  }

  /**
   * @brief 构造指定导数维数的常量自动微分对象。
   *
   * 函数值被设置为 c，所有导数初始化为 0。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] numDerivatives 导数数量，必须大于或等于 0。
   * @param[in] c 函数值。
   * @throws std::invalid_argument 当 numDerivatives 小于 0 时抛出异常。
   */
  template <class RhsValueType>
  Evaluation(int numDerivatives, const RhsValueType &c)
      : data_(storageSize_(numDerivatives), ValueType{})
  {
    setValue(c);
    checkDefined_();
  }

  /**
   * @brief 构造指定独立变量对应的自动微分对象。
   *
   * 函数值设置为 c，第 varPos 个导数设置为 1，其余导数设置为 0。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] numDerivatives 总导数数量。
   * @param[in] c 函数值。
   * @param[in] varPos 当前独立变量在导数向量中的位置。
   * @throws std::invalid_argument 当 numDerivatives 小于 0 时抛出异常。
   * @throws std::out_of_range 当 varPos 超出 [0, numDerivatives) 时抛出异常。
   */
  template <class RhsValueType>
  Evaluation(int numDerivatives, const RhsValueType &c, int varPos)
      : data_(storageSize_(numDerivatives), ValueType{})
  {
    validateDerivativeIndex_(varPos);

    setValue(c);
    data_[derivativeStart_() + varPos] = ValueType(1.0);
    checkDefined_();
  }

  /**
   * @brief 返回当前对象的导数维数。
   *
   * @return 独立变量数量，即导数数量。
   */
  int size() const
  {
    return static_cast<int>(data_.size()) - 1;
  }

  /**
   * @brief 将所有一阶导数清零，同时保留函数值不变。
   */
  void clearDerivatives()
  {
    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      data_[i] = ValueType{};
    }
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的零值对象。
   *
   * @param[in] x 用于提供导数维数的参考对象。
   * @return 新创建的自动微分对象。
   */
  static Evaluation createBlank(const Evaluation &x)
  {
    return Evaluation(x.size());
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的常量 0。
   *
   * @param[in] x 用于提供导数维数的参考对象。
   * @return 函数值为 0、全部导数为 0 的对象。
   */
  static Evaluation createConstantZero(const Evaluation &x)
  {
    return Evaluation(x.size(), ValueType{});
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的常量 1。
   *
   * @param[in] x 用于提供导数维数的参考对象。
   * @return 函数值为 1、全部导数为 0 的对象。
   */
  static Evaluation createConstantOne(const Evaluation &x)
  {
    return Evaluation(x.size(), ValueType(1.0));
  }

  /**
   * @brief 禁止在未提供导数总数时创建动态维数独立变量。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] value 独立变量函数值，仅用于保持统一接口。
   * @param[in] varPos 独立变量位置，仅用于保持统一接口。
   * @return 不返回；该重载始终抛出异常。
   * @throws std::logic_error 动态维数对象必须显式提供导数总数。
   */
  template <class RhsValueType>
  static Evaluation createVariable(const RhsValueType &value, int varPos)
  {
    static_cast<void>(value);
    static_cast<void>(varPos);
    throw std::logic_error(
        "Dynamically sized evaluations require the number of derivatives "
        "when creating a variable.");
  }

  /**
   * @brief 创建具有指定导数总数的独立变量。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] numDerivatives 总导数数量。
   * @param[in] value 独立变量函数值。
   * @param[in] varPos 当前独立变量在导数向量中的位置。
   * @return 新创建的独立变量对象。
   */
  template <class RhsValueType>
  static Evaluation createVariable(int numDerivatives,
                                   const RhsValueType &value,
                                   int varPos)
  {
    return Evaluation(numDerivatives, value, varPos);
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的独立变量。
   *
   * @tparam RhsValueType 输入函数值类型。
   * @param[in] x 用于提供导数维数的参考对象。
   * @param[in] value 独立变量函数值。
   * @param[in] varPos 当前独立变量在导数向量中的位置。
   * @return 新创建的独立变量对象。
   */
  template <class RhsValueType>
  static Evaluation createVariable(const Evaluation &x,
                                   const RhsValueType &value,
                                   int varPos)
  {
    return Evaluation(x.size(), value, varPos);
  }

  /**
   * @brief 创建具有指定导数总数的常量对象。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] numDerivatives 总导数数量。
   * @param[in] value 常量函数值。
   * @return 新创建的常量对象。
   */
  template <class RhsValueType>
  static Evaluation createConstant(int numDerivatives,
                                   const RhsValueType &value)
  {
    return Evaluation(numDerivatives, value);
  }

  /**
   * @brief 禁止在未提供导数总数时创建动态维数常量。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] value 常量函数值，仅用于保持统一接口。
   * @return 不返回；该重载始终抛出异常。
   * @throws std::logic_error 动态维数对象必须显式提供导数总数。
   */
  template <class RhsValueType>
  static Evaluation createConstant(const RhsValueType &value)
  {
    static_cast<void>(value);
    throw std::logic_error(
        "Dynamically sized evaluations require the number of derivatives "
        "when creating a constant.");
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的常量对象。
   *
   * @tparam RhsValueType 输入常量类型。
   * @param[in] x 用于提供导数维数的参考对象。
   * @param[in] value 常量函数值。
   * @return 新创建的常量对象。
   */
  template <class RhsValueType>
  static Evaluation createConstant(const Evaluation &x,
                                   const RhsValueType &value)
  {
    return Evaluation(x.size(), value);
  }

  /**
   * @brief 从另一个对象复制全部导数，但保持当前对象的函数值不变。
   *
   * @param[in] other 导数来源对象。
   * @throws std::invalid_argument 当当前对象与 other 的导数维数不一致时抛出异常。
   */
  void copyDerivatives(const Evaluation &other)
  {
    ensureCompatibleSize_(other);

    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      data_[i] = other.data_[i];
    }
  }

  /**
   * @brief 与另一个自动微分对象执行原位加法。
   *
   * @param[in] other 加数。
   * @return 当前对象引用。
   * @throws std::invalid_argument 当当前对象与 other 的导数维数不一致时抛出异常。
   */
  Evaluation &operator+=(const Evaluation &other)
  {
    ensureCompatibleSize_(other);

    for (int i = 0; i < length_(); ++i)
    {
      data_[i] += other.data_[i];
    }
    return *this;
  }

  /**
   * @brief 与普通标量执行原位加法。
   *
   * 标量仅改变函数值，不改变导数。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 加数。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator+=(const RhsValueType &other)
  {
    data_[valuePosition_()] += other;
    return *this;
  }

  /**
   * @brief 与另一个自动微分对象执行原位减法。
   *
   * @param[in] other 减数。
   * @return 当前对象引用。
   * @throws std::invalid_argument 当当前对象与 other 的导数维数不一致时抛出异常。
   */
  Evaluation &operator-=(const Evaluation &other)
  {
    ensureCompatibleSize_(other);

    for (int i = 0; i < length_(); ++i)
    {
      data_[i] -= other.data_[i];
    }
    return *this;
  }

  /**
   * @brief 与普通标量执行原位减法。
   *
   * 标量仅改变函数值，不改变导数。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 减数。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator-=(const RhsValueType &other)
  {
    data_[valuePosition_()] -= other;
    return *this;
  }

  /**
   * @brief 与另一个自动微分对象执行原位乘法。
   *
   * 导数按照乘积法则计算：d(u v) = v du + u dv。
   *
   * @param[in] other 乘数。
   * @return 当前对象引用。
   * @throws std::invalid_argument 当当前对象与 other 的导数维数不一致时抛出异常。
   */
  Evaluation &operator*=(const Evaluation &other)
  {
    ensureCompatibleSize_(other);

    const ValueType u = value();
    const ValueType v = other.value();

    data_[valuePosition_()] = u * v;
    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      const ValueType uPrime = data_[i];
      const ValueType vPrime = other.data_[i];
      data_[i] = uPrime * v + vPrime * u;
    }
    return *this;
  }

  /**
   * @brief 与普通标量执行原位乘法。
   *
   * 函数值及全部导数同时乘以该标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 乘数。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator*=(const RhsValueType &other)
  {
    for (int i = 0; i < length_(); ++i)
    {
      data_[i] *= other;
    }
    return *this;
  }

  /**
   * @brief 与另一个自动微分对象执行原位除法。
   *
   * 导数按照商法则计算：d(u / v) = (v du - u dv) / v^2。
   *
   * @param[in] other 除数。
   * @return 当前对象引用。
   * @throws std::invalid_argument 当当前对象与 other 的导数维数不一致时抛出异常。
   * @note 本函数不额外检查除数是否为 0，其数值行为由 ValueType 的除法规则决定。
   */
  Evaluation &operator/=(const Evaluation &other)
  {
    ensureCompatibleSize_(other);

    const ValueType u = value();
    const ValueType v = other.value();
    const ValueType denominator = v * v;

    for (int i = derivativeStart_(); i < derivativeEnd_(); ++i)
    {
      const ValueType uPrime = data_[i];
      const ValueType vPrime = other.data_[i];
      data_[i] = (v * uPrime - u * vPrime) / denominator;
    }

    data_[valuePosition_()] = u / v;
    return *this;
  }

  /**
   * @brief 与普通标量执行原位除法。
   *
   * 函数值及全部导数同时除以该标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 除数。
   * @return 当前对象引用。
   */
  template <class RhsValueType>
  Evaluation &operator/=(const RhsValueType &other)
  {
    const ValueType reciprocal = ValueType(1.0) / other;
    for (int i = 0; i < length_(); ++i)
    {
      data_[i] *= reciprocal;
    }
    return *this;
  }

  /**
   * @brief 计算两个自动微分对象之和。
   *
   * @param[in] other 加数。
   * @return 运算结果。
   */
  Evaluation operator+(const Evaluation &other) const
  {
    Evaluation result(*this);
    result += other;
    return result;
  }

  /**
   * @brief 计算自动微分对象与普通标量之和。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 加数。
   * @return 运算结果。
   */
  template <class RhsValueType>
  Evaluation operator+(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result += other;
    return result;
  }

  /**
   * @brief 计算两个自动微分对象之差。
   *
   * @param[in] other 减数。
   * @return 运算结果。
   */
  Evaluation operator-(const Evaluation &other) const
  {
    Evaluation result(*this);
    result -= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象减去普通标量的结果。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 减数。
   * @return 运算结果。
   */
  template <class RhsValueType>
  Evaluation operator-(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result -= other;
    return result;
  }

  /**
   * @brief 返回当前自动微分对象的相反数。
   *
   * @return 函数值和全部导数均取反后的对象。
   */
  Evaluation operator-() const
  {
    Evaluation result(*this);
    for (int i = 0; i < length_(); ++i)
    {
      result.data_[i] = -data_[i];
    }
    return result;
  }

  /**
   * @brief 计算两个自动微分对象之积。
   *
   * @param[in] other 乘数。
   * @return 运算结果。
   */
  Evaluation operator*(const Evaluation &other) const
  {
    Evaluation result(*this);
    result *= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象与普通标量之积。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 乘数。
   * @return 运算结果。
   */
  template <class RhsValueType>
  Evaluation operator*(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result *= other;
    return result;
  }

  /**
   * @brief 计算两个自动微分对象之商。
   *
   * @param[in] other 除数。
   * @return 运算结果。
   */
  Evaluation operator/(const Evaluation &other) const
  {
    Evaluation result(*this);
    result /= other;
    return result;
  }

  /**
   * @brief 计算自动微分对象除以普通标量的结果。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 除数。
   * @return 运算结果。
   */
  template <class RhsValueType>
  Evaluation operator/(const RhsValueType &other) const
  {
    Evaluation result(*this);
    result /= other;
    return result;
  }

  /**
   * @brief 将普通标量赋值给当前自动微分对象。
   *
   * 保留当前导数维数，将函数值设置为 other，并将全部导数清零。
   *
   * @tparam RhsValueType 右端标量类型。
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
   * @brief 判断当前对象的函数值是否等于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 若函数值相等则返回 true，否则返回 false。
   */
  template <class RhsValueType>
  bool operator==(const RhsValueType &other) const
  {
    return value() == other;
  }

  /**
   * @brief 判断两个自动微分对象是否完全相同。
   *
   * 函数值和全部导数均参与比较。
   *
   * @param[in] other 待比较对象。
   * @return 若函数值和全部导数均相等则返回 true，否则返回 false。
   */
  bool operator==(const Evaluation &other) const
  {
    if (size() != other.size())
    {
      return false;
    }

    for (int i = 0; i < length_(); ++i)
    {
      if (data_[i] != other.data_[i])
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief 判断两个自动微分对象是否不完全相同。
   *
   * @param[in] other 待比较对象。
   * @return 若存在函数值或导数差异则返回 true，否则返回 false。
   */
  bool operator!=(const Evaluation &other) const
  {
    return !(*this == other);
  }

  /**
   * @brief 判断当前对象的函数值是否不等于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  bool operator!=(const RhsValueType &other) const
  {
    return !(*this == other);
  }

  /**
   * @brief 判断当前函数值是否大于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  bool operator>(const RhsValueType &other) const
  {
    return value() > other;
  }

  /**
   * @brief 按函数值比较两个自动微分对象的大小。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值更大时返回 true。
   */
  bool operator>(const Evaluation &other) const
  {
    ensureCompatibleSize_(other);
    return value() > other.value();
  }

  /**
   * @brief 判断当前函数值是否小于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  bool operator<(const RhsValueType &other) const
  {
    return value() < other;
  }

  /**
   * @brief 按函数值比较两个自动微分对象的大小。
   *
   * @param[in] other 待比较对象。
   * @return 当前函数值更小时返回 true。
   */
  bool operator<(const Evaluation &other) const
  {
    ensureCompatibleSize_(other);
    return value() < other.value();
  }

  /**
   * @brief 判断当前函数值是否大于等于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  bool operator>=(const RhsValueType &other) const
  {
    return value() >= other;
  }

  /**
   * @brief 按函数值判断当前对象是否大于等于另一个自动微分对象。
   *
   * @param[in] other 待比较对象。
   * @return 比较结果。
   */
  bool operator>=(const Evaluation &other) const
  {
    ensureCompatibleSize_(other);
    return value() >= other.value();
  }

  /**
   * @brief 判断当前函数值是否小于等于普通标量。
   *
   * @tparam RhsValueType 右端标量类型。
   * @param[in] other 待比较标量。
   * @return 比较结果。
   */
  template <class RhsValueType>
  bool operator<=(const RhsValueType &other) const
  {
    return value() <= other;
  }

  /**
   * @brief 按函数值判断当前对象是否小于等于另一个自动微分对象。
   *
   * @param[in] other 待比较对象。
   * @return 比较结果。
   */
  bool operator<=(const Evaluation &other) const
  {
    ensureCompatibleSize_(other);
    return value() <= other.value();
  }

  /**
   * @brief 返回函数值的只读引用。
   *
   * @return data_[0] 的只读引用。
   */
  const ValueType &value() const
  {
    return data_[valuePosition_()];
  }

  /**
   * @brief 设置函数值，同时保持全部导数不变。
   *
   * @tparam RhsValueType 输入函数值类型。
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
   * @param[in] varIdx 导数索引，范围为 [0, size())。
   * @return 对应导数的只读引用。
   * @throws std::out_of_range 当 varIdx 超出合法范围时抛出异常。
   */
  const ValueType &derivative(int varIdx) const
  {
    validateDerivativeIndex_(varIdx);
    return data_[derivativeStart_() + varIdx];
  }

  /**
   * @brief 设置指定独立变量对应的一阶导数。
   *
   * @param[in] varIdx 导数索引，范围为 [0, size())。
   * @param[in] derivativeValue 新导数值。
   * @throws std::out_of_range 当 varIdx 超出合法范围时抛出异常。
   */
  void setDerivative(int varIdx, const ValueType &derivativeValue)
  {
    validateDerivativeIndex_(varIdx);
    data_[derivativeStart_() + varIdx] = derivativeValue;
  }

  /**
   * @brief 将内部存储交给序列化器处理。
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
   * @brief 返回函数值与全部导数的副本。
   *
   * 返回顺序为 [value, d0, d1, ..., dn]。该接口用于与固定维数 Evaluation 的 data()
   * 语义保持一致，并为后续通用数组写入逻辑提供统一入口。
   *
   * @return 包含函数值和全部导数的动态数组副本。
   *
   * @note 该函数需要一次动态内存分配，不应在高频 Jacobian 装配路径中重复调用。
   *       性能敏感场景应优先使用 rawData()。
   */
  std::vector<ValueType> data() const
  {
    std::vector<ValueType> result(static_cast<std::size_t>(length_()));
    for (int i = 0; i < length_(); ++i)
    {
      result[static_cast<std::size_t>(i)] = data_[i];
    }
    return result;
  }

  /**
   * @brief 使用动态数组同时设置函数值和全部导数。
   *
   * 输入顺序必须为 [value, d0, d1, ..., dn]。对象的导数维数会自动调整为
   * data.size() - 1。
   *
   * @param[in] data 新的函数值与导数数据。
   * @throws std::invalid_argument 当输入数组为空时抛出异常。
   */
  void setData(const std::vector<ValueType> &data)
  {
    if (data.empty())
    {
      throw std::invalid_argument(
          "Dynamic evaluation data must contain at least one value entry.");
    }

    StorageType newData(data.size());
    for (std::size_t i = 0; i < data.size(); ++i)
    {
      newData[i] = data[i];
    }
    data_ = std::move(newData);
    checkDefined_();
  }

  /**
   * @brief 返回全部导数的副本。
   *
   * @return 按独立变量索引排列的一阶导数动态数组。
   *
   * @note derivatives() 返回独立导数数组的副本。该操作会产生动态内存分配；
   *       高频装配路径建议使用 derivativeData()。
   */
  std::vector<ValueType> derivatives() const
  {
    std::vector<ValueType> result(static_cast<std::size_t>(size()));
    for (int i = 0; i < size(); ++i)
    {
      result[static_cast<std::size_t>(i)] = derivative(i);
    }
    return result;
  }

  /**
   * @brief 返回内部连续数据首地址。
   *
   * 数据布局为 [value, d0, d1, ..., dn]，可直接用于需要连续内存的底层接口。
   *
   * @return 内部数据首地址。
   * @warning 返回指针在当前对象发生赋值、移动或改变导数维数后可能失效。
   */
  ValueType *rawData()
  {
    return &data_[0];
  }

  /**
   * @brief 返回内部连续数据首地址的只读指针。
   *
   * @return 内部数据首地址。
   * @warning 返回指针在当前对象发生赋值、移动或改变导数维数后可能失效。
   */
  const ValueType *rawData() const
  {
    return &data_[0];
  }

  /**
   * @brief 返回导数连续存储区首地址。
   *
   * @return 当导数数量大于 0 时返回第一个导数的地址，否则返回 nullptr。
   * @warning 返回指针在当前对象发生赋值、移动或改变导数维数后可能失效。
   */
  ValueType *derivativeData()
  {
    return size() > 0 ? &data_[derivativeStart_()] : nullptr;
  }

  /**
   * @brief 返回导数连续存储区首地址的只读指针。
   *
   * @return 当导数数量大于 0 时返回第一个导数的地址，否则返回 nullptr。
   * @warning 返回指针在当前对象发生赋值、移动或改变导数维数后可能失效。
   */
  const ValueType *derivativeData() const
  {
    return size() > 0 ? &data_[derivativeStart_()] : nullptr;
  }

protected:
  /**
   * @brief 返回内部总存储长度。
   *
   * @return 函数值槽位与全部导数槽位的总数量。
   */
  int length_() const
  {
    return static_cast<int>(data_.size());
  }

  /**
   * @brief 返回函数值在内部数组中的位置。
   *
   * @return 固定返回 0。
   */
  static constexpr int valuePosition_()
  {
    return 0;
  }

  /**
   * @brief 返回第一个导数在内部数组中的位置。
   *
   * @return 固定返回 1。
   */
  static constexpr int derivativeStart_()
  {
    return 1;
  }

  /**
   * @brief 返回导数区间的右开端点。
   *
   * @return 内部总存储长度。
   */
  int derivativeEnd_() const
  {
    return length_();
  }

  /**
   * @brief 在调试模式下检查内部数据是否均已定义。
   *
   * 本函数同时检查函数值和所有导数。
   */
  void checkDefined_() const
  {
#ifndef NDEBUG
    for (int i = 0; i < length_(); ++i)
    {
      valgrind::checkDefined(data_[i]);
    }
#endif
  }

private:
  /**
   * @brief 检查两个动态 Evaluation 的导数维数是否一致。
   *
   * 动态维数属于运行时状态，不能只依赖 assert 进行保护，因为 Release 构建会移除
   * assert。所有需要逐项访问另一个 Evaluation 存储区的操作都必须先调用本函数，
   * 从而避免维数不匹配导致越界访问。
   *
   * @param[in] other 待检查的另一自动微分对象。
   * @throws std::invalid_argument 当两侧导数维数不一致时抛出异常。
   */
  void ensureCompatibleSize_(const Evaluation &other) const
  {
    if (size() != other.size())
    {
      throw std::invalid_argument(
          "Dynamic evaluations require matching derivative dimensions.");
    }
  }

  /**
   * @brief 检查导数索引是否合法。
   *
   * @param[in] varIdx 待访问的导数索引。
   * @throws std::out_of_range 当索引不位于 [0, size()) 时抛出异常。
   */
  void validateDerivativeIndex_(int varIdx) const
  {
    if (varIdx < 0 || varIdx >= size())
    {
      throw std::out_of_range(
          "Dynamic evaluation derivative index is out of range.");
    }
  }

  /**
   * @brief 根据导数数量计算内部存储长度并检查输入合法性。
   *
   * @param[in] numDerivatives 导数数量。
   * @return 函数值槽位与导数槽位的总数量。
   * @throws std::invalid_argument 当 numDerivatives 小于 0 时抛出异常。
   */
  static std::size_t storageSize_(int numDerivatives)
  {
    if (numDerivatives < 0)
    {
      throw std::invalid_argument(
          "The number of derivatives must be non-negative.");
    }
    return static_cast<std::size_t>(numDerivatives) + 1U;
  }

  /**
   * @brief 将当前对象重置为不包含导数的零值对象。
   *
   * 主要用于恢复移动操作后的源对象，使其继续满足内部数据至少包含一个函数值槽位的不变量。
   */
  void resetToScalarZero_()
  {
    data_ = StorageType(1, ValueType{});
  }

  /** 函数值和全部一阶导数的连续存储区。 */
  StorageType data_;
};

/**
 * @brief 动态导数维数自动微分类型的便捷别名。
 *
 * @tparam Scalar 函数值和导数的标量类型。
 * @tparam staticSize SmallVector 的栈上静态缓冲区容量。
 */
template <class Scalar, unsigned staticSize = 0>
using DynamicEvaluation = Evaluation<Scalar, DynamicSize, staticSize>;

} // namespace DenseAd


} // namespace MPMC
