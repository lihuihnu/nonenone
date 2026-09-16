/**
 * @file Math.hpp
 * @brief 自动微分数学函数及一阶导数传播规则。
 */
#pragma once

#include <ad/Evaluation.hpp>
#include <common/math.hpp>

#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace MPMC
{
namespace DenseAd
{
namespace detail
{

/**
 * @brief 检查两个 Evaluation 的运行时导数维数是否一致。
 *
 * 固定维数 Evaluation 的维数由模板参数保证一致，因此编译器可以消除该检查；
 * 动态维数 Evaluation 则必须在运行时确认两侧导数向量具有相同长度。
 *
 * @tparam ValueType Evaluation 的标量类型。
 * @tparam numVars 编译期导数维数；DynamicSize 表示运行时动态维数。
 * @tparam staticSize 动态存储的小对象缓冲区容量。
 * @param[in] lhs 左侧自动微分对象。
 * @param[in] rhs 右侧自动微分对象。
 * @throws std::invalid_argument 动态维数对象的导数数量不一致时抛出异常。
 */
template <class ValueType, int numVars, unsigned staticSize>
void ensureCompatibleDerivativeSize(
    const Evaluation<ValueType, numVars, staticSize> &lhs,
    const Evaluation<ValueType, numVars, staticSize> &rhs)
{
  if constexpr (numVars == DynamicSize)
  {
    if (lhs.size() != rhs.size())
    {
      throw std::invalid_argument(
          "DenseAd math operation requires matching derivative dimensions.");
    }
  }
}

/**
 * @brief 将无符号导数数量安全转换为 Evaluation 使用的 int 类型。
 *
 * @param[in] numDerivatives 导数数量。
 * @return 转换后的导数数量。
 * @throws std::overflow_error 当导数数量超过 int 可表示范围时抛出异常。
 */
inline int checkedDerivativeCount(unsigned numDerivatives)
{
  if (numDerivatives >
      static_cast<unsigned>(std::numeric_limits<int>::max()))
  {
    throw std::overflow_error(
        "Derivative dimension exceeds the supported integer range.");
  }

  return static_cast<int>(numDerivatives);
}


/**
 * @brief 判断一个标量或嵌套 Evaluation 是否严格为零。
 *
 * 对普通数值类型直接与零比较；对嵌套 Evaluation 则递归检查函数值和全部导数。
 * 该辅助函数用于幂函数的边界分支：只有当某一导数方向严格为零时，才可以安全跳过
 * 该方向对应的偏导项，从而避免出现数学上应为 0、数值实现却形成 0*Inf 或 0*NaN
 * 的伪 NaN。
 *
 * @tparam T 待检查类型。
 * @param[in] value 待检查对象。
 * @return 函数值及全部嵌套导数均严格为零时返回 true。
 */
template <class T>
[[nodiscard]] bool isExactlyZero(const T &value)
{
  if constexpr (is_evaluation_v<T>)
  {
    if (!isExactlyZero(value.value()))
    {
      return false;
    }

    for (int i = 0; i < value.size(); ++i)
    {
      if (!isExactlyZero(value.derivative(i)))
      {
        return false;
      }
    }
    return true;
  }
  else
  {
    return value == T{};
  }
}

/**
 * @brief 创建与输入对象具有相同内部结构的严格零值。
 *
 * 对普通标量等价于返回 0；对嵌套或动态 Evaluation，通过复制后乘以 0 保留其
 * 内部导数维数，避免使用默认构造零值改变嵌套动态 AD 的存储结构。
 *
 * @tparam T 输入类型。
 * @param[in] value 用于提供内部结构的参考对象。
 * @return 与 value 结构一致的零值对象。
 */
template <class T>
[[nodiscard]] T zeroLike(const T &value)
{
  T result(value);
  result *= 0.0;
  return result;
}

} // namespace detail

/**
 * @brief 计算 Evaluation 的绝对值。
 *
 * 绝对值在零点不可微。这里在 x == 0 时选择非负分支，从而保留输入对象当前的导数方向；
 * 该选择主要用于保证分段模型在数值迭代中的确定性，并不代表零点存在唯一数学导数。
 *
 * @tparam ValueType Evaluation 的标量类型。
 * @tparam numVars 导数维数。
 * @tparam staticSize 动态存储的小对象缓冲区容量。
 * @param[in] x 输入自动微分对象。
 * @return x 的绝对值及所选分支对应的一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
abs(const Evaluation<ValueType, numVars, staticSize> &x)
{
  return (x >= 0.0) ? x : -x;
}

/**
 * @brief 返回两个 Evaluation 中函数值较小的对象。
 *
 * 对动态维数对象会首先检查导数维数是否一致。函数值相等时选择第二个参数，
 * 因而在不可微的分支交点处采用第二个分支的导数。
 *
 * @param[in] x1 第一个自动微分对象。
 * @param[in] x2 第二个自动微分对象。
 * @return 函数值较小的对象及其导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
min(const Evaluation<ValueType, numVars, staticSize> &x1,
    const Evaluation<ValueType, numVars, staticSize> &x2)
{
  detail::ensureCompatibleDerivativeSize(x1, x2);
  return (x1 < x2) ? x1 : x2;
}

/**
 * @brief 返回普通标量与 Evaluation 中较小的值。
 *
 * 当普通标量分支被选中时，返回与 x2 具有相同导数维数的常量对象，所有导数均为 0。
 *
 * @tparam Scalar 左侧普通标量类型。
 * @param[in] x1 普通标量。
 * @param[in] x2 自动微分对象。
 * @return 两者中的较小值。
 */
template <class Scalar, class ValueType, int numVars, unsigned staticSize,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
min(const Scalar &x1, const Evaluation<ValueType, numVars, staticSize> &x2)
{
  if (x1 < x2)
  {
    return Evaluation<ValueType, numVars, staticSize>::createConstant(x2, x1);
  }

  return x2;
}

/**
 * @brief 返回 Evaluation 与普通标量中较小的值。
 *
 * @tparam Scalar 右侧普通标量类型。
 * @param[in] x1 自动微分对象。
 * @param[in] x2 普通标量。
 * @return 两者中的较小值。
 */
template <class ValueType, int numVars, unsigned staticSize, class Scalar,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
min(const Evaluation<ValueType, numVars, staticSize> &x1, const Scalar &x2)
{
  return min(x2, x1);
}

/**
 * @brief 返回两个 Evaluation 中函数值较大的对象。
 *
 * 对动态维数对象会首先检查导数维数是否一致。函数值相等时选择第二个参数，
 * 因而在不可微的分支交点处采用第二个分支的导数。
 *
 * @param[in] x1 第一个自动微分对象。
 * @param[in] x2 第二个自动微分对象。
 * @return 函数值较大的对象及其导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
max(const Evaluation<ValueType, numVars, staticSize> &x1,
    const Evaluation<ValueType, numVars, staticSize> &x2)
{
  detail::ensureCompatibleDerivativeSize(x1, x2);
  return (x1 > x2) ? x1 : x2;
}

/**
 * @brief 返回普通标量与 Evaluation 中较大的值。
 *
 * 当普通标量分支被选中时，返回与 x2 具有相同导数维数的常量对象，所有导数均为 0。
 *
 * @tparam Scalar 左侧普通标量类型。
 * @param[in] x1 普通标量。
 * @param[in] x2 自动微分对象。
 * @return 两者中的较大值。
 */
template <class Scalar, class ValueType, int numVars, unsigned staticSize,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
max(const Scalar &x1, const Evaluation<ValueType, numVars, staticSize> &x2)
{
  if (x1 > x2)
  {
    return Evaluation<ValueType, numVars, staticSize>::createConstant(x2, x1);
  }

  return x2;
}

/**
 * @brief 返回 Evaluation 与普通标量中较大的值。
 *
 * @tparam Scalar 右侧普通标量类型。
 * @param[in] x1 自动微分对象。
 * @param[in] x2 普通标量。
 * @return 两者中的较大值。
 */
template <class ValueType, int numVars, unsigned staticSize, class Scalar,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
max(const Evaluation<ValueType, numVars, staticSize> &x1, const Scalar &x2)
{
  return max(x2, x1);
}

/**
 * @brief 计算正切函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return tan(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
tan(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  const ValueType tanValue = Toolbox::tan(x.value());
  result.setValue(tanValue);

  const ValueType derivativeFactor = 1.0 + tanValue * tanValue;
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算反正切函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return atan(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
atan(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::atan(x.value()));

  const ValueType derivativeFactor =
      1.0 / (1.0 + x.value() * x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算 atan2(y, x) 并传播两个 Evaluation 的一阶导数。
 *
 * 导数采用稳定形式
 * d atan2(y,x) = (x dy - y dx) / (x^2 + y^2)，
 * 该形式只在 `(x,y)=(0,0)` 这一数学奇点处无定义。
 *
 * @param[in] y atan2 的第一参数。
 * @param[in] x atan2 的第二参数。
 * @return atan2(y, x) 及其一阶导数。
 * @note 当 x == 0 且 y == 0 时 atan2 的导数本身未定义。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
atan2(const Evaluation<ValueType, numVars, staticSize> &y,
      const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  detail::ensureCompatibleDerivativeSize(y, x);

  Evaluation<ValueType, numVars, staticSize> result(y);
  result.setValue(Toolbox::atan2(y.value(), x.value()));

  const ValueType denominator =
      x.value() * x.value() + y.value() * y.value();
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(
        i, (x.value() * y.derivative(i) - y.value() * x.derivative(i)) /
               denominator);
  }

  return result;
}

/**
 * @brief 计算 atan2(y, x)；其中第一参数为 Evaluation，第二参数为普通标量。
 *
 * @tparam Scalar 第二参数的普通标量类型。
 * @param[in] y atan2 的第一参数。
 * @param[in] x atan2 的第二参数。
 * @return atan2(y, x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize, class Scalar,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
atan2(const Evaluation<ValueType, numVars, staticSize> &y, const Scalar &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(y);
  result.setValue(Toolbox::atan2(y.value(), x));

  const ValueType denominator = x * x + y.value() * y.value();
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, x * y.derivative(i) / denominator);
  }

  return result;
}

/**
 * @brief 计算 atan2(y, x)；其中第一参数为普通标量，第二参数为 Evaluation。
 *
 * @tparam Scalar 第一参数的普通标量类型。
 * @param[in] y atan2 的第一参数。
 * @param[in] x atan2 的第二参数。
 * @return atan2(y, x) 及其一阶导数。
 */
template <class Scalar, class ValueType, int numVars, unsigned staticSize,
          std::enable_if_t<!is_evaluation_v<Scalar>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
atan2(const Scalar &y, const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::atan2(y, x.value()));

  const ValueType denominator =
      x.value() * x.value() + y * y;
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, -y * x.derivative(i) / denominator);
  }

  return result;
}

/**
 * @brief 计算正弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return sin(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
sin(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::sin(x.value()));

  const ValueType derivativeFactor = Toolbox::cos(x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算反正弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return asin(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
asin(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::asin(x.value()));

  const ValueType derivativeFactor =
      1.0 / Toolbox::sqrt(1.0 - x.value() * x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算双曲正弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return sinh(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
sinh(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::sinh(x.value()));

  const ValueType derivativeFactor = Toolbox::cosh(x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算反双曲正弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return asinh(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
asinh(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::asinh(x.value()));

  const ValueType derivativeFactor =
      1.0 / Toolbox::sqrt(x.value() * x.value() + 1.0);
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算余弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return cos(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
cos(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::cos(x.value()));

  const ValueType derivativeFactor = -Toolbox::sin(x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算反余弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return acos(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
acos(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::acos(x.value()));

  const ValueType derivativeFactor =
      -1.0 / Toolbox::sqrt(1.0 - x.value() * x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算双曲余弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return cosh(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
cosh(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::cosh(x.value()));

  const ValueType derivativeFactor = Toolbox::sinh(x.value());
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算反双曲余弦函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return acosh(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
acosh(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::acosh(x.value()));

  const ValueType derivativeFactor =
      1.0 / Toolbox::sqrt(x.value() * x.value() - 1.0);
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算平方根并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return sqrt(x) 及其一阶导数。
 * @note x == 0 时平方根导数趋于无穷；x < 0 时实数域结果无效。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
sqrt(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  const ValueType sqrtValue = Toolbox::sqrt(x.value());
  result.setValue(sqrtValue);

  const ValueType derivativeFactor = 0.5 / sqrtValue;
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算指数函数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return exp(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
exp(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  const ValueType expValue = Toolbox::exp(x.value());
  result.setValue(expValue);

  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, expValue * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算 Evaluation 底数的标量幂并传播一阶导数。
 *
 * 使用 d(a^p)/da = p*a^(p-1)。指数为 0 时显式返回零导数，
 * 避免通过 `a^p/a` 引入不必要的 0/0。
 *
 * @tparam Exponent 指数的普通标量类型。
 * @param[in] base 自动微分底数。
 * @param[in] exponent 普通标量指数。
 * @return base^exponent 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize, class Exponent,
          std::enable_if_t<!is_evaluation_v<Exponent>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
pow(const Evaluation<ValueType, numVars, staticSize> &base,
    const Exponent &exponent)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(base);
  const ValueType powValue = Toolbox::pow(base.value(), exponent);
  result.setValue(powValue);

  if (exponent == Exponent(0))
  {
    result.clearDerivatives();
    return result;
  }

  const ValueType derivativeFactor =
      exponent * Toolbox::pow(base.value(), exponent - Exponent(1));
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * base.derivative(i));
  }

  return result;
}

/**
 * @brief 计算普通标量底数的 Evaluation 指数幂并传播一阶导数。
 *
 * 对 y = a^g，指数方向的一阶导数为 y*ln(a)*dg。在实数连续可微意义下，
 * 当 dg 非零时通常要求 a > 0。实现中仅在当前导数方向确实非零时计算 ln(a)：
 *
 * - dg == 0：该方向导数严格为 0，不计算 ln(a)；
 * - a == 0 且 g > 0：0^g 在指数方向局部恒为 0，因此导数取 0；
 * - 其余不满足实数可微条件的情况继续交由底层数学函数产生 NaN/Inf，避免伪造导数。
 *
 * 该分支可正确处理“负底数 + 常量整数 Evaluation 指数”，不会再因为无意义地计算
 * log(负数) 而把本应为 0 的指数方向导数污染为 NaN。
 *
 * @tparam Base 普通标量底数类型。
 * @param[in] base 普通标量底数。
 * @param[in] exponent 自动微分指数。
 * @return base^exponent 及其一阶导数。
 */
template <class Base, class ValueType, int numVars, unsigned staticSize,
          std::enable_if_t<!is_evaluation_v<Base>, int> = 0>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
pow(const Base &base,
    const Evaluation<ValueType, numVars, staticSize> &exponent)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(exponent);
  const ValueType powValue = Toolbox::pow(base, exponent.value());
  result.setValue(powValue);

  const bool baseIsZero = (base == Base{});
  const auto scalarExponent = Toolbox::scalarValue(exponent.value());

  for (int i = 0; i < result.size(); ++i)
  {
    const ValueType &exponentDerivative = exponent.derivative(i);

    if (detail::isExactlyZero(exponentDerivative))
    {
      result.setDerivative(i, detail::zeroLike(exponentDerivative));
      continue;
    }

    if (baseIsZero && scalarExponent > 0)
    {
      result.setDerivative(i, detail::zeroLike(exponentDerivative));
      continue;
    }

    const ValueType derivativeFactor = Toolbox::log(base) * powValue;
    result.setDerivative(i, derivativeFactor * exponentDerivative);
  }

  return result;
}

/**
 * @brief 计算 Evaluation 底数与 Evaluation 指数的幂并传播一阶导数。
 *
 * 对 y = f^g，将每个独立变量方向拆成两个偏导贡献：
 *
 * d y = g*f^(g-1)*d f + f^g*ln(f)*d g。
 *
 * 与直接使用 y*(g/f*df + ln(f)*dg) 相比，该写法不会在 f == 0 时无条件引入
 * g/f。实现还会跳过严格为零的导数贡献，因此能够正确处理：
 *
 * - f == 0、g 为正常正指数时的有限底数方向导数；
 * - 负底数 + 常量整数 Evaluation 指数；
 * - 某一方向 dg == 0 时无需计算 ln(f)，避免 0*NaN 污染结果。
 *
 * 当 dg 非零且 f < 0 时，实数域下一般不存在关于连续指数的导数；此时底层 log/pow
 * 仍会自然产生 NaN/Inf。本函数只消除“数学导数本来有限，却因无条件计算无效项产生
 * NaN”的情况，不人为定义不存在的导数。
 *
 * @param[in] base 自动微分底数。
 * @param[in] exponent 自动微分指数。
 * @return base^exponent 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
pow(const Evaluation<ValueType, numVars, staticSize> &base,
    const Evaluation<ValueType, numVars, staticSize> &exponent)
{
  using Toolbox = MathToolbox<ValueType>;

  detail::ensureCompatibleDerivativeSize(base, exponent);

  Evaluation<ValueType, numVars, staticSize> result(base);
  const ValueType powValue = Toolbox::pow(base.value(), exponent.value());
  result.setValue(powValue);

  const auto scalarBase = Toolbox::scalarValue(base.value());
  const auto scalarExponent = Toolbox::scalarValue(exponent.value());

  for (int i = 0; i < result.size(); ++i)
  {
    const ValueType &baseDerivative = base.derivative(i);
    const ValueType &exponentDerivative = exponent.derivative(i);
    ValueType derivative = detail::zeroLike(baseDerivative);

    if (!detail::isExactlyZero(baseDerivative) &&
        !detail::isExactlyZero(exponent.value()))
    {
      const ValueType baseDerivativeFactor =
          exponent.value() *
          Toolbox::pow(base.value(), exponent.value() - 1.0);
      derivative += baseDerivativeFactor * baseDerivative;
    }

    if (!detail::isExactlyZero(exponentDerivative))
    {
      if (!(scalarBase == 0 && scalarExponent > 0))
      {
        derivative +=
            powValue * Toolbox::log(base.value()) * exponentDerivative;
      }
    }

    result.setDerivative(i, derivative);
  }

  return result;
}

/**
 * @brief 计算自然对数并传播一阶导数。
 *
 * @param[in] x 输入自动微分对象。
 * @return log(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
log(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::log(x.value()));

  const ValueType derivativeFactor = 1.0 / x.value();
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

/**
 * @brief 计算以 10 为底的对数并传播一阶导数。
 *
 * 导数为 1 / (x ln(10))。通过底层 MathToolbox 计算换底系数，以保持 ValueType
 * 可能为嵌套自动微分类型时的泛型能力。
 *
 * @param[in] x 输入自动微分对象。
 * @return log10(x) 及其一阶导数。
 */
template <class ValueType, int numVars, unsigned staticSize>
[[nodiscard]] Evaluation<ValueType, numVars, staticSize>
log10(const Evaluation<ValueType, numVars, staticSize> &x)
{
  using Toolbox = MathToolbox<ValueType>;

  Evaluation<ValueType, numVars, staticSize> result(x);
  result.setValue(Toolbox::log10(x.value()));

  const ValueType derivativeFactor =
      (1.0 / x.value()) * Toolbox::log10(Toolbox::exp(1.0));
  for (int i = 0; i < result.size(); ++i)
  {
    result.setDerivative(i, derivativeFactor * x.derivative(i));
  }

  return result;
}

} // namespace DenseAd

/**
 * @brief DenseAd::Evaluation 对 MathToolbox 的特化。
 *
 * 该特化是通用物性、EOS 和离散代码访问自动微分数学运算的统一入口。
 * 所有导数遍历均基于 eval.size()，因此同时支持固定导数维数 Evaluation 和
 * 运行时动态维数 DynamicEvaluation。
 *
 * @tparam ValueT Evaluation 内部标量类型。
 * @tparam numVars 编译期导数维数；DynamicSize 表示动态维数。
 * @tparam staticSize 动态存储的小对象缓冲区容量。
 */
template <class ValueT, int numVars, unsigned staticSize>
struct MathToolbox<DenseAd::Evaluation<ValueT, numVars, staticSize>>
{
  using ValueType = ValueT;
  using InnerToolbox = MathToolbox<ValueType>;
  using Scalar = typename InnerToolbox::Scalar;
  using Evaluation = DenseAd::Evaluation<ValueType, numVars, staticSize>;

  /**
   * @brief 返回 Evaluation 的直接函数值。
   *
   * @param[in] eval 输入自动微分对象。
   * @return 外层 Evaluation 的 value 部分。
   */
  [[nodiscard]] static ValueType value(const Evaluation &eval)
  {
    return eval.value();
  }

  /**
   * @brief 递归提取最内层普通标量值。
   *
   * @param[in] eval 输入自动微分对象。
   * @return 最内层普通标量值。
   */
  [[nodiscard]] static Scalar scalarValue(const Evaluation &eval)
  {
    return InnerToolbox::scalarValue(eval.value());
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的空白对象。
   *
   * @param[in] x 参考对象。
   * @return 新的空白 Evaluation。
   */
  [[nodiscard]] static Evaluation createBlank(const Evaluation &x)
  {
    return Evaluation::createBlank(x);
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的常量 0。
   *
   * @param[in] x 参考对象。
   * @return 常量 0 Evaluation。
   */
  [[nodiscard]] static Evaluation createConstantZero(const Evaluation &x)
  {
    return Evaluation::createConstantZero(x);
  }

  /**
   * @brief 创建与参考对象具有相同导数维数的常量 1。
   *
   * @param[in] x 参考对象。
   * @return 常量 1 Evaluation。
   */
  [[nodiscard]] static Evaluation createConstantOne(const Evaluation &x)
  {
    return Evaluation::createConstantOne(x);
  }

  /**
   * @brief 创建常量 Evaluation。
   *
   * 对动态维数 Evaluation，该接口没有参考维数可供推断，因此底层实现会拒绝调用；
   * 动态模型应使用带 numDeriv 或参考对象的重载。
   *
   * @param[in] inputValue 常量函数值。
   * @return 常量 Evaluation。
   */
  [[nodiscard]] static Evaluation createConstant(const ValueType &inputValue)
  {
    return Evaluation::createConstant(inputValue);
  }

  /**
   * @brief 创建具有显式导数维数的常量 Evaluation。
   *
   * @param[in] numDeriv 导数数量。
   * @param[in] inputValue 常量函数值。
   * @return 常量 Evaluation。
   */
  [[nodiscard]] static Evaluation createConstant(unsigned numDeriv,
                                                 const ValueType &inputValue)
  {
    const int derivativeCount =
        DenseAd::detail::checkedDerivativeCount(numDeriv);

    if constexpr (numVars == DenseAd::DynamicSize)
    {
      return Evaluation::createConstant(derivativeCount, inputValue);
    }
    else
    {
      if (derivativeCount != numVars)
      {
        throw std::logic_error(
            "Requested derivative dimension does not match the fixed-size "
            "Evaluation type.");
      }

      // 固定维数由类型本身确定，直接使用无维数参数的工厂接口。
      // 固定维数由 Evaluation 类型本身决定，因此直接使用无维数参数的工厂接口。
      return Evaluation::createConstant(inputValue);
    }
  }

  /**
   * @brief 根据参考对象创建同维数常量 Evaluation。
   *
   * @param[in] x 参考 Evaluation。
   * @param[in] inputValue 常量函数值。
   * @return 与参考对象导数维数一致的常量 Evaluation。
   */
  [[nodiscard]] static Evaluation createConstant(
      const Evaluation &x, const ValueType &inputValue)
  {
    return Evaluation::createConstant(x, inputValue);
  }

  /**
   * @brief 创建固定维数独立变量。
   *
   * 动态维数 Evaluation 无法仅由该接口推断总导数数量，因此动态模型应使用另外两个重载。
   *
   * @param[in] inputValue 独立变量函数值。
   * @param[in] varIdx 导数索引。
   * @return 独立变量 Evaluation。
   */
  [[nodiscard]] static Evaluation createVariable(const ValueType &inputValue,
                                                 unsigned varIdx)
  {
    return Evaluation::createVariable(inputValue, static_cast<int>(varIdx));
  }

  /**
   * @brief 创建具有显式导数数量的独立变量。
   *
   * 该重载补齐 MathToolbox 通用 variable(numDeriv, value, idx) 接口，
   * 是动态维数 Evaluation 正确工作的必要入口。
   *
   * @param[in] numDeriv 总导数数量。
   * @param[in] inputValue 独立变量函数值。
   * @param[in] varIdx 当前独立变量对应的导数索引。
   * @return 独立变量 Evaluation。
   */
  [[nodiscard]] static Evaluation createVariable(unsigned numDeriv,
                                                 const ValueType &inputValue,
                                                 unsigned varIdx)
  {
    const int derivativeCount =
        DenseAd::detail::checkedDerivativeCount(numDeriv);
    const int derivativeIndex = static_cast<int>(varIdx);

    if constexpr (numVars == DenseAd::DynamicSize)
    {
      return Evaluation::createVariable(derivativeCount, inputValue,
                                        derivativeIndex);
    }
    else
    {
      if (derivativeCount != numVars)
      {
        throw std::logic_error(
            "Requested derivative dimension does not match the fixed-size "
            "Evaluation type.");
      }

      // 固定维数由类型本身确定，避免依赖旧手工特化中不一致的三参数构造路径。
      return Evaluation::createVariable(inputValue, derivativeIndex);
    }
  }

  /**
   * @brief 根据参考对象创建具有相同导数维数的独立变量。
   *
   * @param[in] x 参考 Evaluation。
   * @param[in] inputValue 独立变量函数值。
   * @param[in] varIdx 当前独立变量对应的导数索引。
   * @return 与参考对象导数维数一致的独立变量 Evaluation。
   */
  [[nodiscard]] static Evaluation createVariable(const Evaluation &x,
                                                 const ValueType &inputValue,
                                                 unsigned varIdx)
  {
    return Evaluation::createVariable(x, inputValue,
                                      static_cast<int>(varIdx));
  }

  /**
   * @brief 将 Evaluation 保持为相同 Evaluation 类型。
   *
   * @tparam LhsEval 目标类型。
   * @param[in] eval 输入对象。
   * @return 输入对象副本。
   */
  template <class LhsEval>
  [[nodiscard]] static std::enable_if_t<std::is_same_v<Evaluation, LhsEval>,
                                        LhsEval>
  decay(const Evaluation &eval)
  {
    return eval;
  }

  /**
   * @brief 将 Evaluation 的函数值衰减为普通浮点类型。
   *
   * 该转换显式丢弃外层导数信息，仅用于调用方明确需要普通数值的场景。
   *
   * @tparam LhsEval 目标普通浮点类型。
   * @param[in] eval 输入自动微分对象。
   * @return Evaluation 的 value 部分转换后的结果。
   */
  template <class LhsEval>
  [[nodiscard]] static std::enable_if_t<std::is_floating_point_v<LhsEval>,
                                        LhsEval>
  decay(const Evaluation &eval)
  {
    return static_cast<LhsEval>(eval.value());
  }

  /**
   * @brief 按给定容差比较两个 Evaluation 的函数值和全部导数。
   *
   * 动态维数对象的导数数量不一致时直接返回 false，而不是遗漏导数比较。
   *
   * @param[in] a 第一个对象。
   * @param[in] b 第二个对象。
   * @param[in] tolerance 比较容差。
   * @return 函数值和全部导数均在容差内时返回 true。
   */
  [[nodiscard]] static bool isSame(const Evaluation &a, const Evaluation &b,
                                   Scalar tolerance)
  {
    if (a.size() != b.size())
    {
      return false;
    }

    if (!InnerToolbox::isSame(a.value(), b.value(), tolerance))
    {
      return false;
    }

    for (int i = 0; i < a.size(); ++i)
    {
      if (!InnerToolbox::isSame(a.derivative(i), b.derivative(i), tolerance))
      {
        return false;
      }
    }

    return true;
  }

  /** @brief 返回两个输入中的较大值。 */
  template <class Arg1, class Arg2>
  [[nodiscard]] static Evaluation max(const Arg1 &arg1, const Arg2 &arg2)
  {
    return DenseAd::max(arg1, arg2);
  }

  /** @brief 返回两个输入中的较小值。 */
  template <class Arg1, class Arg2>
  [[nodiscard]] static Evaluation min(const Arg1 &arg1, const Arg2 &arg2)
  {
    return DenseAd::min(arg1, arg2);
  }

  /** @brief 计算绝对值。 */
  [[nodiscard]] static Evaluation abs(const Evaluation &arg)
  {
    return DenseAd::abs(arg);
  }

  /** @brief 计算正切函数。 */
  [[nodiscard]] static Evaluation tan(const Evaluation &arg)
  {
    return DenseAd::tan(arg);
  }

  /** @brief 计算反正切函数。 */
  [[nodiscard]] static Evaluation atan(const Evaluation &arg)
  {
    return DenseAd::atan(arg);
  }

  /** @brief 计算两个 Evaluation 的 atan2。 */
  [[nodiscard]] static Evaluation atan2(const Evaluation &arg1,
                                        const Evaluation &arg2)
  {
    return DenseAd::atan2(arg1, arg2);
  }

  /** @brief 计算 Evaluation 与普通值组合的 atan2。 */
  template <class Arg2>
  [[nodiscard]] static Evaluation atan2(const Evaluation &arg1,
                                        const Arg2 &arg2)
  {
    return DenseAd::atan2(arg1, arg2);
  }

  /** @brief 计算普通值与 Evaluation 组合的 atan2。 */
  template <class Arg1>
  [[nodiscard]] static Evaluation atan2(const Arg1 &arg1,
                                        const Evaluation &arg2)
  {
    return DenseAd::atan2(arg1, arg2);
  }

  /** @brief 计算正弦函数。 */
  [[nodiscard]] static Evaluation sin(const Evaluation &arg)
  {
    return DenseAd::sin(arg);
  }

  /** @brief 计算反正弦函数。 */
  [[nodiscard]] static Evaluation asin(const Evaluation &arg)
  {
    return DenseAd::asin(arg);
  }

  /** @brief 计算双曲正弦函数。 */
  [[nodiscard]] static Evaluation sinh(const Evaluation &arg)
  {
    return DenseAd::sinh(arg);
  }

  /** @brief 计算反双曲正弦函数。 */
  [[nodiscard]] static Evaluation asinh(const Evaluation &arg)
  {
    return DenseAd::asinh(arg);
  }

  /** @brief 计算余弦函数。 */
  [[nodiscard]] static Evaluation cos(const Evaluation &arg)
  {
    return DenseAd::cos(arg);
  }

  /** @brief 计算反余弦函数。 */
  [[nodiscard]] static Evaluation acos(const Evaluation &arg)
  {
    return DenseAd::acos(arg);
  }

  /** @brief 计算双曲余弦函数。 */
  [[nodiscard]] static Evaluation cosh(const Evaluation &arg)
  {
    return DenseAd::cosh(arg);
  }

  /** @brief 计算反双曲余弦函数。 */
  [[nodiscard]] static Evaluation acosh(const Evaluation &arg)
  {
    return DenseAd::acosh(arg);
  }

  /** @brief 计算平方根。 */
  [[nodiscard]] static Evaluation sqrt(const Evaluation &arg)
  {
    return DenseAd::sqrt(arg);
  }

  /** @brief 计算指数函数。 */
  [[nodiscard]] static Evaluation exp(const Evaluation &arg)
  {
    return DenseAd::exp(arg);
  }

  /** @brief 计算自然对数。 */
  [[nodiscard]] static Evaluation log(const Evaluation &arg)
  {
    return DenseAd::log(arg);
  }

  /** @brief 计算以 10 为底的对数。 */
  [[nodiscard]] static Evaluation log10(const Evaluation &arg)
  {
    return DenseAd::log10(arg);
  }

  /** @brief 计算 Evaluation 底数与普通指数的幂。 */
  template <class RhsValueType>
  [[nodiscard]] static Evaluation pow(const Evaluation &arg1,
                                      const RhsValueType &arg2)
  {
    return DenseAd::pow(arg1, arg2);
  }

  /** @brief 计算普通底数与 Evaluation 指数的幂。 */
  template <class LhsValueType>
  [[nodiscard]] static Evaluation pow(const LhsValueType &arg1,
                                      const Evaluation &arg2)
  {
    return DenseAd::pow(arg1, arg2);
  }

  /** @brief 计算两个 Evaluation 的幂。 */
  [[nodiscard]] static Evaluation pow(const Evaluation &arg1,
                                      const Evaluation &arg2)
  {
    return DenseAd::pow(arg1, arg2);
  }

  /**
   * @brief 判断 Evaluation 的函数值和全部导数是否均为有限数。
   *
   * @param[in] arg 输入自动微分对象。
   * @return 全部存储量均有限时返回 true。
   */
  [[nodiscard]] static bool isfinite(const Evaluation &arg)
  {
    if (!InnerToolbox::isfinite(arg.value()))
    {
      return false;
    }

    for (int i = 0; i < arg.size(); ++i)
    {
      if (!InnerToolbox::isfinite(arg.derivative(i)))
      {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief 判断 Evaluation 的函数值或任一导数是否为 NaN。
   *
   * @param[in] arg 输入自动微分对象。
   * @return 存在 NaN 时返回 true。
   */
  [[nodiscard]] static bool isnan(const Evaluation &arg)
  {
    if (InnerToolbox::isnan(arg.value()))
    {
      return true;
    }

    for (int i = 0; i < arg.size(); ++i)
    {
      if (InnerToolbox::isnan(arg.derivative(i)))
      {
        return true;
      }
    }

    return false;
  }
};

} // namespace MPMC
