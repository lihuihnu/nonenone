/**
 * @file natural_runtime_options.hpp
 * @brief Natural/PETSc 运行时数值选项和输出控制参数。
 */
#pragma once

#include <natural/petsc/natural_scaling.hpp>

#include <array>
#include <optional>

namespace MPMC
{

/**
 * @brief Natural PETSc 运行时参数，与具体网格类型无关。
 *
 * 网格只负责 DOF、ghost、几何与岩石属性访问；时间离散、EOS/flash、
 * Land、吸附、溶解等模型参数统一由 Natural runtime 管理。
 */
template <class Indices>
struct NaturalRuntimeOptions final
{
    /** 当前隐式时间步长。Natural 残差中的积累项使用 `(M^{n+1}-M^n)/timeStep`。 */
    double timeStep{1.0};

    /** 逸度平衡残差的无量纲缩放因子，只改变方程尺度，不改变热力学平衡条件。 */
    double fugacityScalingFactor{1.0};

    /** Land 滞留模型常数 `C_L`，满足 `S_gr=S_g,max/(1+C_L S_g,max)`。 */
    double landConstant{2.0};

    /** 岩石骨架密度 [kg/m^3]，用于吸附组分质量 `m_ads=(1-phi) rho_r theta V_b`。 */
    double rockDensity{2650.0};

    /** 水相溶解 CO2 对应的组分索引；-1 表示使用 FluidSystem::aqueousCO2Component。 */
    int dissolvedCO2Component{-1};

    /** 是否由 PETSc SNESVI 对主变量施加运行时上下界。 */
    bool useVariableBounds{false};

    /**
     * @brief 全组分水相 Newton 候选态允许的非 H2O/CO2 总摩尔分数。
     *
     * 未设置时不增加模型相关限制。IAPWS+Garcia/McBride-Wright 这类仅适用于
     * H2O+CO2 富水相的闭包可显式设置该值，使 line-search 候选态保持在闭包
     * 的适用域内，而不改变收敛解或放宽物性模型的适用范围。
     */
    std::optional<double> maximumAqueousUnsupportedNewtonMoleFraction;

    /** 上述水相 Newton 适用域中 H2O 与 CO2 的组分索引。 */
    int aqueousNewtonWaterComponent{-1};
    int aqueousNewtonCo2Component{-1};

    /** 一致的残差行/变量列尺度化；只改变数值尺度，不改变非线性零点。 */
    NaturalScalingOptions scaling{};

    /**
     * @brief 各组分标准态气体密度 [kg/m^3]。
     *
     * 仅在需要由储层摩尔/质量流量换算标准态体积流量时使用；未设置时
     * runtime 采用模型内可获得的默认标准态密度。
     */
    std::optional<std::array<double, Indices::numComponents>> standardGasDensity;
};

} // namespace MPMC
