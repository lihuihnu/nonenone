# H2O–CO2–CH4–nC16 three-EOS flash-plot example

这个 example 演示 `tools` 如何在**不进入 reservoir time loop** 的情况下，直接复用生产 `CubicThreePhaseFlash` 对同一套 H2O–CO2–CH4–nC16 流体分别运行：

- ordinary Peng–Robinson (PR)；
- Søreide–Whitson (SW)；
- Cubic-Plus-Association (CPA)。

## 1. 两种分析视角

### Unrestricted O/G/W

主 P–T 图和 P/T 路径图允许油富集、气富集、水富集三相自由出现，用于研究真实三相相态、phase fraction、saturation、Z、density、K-value 和 phase composition。

参考扫描覆盖：

- `T = 300–850 K`；
- `P = 1–600 bar`。

它跨过 ordinary water critical reference `Tc=647.096 K`, `Pc=220.64 bar`，并在两侧保留超临界余量。该参考点只用于画图标记，不会用 IAPWS 方程替换 PR/SW/CPA。

### Restricted O/G VLE projection

论文中常见的 bubble/dew envelope 是 conventional vapor–liquid 概念。为了不让第三个 water-rich phase 改变 bubble/dew 的定义，example 另外使用同一 EOS / 同一四组分流体，但将允许相集合限制为 O+G，输出：

- continuous bubble line；
- continuous dew line；
- envelope-coalescence critical estimate；
- composition-path critical-locus estimate。

这些 critical 输出明确是 **estimate**。CSV 同时给出 `critical_relative_envelope_width`、`critical_oil_gas_composition_L1` 和 `critical_estimate_score`，用于判断候选点离真正相合并有多远；它不是严格的 Michelsen criticality-equation solver。

## 2. 输出

每个 `pr/`, `sw/`, `cpa/` 目录生成：

- `pt_three_phase.csv` + `plot_pt_three_phase.m`
- `plot_pt_phase_count_quality.m`
- `pt_phase_onset_boundaries.csv` + `plot_pt_phase_onset_boundaries.m`
- `pt_vle_envelope.csv` + `plot_pt_vle_envelope.m`
- `pressure_profile_673K.csv` + `plot_pressure_profile_673K.m`
- `temperature_profile_300bar.csv` + `plot_temperature_profile_300bar.m`
- `critical_locus_estimate.csv` + `plot_critical_locus_estimate.m`

根目录另外生成：

- `plot_compare_vle_envelopes.m`
- `plot_compare_critical_loci.m`

profile CSV 可直接绘制：

- `beta_o/beta_g/beta_w`
- `So/Sg/Sw`
- `rho_m_o/rho_m_g/rho_m_w`
- `Z_o/Z_g/Z_w`
- `Kg_<component>` / `Kw_<component>`
- `xo/yg/xw` phase compositions

## 3. 构建与运行

```bash
cd tools
make examples
make run-three-eos
```

重建随工程跟踪的参考输出：

```bash
make reference-three-eos
```

MATLAB 中例如：

```matlab
cd('.../tools/example/h2o_co2_ch4_nc16_three_eos/results')
run('plot_compare_vle_envelopes.m')
run('plot_compare_critical_loci.m')
run('cpa/plot_pt_three_phase.m')
run('cpa/plot_pressure_profile_673K.m')
```

## 4. 参数与适用域说明

本例的纯组分参数/BIP 沿用项目 3p4c H2O–CO2–CH4–nC16 benchmark 风格。SW 对 CH4/CO2 使用生产 correlation；原始 SW hydrocarbon fit 主要针对轻烃，因此本例**不把该 correlation 外推到 nC16**，而对 nC16–water 采用显式 `0.5` engineering BIP。CPA 使用项目现有 Standard 4C water + non-associating SRK fallback；这不是重新拟合的 H2O–CO2–CH4–nC16 CPA 参数集。

因此 `300–850 K, 1–600 bar` 表示**数值扫描范围**，不表示三个 backend 在该整个区域都拥有同等级的实验拟合精度。失败 flash 保留在 CSV 中（`phase_code=0`, `status_code!=0`），不会被 MATLAB 静默删除成“看似连续”的边界。

参考来源：

- IAPWS-95 (2018 revised release), ordinary-water critical/reference properties and water EOS validity range.
- Søreide & Whitson, *Fluid Phase Equilibria* 77 (1992) 217–240, DOI `10.1016/0378-3812(92)85105-H`.
- Kontogeorgis et al., *Industrial & Engineering Chemistry Research* 35 (1996) 4310–4318, DOI `10.1021/ie9600203`.
