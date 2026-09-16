# MPMC_SCW Tools

`tools/` 是不进入 reservoir solver 主循环的离线扩展层。它只消费已经验证的生产 API，适合相图/PVT/批量闪蒸/参数研究/格式转换等功能；production modules 不得反向依赖 `tools`。

当前 phase-diagram / flash-plot 能力包括：

- P–T / P–composition / ternary 相区扫描；
- continuous bubble/dew envelope；
- composition-path critical-locus **estimate**；
- unrestricted O/G/W oil/gas/water onset boundaries；
- phase count + gas-quality contours；
- pressure/temperature flash profiles：phase fraction、saturation、molar density、Z、K-value、phase composition；
- PR / Søreide–Whitson / CPA 对比 MATLAB 脚本；
- MATLAB-friendly CSV writer。

## 快速开始

```bash
cd tools
make examples
```

普通 PR 烃类示例：

```bash
make run-example
```

H2O–CO2–CH4–nC16 三 EOS / 超临界水范围示例：

```bash
make run-three-eos
```

公开 CO2–DME–H2O CPA 汽–液–液三相闪蒸示例：

```bash
make run-cpa-dme-vlle
```

该算例的公开参数、实验条件和计算结果见
[cpa_dme_co2_water_vlle/README.md](example/cpa_dme_co2_water_vlle/README.md)。

需要生成本地参考输出时，目录分别为：

```text
example/pr_c2_c3_nc5/reference_output/
example/h2o_co2_ch4_nc16_three_eos/reference_output/
```

这些目录不进入 Git；重新生成：

```bash
make reference
make reference-three-eos
```

完整 API、critical-estimate 边界和扩展规则见 [MODULE.md](MODULE.md)。
