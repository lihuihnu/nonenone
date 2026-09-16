# Tools/example 子模块说明书

## 1. 职责

`tools/example/` 只存放**具体研究/应用代码**。它类似 `case/` 中的算例，但不启动 reservoir time loop；每个 example 应直接使用 `tools` 公共 API，并将自己的组分参数、扫描范围和输出文件组织在独立目录中。

## 2. 当前示例

`pr_c2_c3_nc5/`：C2/C3/nC5 ordinary Peng-Robinson 示例。

生成：

1. 固定 `z=(0.3,0.4,0.3)` 的 P-T 相区图数据；
2. 320 K 下 C2-rich 到 nC5-rich 组成路径的 P-composition 数据；
3. 320 K、8 bar 下 C2/C3/nC5 三元组成图数据与 O-G tie lines。

该示例显式使用 restricted O+G flash，因此它验证的是 conventional hydrocarbon VLE，而不是含水三相问题。

`h2o_co2_ch4_nc16_three_eos/`：H2O/CO2/CH4/nC16 三 backend 示例。它同时使用：

- unrestricted O/G/W：真实三相 P-T map、oil/gas/water onset、673.15 K pressure profile、300 bar temperature profile；
- restricted O/G：standard bubble/dew envelope 与 composition-path critical-locus estimate；
- PR / Søreide-Whitson / CPA：同一组分和扫描坐标下分别输出，并生成 backend-comparison MATLAB scripts。

`three_eos_public_benchmark/`：不复制生产公式的独立验证入口。C++ 批处理端只调用生产 EOS/Flash；Python 驱动负责 NIST ThermoML 解析、ThermoPack/NeqSim oracle、冻结门禁、逐点 CSV 和独立工程图。`reference_data/` 保留带哈希的公开原始输入，`results/` 仍为不跟踪的可重生成输出。

`scw_binary_calibration/`：不进入流动时间循环的 H2O–squalane 二元标定入口。它顺序回归纯角鲨烷体积平移、LBC 有效临界体积和生产 flash 的温度相关 BIP，并以 IAPWS-2008 检查纯水黏度。拟合与留出误差分开输出；未通过数据不确定度/组成精度门禁的参数不得自动写回流动算例。

同目录的 `run_phase_diagrams.py` 使用 May CH4-C2H6、Messabeb CO2-NaCl-H2O 与 Kurihara 甲醇-水公开数据分别重建 PR、SW、CPA 相图，并调用 `h2o_co2_ch4_nc16_three_eos` 生成同流体 P-T 包络对比。六类图分别导出单坐标轴英文版和中文版，中文版坐标文字为中文、图例为英文；公开数据验证与工程参数敏感性比较在报告中明确分开。`--redraw-only` 复用保存的逐点生产 CSV，只重建图形与报告，不重新执行热力学计算。

参考 P-T 范围为 300–850 K、1–600 bar，跨越 IAPWS ordinary-water critical reference；失败点原样保留，不做跨 failure 插值。

## 3. 新增具体应用

建议结构：

```text
tools/example/my_study/
  main.cpp
  README.md
  results/            # 运行生成，Git 忽略
  reference_output/   # 可选本地重生成目录，Git 忽略
```

通用算法不要放入 example；发现可复用逻辑后应上移到 `tools/include/tools/` 并增加 unit test。

## 4. 不可破坏约束

- example 参数必须写明来源/用途。
- `results/` 永不进入 Git。
- `reference_output/` 与 `results/` 均为可重生成运行产物，不进入 Git；真正的外部 reference input 放在 `reference_data/`。
- 具体应用不得复制 `CubicThreePhaseFlash` 或 EOS 公式。
