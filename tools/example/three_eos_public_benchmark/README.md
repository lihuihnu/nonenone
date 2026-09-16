# 三状态方程与 Flash 独立公开基准

该工具从当前生产 `CubicEquationOfState` 与 `CubicThreePhaseFlash` API 直接取结果，不包含另一份 MPMC EOS/Flash 公式。验证证据分为三层：

1. **公开方程/同参数恒等性**：PR、SRK-CPA 与 ThermoPack 2.2.3 比较 `Z`、摩尔密度和 `ln(phi_i)`；SW 水 alpha、1992 原始 CO2-H2O 水相 BIP 与 Chabab 2019 改进 BIP 分别和公开方程的独立转录比较；
2. **独立 Flash**：PR/CPA 与 ThermoPack、SW 与 NeqSim 3.18.0 比较相数、相分率和相组成，并独立重算 MPMC 物料闭合；
3. **公开实验预测**：PR 对 May et al. 甲烷–乙烷 VLE，SW 对 Messabeb et al. CO2–NaCl–H2O 溶解度，CPA 对 Soujanya et al. 甲醇–水 VLE。

本程序不使用报告中的实验点现场拟合参数。SW 同时展示 1992 原始参数化和 Chabab 2019 公开再拟合参数；后者并非由本程序生成，因此修复前后比较可复现，但不把该实验层冒充完全独立的留出验证。实验误差不混入程序实现门禁。

## 冻结协议

- Latin hypercube 固定种子：`20260822`；
- 方程门禁：`|ΔZ| <= 1e-9`、`max|Δln(phi)| <= 1e-7`、`|Δrho|/rho <= 1e-8`；
- 同参数 Flash：相数一致、`max|Δbeta| <= 1e-6`、`max|Δx| <= 1e-6`；
- MPMC 物料闭合：`max_i |z_i-sum_p beta_p x_pi| <= 1e-10`；
- SW 公开关联式：`max|Δcorrelation| <= 1e-12`；
- SW 与 NeqSim 的完整 Flash 比较为同模型家族的独立实现证据，不宣称两套软件的隐藏数据库参数逐位相同，因此不进入严格 Flash 门禁。

这些标准在运行结果产生前写入代码和文档，不允许根据结果放宽。

空间填充矩阵覆盖：PR 方程 14 个状态×液/汽根、PR Flash 18 个状态；CPA 方程 12 个状态×液/汽根、CPA Flash 16 个状态；SW Flash 12 个状态；SW 公开关联式 16 个状态×2 个公式。压力采用对数分层，温度、组成和盐度采用均匀分层。逐点输入与结果都写入 CSV，失败不会被静默删除。

## 同参数锁定

- PR：C1/C2 的 ThermoPack 组分数据库参数、`k12=0`，原 PR 全精度 `Ωa=0.4572355289213822`、`Ωb=0.07779607390388846`；
- CPA：H2O/MEOH 的 sCPA 4C/2B 参数、`k12=-0.09`，其中 Classic alpha 使用 ThermoPack 数据库的 `Tc(H2O)=647.3 K`，而非另一套参考物性常数；
- SW：兼容路径保留 Søreide–Whitson 1992 原关联式；CO2-NaCl-H2O 实验路径明确选择 Chabab et al. 2019 改进 aqueous CO2-H2O BIP。NeqSim 同样显式选择 `CHABAB_2019`，但内部数据库和盐度质量基准不能全部逐项锁定，所以完整 Flash 结果仍只作为独立实现和趋势证据。

预运行发现水–甲醇液相的阻尼位点方程可能超过历史默认 100 次迭代。生产默认安全上限已提高到 2000，并加入单元回归；已收敛状态仍提前退出，公式、容差和冻结门槛没有变化。最终基准不设置专属迭代覆盖，直接验证生产默认值。

## 数据与来源

原始 JSON 位于 `reference_data/`，由 `fetch_reference_data.py` 下载并校验 SHA-256：

- [May et al. 2015, reference-quality methane binary VLE](https://trc.nist.gov/ThermoML/10.1021/acs.jced.5b00610.html)；
- [Messabeb et al. 2016, CO2 solubility in aqueous NaCl](https://trc.nist.gov/ThermoML/10.1021/acs.jced.6b00505.html)；
- [Soujanya et al. 2010, methanol-water VLE](https://trc.nist.gov/ThermoML/10.1016/j.jct.2009.11.020.html)。

理论与算法来源：

- Peng & Robinson, 1976, DOI `10.1021/i160057a011`；
- Soreide & Whitson, 1992, DOI `10.1016/0378-3812(92)85105-H`；
- Chabab et al., 2019, DOI `10.1016/j.ijggc.2019.102825`；
- Kontogeorgis et al., 1996, DOI `10.1021/ie9600203`；
- Michelsen, 1982, DOI `10.1016/0378-3812(82)85001-2` 与 `10.1016/0378-3812(82)85002-4`；
- [ThermoPack](https://github.com/thermotools/thermopack)，Apache-2.0；
- [NeqSim](https://github.com/equinor/neqsim)，Apache-2.0。

CPA 参数可直接追溯到 ThermoPack 官方 [Water.json](https://github.com/thermotools/thermopack/blob/main/fluids/Water.json) 与 [Methanol.json](https://github.com/thermotools/thermopack/blob/main/fluids/Methanol.json)。

## 运行

Oracle Python 环境与项目运行依赖隔离：

```bash
python -m venv .oracle-venv
.oracle-venv/bin/python -m pip install -r tools/example/three_eos_public_benchmark/requirements-oracle.txt
make -C tools bin/three_eos_public_benchmark
JAVA_HOME=/path/to/jdk17 .oracle-venv/bin/python \
  tools/example/three_eos_public_benchmark/run_benchmark.py \
  --mpmc tools/bin/three_eos_public_benchmark
```

生成物位于 `results/` 且不进入 Git：逐点 CSV、`summary.json`、`benchmark_report.md`、五张独立 PNG 和汇总图册 PDF。真正的公开原始输入保留在 `reference_data/`。

## 公开相图重建

`run_phase_diagrams.py` 另行生成三种 EOS 的公开相平衡图，不读取历史图片或报告：

- PR：May et al. 2015 的 CH4-C2H6 243.60 K 等温 P-x-y 图；
- SW-2019：Messabeb et al. 2016 的 CO2-NaCl-H2O 平衡溶解度图，分别绘制 323.15、373.15、423.15 K；
- SRK-CPA：Kurihara et al. 1995 表 2 的甲醇-水 333.15 K 等温 P-x-y 图；
- 同一 H2O-CO2-CH4-nC16 工程参数集的 PR/SW/CPA restricted O/G P-T 包络，用于模型敏感性比较，不作为公开实验验证。

Kurihara 数据保存在 `reference_data/kurihara1995_cpa_meoh_water_333K.json`，包含原论文 DOI、表格位置和 12 个公开 P-x-y 点。PR/SW 继续读取经 SHA-256 固定的 NIST ThermoML 原始 JSON。曲线来自生产 EOS/Flash 的逐点结果；失败点保留在 CSV 并在图中断开，不做跨失败点插值。

```bash
make -C tools \
  PYTHON=../.oracle-venv/bin/python \
  run-public-phase-diagrams
```

仅调整图形呈现而不重复 EOS/Flash 计算时，可复用已保存的逐点生产结果：

```bash
make -C tools \
  PYTHON=../.oracle-venv/bin/python \
  redraw-public-phase-diagrams
```

输出位于 `results/phase_diagrams/`：六类图各含英文版和中文版，共十二张独立 PNG、对应矢量 PDF、英文/中文图册、逐点请求/结果 CSV、指标与来源报告。报告中的数据源论文采用 GB/T 7714—2015 顺序编码格式，并将论文引用与机器可读数据地址分开列示。每张图仅含一个坐标轴，不设大标题和网格，采用完整四边框及框内纯数据图例；中文版坐标文字为中文，图例仍为英文，并以 `Our` 标识自主程序。生成物不进入 Git。

## 已知 Oracle 包装问题

NeqSim 3.18.0 的 `SystemSoreideWhitson` 已创建组分但未把组分计数传播到替换后的 phase。驱动调用 `setNumberOfComponents(2)` 只暴露这些已创建组分，不修改 JAR、参数或方程。ThermoPack 2.2.3 的 sCPA Python convenience wrapper 将布尔参数构造成 `c_int`；驱动以正确的 `c_bool` 调用同一公开库函数。两项补救均在报告中显式记录。
