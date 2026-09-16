# CO2–DME–H2O 公开 CPA 三相闪蒸算例

本算例使用公开的 `CO2 + dimethyl ether (DME) + H2O` 汽–液–液平衡
（VLLE）数据，检查生产 `CubicEquationOfState` 与
`CubicThreePhaseFlash` 在带交叉缔合的 CPA 参数下能否恢复三个平衡相。

## 公开依据

- 实验数据：Laursen 等在 308.15 K 测得的 CO2–DME–H2O VLLE 组成。
- CPA 参数：Folas 的 sCPA 参数、二元相互作用参数，以及 water–DME
  的 mCR-1 交叉缔合规则。

标准引用：

1. LAURSEN T, RASMUSSEN P, ANDERSEN S I. VLE and VLLE measurements of
   dimethyl ether containing systems[J]. Journal of Chemical & Engineering
   Data, 2002, 47(2): 198-202. DOI: [10.1021/je010154+](https://doi.org/10.1021/je010154%2B).
   [包含原始数据表的 DTU 公开学位论文](https://backend.orbit.dtu.dk/ws/portalfiles/portal/5482520/KT2002-Torben%2BLaursen-Measurments%2Band%2BModeling%2Bof%2BVLLE%2Bat%2BElevated%2BPressures.pdf).
2. FOLAS G K. Modeling of complex mixtures containing hydrogen bonding
   molecules[D]. Lyngby: Technical University of Denmark, 2007. ISBN
   978-87-91435-45-5. [DTU 公开全文](https://backend.orbit.dtu.dk/ws/files/123926886/CDocuments_and_SettingsalbDesktopPh.D._Georgios_K._Folas.pdf.pdf).

## 条件与参数

组分内部顺序为 `H2O/DME/CO2`，温度固定为 308.15 K，压力为 19.0、
31.7、46.0 和 58.8 bar。这里只选择三相中三个组分均非零的公开点。

| 组分 | a0 (Pa·m^6/mol^2) | b (m^3/mol) | c1 | donor | acceptor | epsilon (J/mol) | beta |
|---|---:|---:|---:|---:|---:|---:|---:|
| H2O | 0.12277 | 1.4515e-5 | 0.67359 | 2 | 2 | 16655 | 0.0692 |
| DME | 0.84354 | 4.9600e-5 | 0.72125 | 0 | 1 | — | — |
| CO2 | SRK 临界参数式 | SRK 临界参数式 | SRK 偏心因子式 | 0 | 0 | — | — |

二元参数为 `k(H2O,DME)=-0.160`、`k(H2O,CO2)=-0.066`、
`k(DME,CO2)=-0.016`。water–DME 交叉缔合使用
`epsilon_cross=8327.5 J/mol`、`beta_cross=0.2877`；径向分布函数使用
Folas 文中的 simplified CPA 形式。

原实验没有报告进料组成。算例在每个温压点取三个实验相组成的等权凸组合
作为总体组成，使总体组成严格位于实验三相三角形内部；该选择决定相分率，但在
给定温压下不应改变平衡相组成。

## 运行

```bash
cd tools
make run-cpa-dme-vlle
make plot-cpa-dme-vlle
make plot-cpa-dme-vlle-by-pressure
```

绘图命令需要 Matplotlib；若默认 `python3` 未安装该包，可通过
`PYTHON=/path/to/python-with-matplotlib` 指定已有 Python 环境。

生成文件位于 `results/`：

- `flash_results.csv`：收敛、相态、相分率和数值闭合；
- `phase_comparison.csv`：每个相、每个组分的实验值与计算值；
- `aard_by_phase_component.csv`：按相和组分统计的 AARD；
- `metrics.csv`：总体组成 MAE 与最大绝对误差。
- `figures/en/` 和 `figures/zh/`：英文、中文科研风格组成 parity 图的 PNG/PDF；
- `01_cpa_vlle_composition_parity`：全部压力合并的单图版本；
- `02_cpa_vlle_composition_by_pressure`：四个压力分别成图、采用统一坐标尺度的
  `2×2` 对比版本；
- `figures/figure_manifest.json`：合并图的数据哈希、编码方式、尺寸和导出参数；
- `figures/02_cpa_vlle_composition_by_pressure_manifest.json`：四压力分图的
  数据哈希、分图顺序和导出参数。

`reference_data/` 同时保留逐点实验组成和 Folas 表 7.5 的公开 CPA AAD，便于
复核本报告中的两类比较。

完整计算结果和适用边界见 [RESULT_REPORT.md](RESULT_REPORT.md)。
