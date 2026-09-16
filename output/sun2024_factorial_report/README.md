# Sun-2024超临界水-CO2一维砂管报告包

本目录只汇集本轮新生成算例的派生结果。绘图脚本不会搜索或读取仓库中已有的
数值实验结果。

## 交付物

- `report/超临界水_CO2_nC16_一维砂管可复现数值算例报告.docx`：可编辑报告；
- `report/超临界水_CO2_nC16_一维砂管可复现数值算例报告.pdf`：固定版式报告；
- `figures/*.png`：600 dpi科研图；
- `figures/*.pdf`：矢量科研图；
- `data/*.csv`：图件对应的机器可读数据；
- `data/figure_provenance.json`：数据源、变换、绘图环境和模型边界。

## 数值结果复现

在WSL/GNU Make环境中，从仓库根目录执行：

```bash
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh grid
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh time
RATE_BASIS=in_situ_volume \
  bash case/sun2024_scw_co2_nc16_factorial_1d/scripts/run_matrix.sh formal
```

正式矩阵默认使用`nx=96`、`dt_pv=0.01`和`target_pv=4`。严格物理复现应在取得
泵/中间容器的计量状态密度后改用`reference_density`，详见算例目录README。

## 图表复现

```powershell
python case/sun2024_scw_co2_nc16_factorial_1d/scripts/generate_report_artifacts.py `
  --formal-root tmp/sun2024_factorial_report_v1/formal `
  --case-root case/sun2024_scw_co2_nc16_factorial_1d `
  --grid-summary tmp/sun2024_factorial_convergence_v4/grid/convergence_summary.csv `
  --time-summary tmp/sun2024_factorial_convergence_v4/time/convergence_summary.csv `
  --output-root output/sun2024_factorial_report
```

图件不做平滑，采用色盲友好配色、线型和标记冗余编码，同时输出PNG和PDF。

可编辑报告由汇总数据与图件生成：

```powershell
python case/sun2024_scw_co2_nc16_factorial_1d/scripts/build_case_report.py `
  --artifact-root output/sun2024_factorial_report `
  --output "output/sun2024_factorial_report/report/超临界水_CO2_nC16_一维砂管可复现数值算例报告.docx" `
  --git-commit 2e4be8098077
```

该步骤需要`python-docx`。报告中的代码基线应替换为实际运行时的Git提交号。
