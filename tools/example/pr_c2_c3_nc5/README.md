# PR C2/C3/nC5 phase-diagram example

这是 `tools` 的第一套实例。纯组分采用 C2、C3、nC5，`kij=0`。

公开回归点：

```text
T = 320 K
P = 8 bar
z = (0.3, 0.4, 0.3)
```

生产 flash 应得到 O+G，两相摩尔分数约为：

```text
beta_o = 0.2904665759
beta_g = 0.7095334241
```

## 运行

```bash
cd tools
make run-example
```

生成的参考内容：

```text
pt_phase_map.csv           + plot_pt_phase_map.m
pt_envelope.csv            + plot_pt_envelope.m
px_phase_map.csv           + plot_px_phase_map.m
px_envelope.csv            + plot_px_envelope.m
ternary_phase_map.csv      + plot_ternary_phase_map.m
critical_locus.csv         + plot_critical_locus.m
```

`reference_output/` 不随工程跟踪；需要固定参考输出时执行 `make reference` 即可确定性重建。

## MATLAB

```matlab
cd('.../tools/example/pr_c2_c3_nc5/reference_output')
run('plot_pt_phase_map.m')
run('plot_pt_envelope.m')
run('plot_px_phase_map.m')
run('plot_px_envelope.m')
run('plot_ternary_phase_map.m')
run('plot_critical_locus.m')
```

说明：

- `plot_pt_phase_map.m`：P-T 相区散点图。
- `plot_pt_envelope.m`：连续泡点线、露点线和近似临界点。
- `plot_px_phase_map.m`：固定 320 K 的 P-composition 相区图。
- `plot_px_envelope.m`：固定 320 K 的连续 bubble/dew 线。
- `plot_ternary_phase_map.m`：三元相区图 + 代表性 O-G tie lines。
- `plot_critical_locus.m`：沿 C2↔nC5 组成路径的近似临界点曲线。

P-T 与 P-composition 的 envelope 脚本不重新做热力学计算，只读取 CSV 中已经提取好的边界。`critical_locus.csv` 中的临界点曲线是基于固定组成 P-T 包络的最小 bubble/dew 间距近似得到的趋势图，不是独立严格的临界 continuation 结果。
