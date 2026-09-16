# H2O–squalane 60×20×1 均质二维算例

本算例把 `scw_kerogen_squalane_1d` 扩展为 1200 个均质网格，保持 1.20×0.10×0.10 m 外形、0.35 孔隙度、1500 mD 水平渗透率、总孔隙体积和 PVI 速率不变。注井在 `i=0,j=9`，采井在 `i=59,j=9`。由于 `ny=20` 为偶数，两个几何中心候选网格是 `j=9/10`，主算例相对短边中点偏移 2.5 mm。

运行时用 `-eos pr|sw|cpa` 切换模型。三者使用同一总体组成
`z(H2O/squalane)=0.60/0.40`，在 27.74 MPa、653.2 K 均从 Oil 角色单相
起步；不再人为给定初始水相。`zH2O=0.60` 受实验 squalane 富相端点
`xH2O=0.911` 约束，同时位于三模型共同单相区（CPA 在约 0.62 已预测两相）。
PR 非水相使用直接组成拟合的 `kij(653.2 K)=0.0532336595`；SW 的水相
BIP 仍是对 C30 的结构外推。Water 角色只有在非水组分不超过 2 mol% 时
才能调用 IAPWS 黏度，防止高含水但质量上烃富的单相被误标。

```bash
make case CASE=scw_kerogen_squalane_2d -j2
make run CASE=scw_kerogen_squalane_2d NP=2 EOS=pr \
  RESULT_DIR=./results/pr
```

产井见水以 `m_component_H2O` 计算的总产出质量分数为主判据；相体积含水率
仅保留为相槽诊断。当前烃相黏度仍是 LBC 对照路径，不能作为已完成的
f-theory/Pedersen 定量标定。
