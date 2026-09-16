# H2O + 轻/中/重裂解产物 60×20×1 均质二维算例

本算例与 H2O–squalane 二维算例共用网格、岩石、井位、总孔隙体积和 PVI 速率，仅把烃组分扩展为 nC4/nC10/squalane。运行时可用 `-eos pr|sw|cpa` 切换 PR、Søreide–Whitson 或 CPA。

三 EOS 使用完全相同的总体组成
`H2O/nC4/nC10/squalane=0.25/0.078125/0.1953125/0.4765625`，初态均为
Oil 角色单相且 `Sw=0`。653.2 K 的非水相参数为 H2O–nC4 `0.5091`、
H2O–nC10 `0.2618373654`、H2O–squalane `0.0532336595`，烃–烃 BIP 为零。
nC4 是公开先验，nC10 是 573.2–613.2 K 文献参数向目标温度的外推，
squalane 是目标窗组成拟合且仍有模型结构失配；三者证据等级不能混同。

```bash
make case CASE=scw_kerogen_lmh_2d -j2
make run CASE=scw_kerogen_lmh_2d NP=2 EOS=pr \
  RESULT_DIR=./results/pr
```

水相黏度使用 IAPWS-2008 工业项，产井见水用 H2O 组分质量分数。当前烃相
仍走 LBC 对照闭包；在获得能约束 637–653 K 外推的纯烃/混合物数据之前，
不把未标定 f-theory 或 Pedersen 伪装成主模型。
