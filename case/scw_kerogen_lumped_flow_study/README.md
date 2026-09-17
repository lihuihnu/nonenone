# SCW–干酪根裂解产物拟组分流动算例集合

本目录用于后续建立“超临界水与预生成干酪根裂解产物的相平衡—黏度—多相多组分运移”数值实验集合。

## 研究边界

- 主目标工况固定为 **380 °C、25 MPa**，直接对齐 Zhao et al. (2023) 的 Chang 7 超临界水转化实验；后续保留 **360 °C、25 MPa** 作为同压力亚临界水对照候选。
- 研究对象是**预生成的干酪根裂解产物**与注入水之间的相平衡和流动耦合。
- 当前不把干酪根裂解反应、aquathermolysis、焦炭生成或能量方程纳入本算例集合；这些过程若后续实现，应单独扩展并重新验证。
- 组分拟组分化、PR/CPA 参数、二元作用参数和黏度模型必须有可追溯实验或文献依据，不以人为调参替代物性标定。

## 当前数据审计状态

主样品仍以 Zhao et al. (*Sustainable Energy & Fuels*, 2023, DOI `10.1039/D2SE01361D`) 的完整 Chang 7 raw-shale `380 °C / 25 MPa / 4 h` 实验为第一优先级数据源。

另已审计 Zhao et al. (*Industrial & Engineering Chemistry Research*, 2023, DOI `10.1021/acs.iecr.3c02759`) 的原文转写包。该论文研究的是同地质来源、但**酸洗去矿物后的 Type-II 纯干酪根**，因此只作为 `SECONDARY_PAIRED` 约束。

其中 Figure 6 的 380 °C pure-kerogen generated-oil simulated distillation 有直接印刷数值：

- Gasoline, IBP–180 °C: `0.81 wt%`
- Diesel, 180–350 °C: `23.73 wt%`
- Middle, 350–500 °C: `34.11 wt%`
- Heavy, >500 °C: `41.35 wt%`

这组数据可用于设计候选 boiling-range lumping 和代表组分筛选，但不能直接替代完整 raw-shale 主样品最终组成。完整数据与源数据冲突说明见 `fluid_characterization/raw/`。

## 计划的流体表征

第一阶段优先研究：

- H2O；
- Light：候选 C6–C14 lump；
- Middle：候选 C15–C20 lump；
- Heavy：候选 C21+ lump；
- 如实验数据表明气体产物不可忽略，再增加 CH4/C2–C5 gas lump。

上述碳数切分仍只是**候选方案**。ACS Figure 6 给出了另一套可追溯的 boiling-range 分段（IBP–180 / 180–350 / 350–500 / >500 °C），后续 lumping 应同时比较碳数切分和馏程切分对 EOS 参数化、平均分子量和黏度闭合的影响，再决定最终方案。

现有 `nC4 / nC10 / squalane` 体系可作为代码和机理筛选基线，但不直接等同于 380 °C 干酪根裂解产物的最终实验表征。

## 计划的数值实验

基础几何保持伪三维单层结构：

- 网格：`60 x 20 x 1`；
- 均质岩石；
- 等温全组分多相流；
- 左侧注入水、右侧生产；
- 主工况：`380 °C / 25 MPa`；
- 后续控制：`360 °C / 25 MPa`，用于同压力下比较亚临界水与超临界水；
- 重点比较温度跨临界变化，以及 PR 与 CPA 的热力学结构差异。

主要观测量：

1. 各相组成、密度和相态演化；
2. 水富相和烃富相黏度及相流度；
3. 生产井 H2O、Light、Middle、Heavy 的质量分数随 PVI 的变化；
4. 各拟组分累计采出率；
5. 轻/重选择性产出指标；
6. 组分质量守恒和数值收敛性。

## 后续目录规划

待物性调研完成后，再逐步加入：

- `fluid_characterization/`：380 °C 裂解产物原始数据、拟组分规则和参数证据；
- `binary_pvt/`：H2O–各拟组分的 PR/CPA 二元相平衡标定；
- `viscosity/`：IAPWS、f-theory、Pedersen 等黏度模型验证；
- `flow_360c_cpa/`：360 °C / 25 MPa 亚临界水对照；
- `flow_380c_cpa/`：380 °C / 25 MPa 主算例；
- `flow_380c_pr/`：380 °C / 25 MPa PR 结构对照。

上述子目录在对应数据、参数和验收标准明确后再创建，避免提前固化未经验证的模型假设。
