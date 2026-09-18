# 17 — SCW-H01 初版完整算例设计与调研索引

调研基线：54995b96b022fc7cbd3aa30e6fcaf1f758cfe7ec；2026-09-18。
状态：DESIGN_READY。本文没有新增流动计算结果，没有修改EOS/Heavy参数、原门禁或CI。它是实验协议，不是已接入现有解析器的运行输入卡。优先级服从00_SCW_EXPERIMENT_REFOCUS_20260918.md：暂停CPA扩展，不以Jia Figure7完成作为本案前置条件。

## 1. 科学问题和最小体系

在相同温度、初始Heavy库存、实际PVI与井控下，水–Heavy组分交换及组成相关黏度反馈分别怎样改变重组分采出、携带路径和压降？允许无改善或不利影响，不要求CPA>SW>PR，不预设混相或三相。

体系为H2O+OIL_HEAVY。Heavy沿用>500°C裂解产物馏分定义，不是全产油，不是squalane。第一阶段只研究重质组分转移和采出；不能验证轻重选择性。产物预先生成，等温、无反应；这是假设，不是本样品在380°C必然热稳定的证据。

## 2. 调研按实验用途归档

以下相对路径位于本目录，除特别标注外。

|用途|文件|可用结论|限制|
|---|---|---|---|
|主样品/产物|01_RESEARCH_OBJECT.md、02_RAW_PRODUCT_DATASET.md、research_object_baseline.csv、fluid_characterization/raw|原页岩主样品、380°C/25MPa/4h产油及SARA|反应产率不是驱替采收率；原产油SARA不等于独立Heavy馏分SARA|
|馏分表征|03_EXPERIMENT_DRIVEN_LUMPING.md、04_PSEUDOCOMPONENT_CHARACTERIZATION.md、fluid_characterization/pseudo_component_characterization_380c.csv|实验馏程、图5内部权重和拟临界参数|酸洗干酪根是SECONDARY_PAIRED；相关式参数不是实测临界性质|
|相平衡|05、06、16；binary_pr_calibration/source_manifest.csv、h2o_heavy_evidence_manifest.csv|可用代理数据、适用域、缺口|不能把代理kij当真实Heavy参数|
|密度黏度|08_DENSITY_VISCOSITY_VALIDATION.md、property_validation|纯水实现验证、烃类参考点|混合相黏度未闭合；常温产出油黏度不是孔隙内高温黏度|
|介质/井控|11–14；porous_media|均质、Corey敏感性、PVI、注入率+最大BHP、出口固定BHP|旧1.2m、phi=.35、1500mD、.25PV/day不再当实测输入|
|产出|15_PRODUCER_COMPOSITION_AND_SELECTIVITY.md、producer_composition_output_contract.csv|真实守恒组分通量和accepted-step积分|出口单元组成不能代替产出通量；携带相分解需要单独核验|
|历史程序|case/scw_kerogen_common、scw_kerogen_lmh_1d、scw_kerogen_squalane_1d|构建/流动/后处理骨架|历史nC4/nC10/squalane算例不能重命名为本案|
|CPA支线|tools/example/cpa_*|保留已有失败与回归|不继续扩大为本案主任务|

最相关文献按作用排列：

- Zhao et al. Energy & Fuels 2018，DOI 10.1021/acs.energyfuels.7b03839：填砂驱替、采出和黏度响应。出版社摘要代表条件400°C/25MPa；仅指导实验，不作为380°C/28MPa同样品无反应定量验证。
- Matsui et al. JPI 2014，DOI 10.1627/jpi.57.118：出版社摘要范围603–653K、1–17MPa，指导分配与停留时间实验。其约30min内裂解可忽略的样品特定观察，不能直接移作本案结论。仓库另一清单区分定量/视觉范围，正式摄取应回原表逐行核对。
- Sato et al. JPI 2018，DOI 10.1627/jpi.61.256：仓库登记603–643K、2–10.2MPa常压渣油研究；不覆盖28MPa。
- RSC 2023，DOI 10.1039/D2SE01361D：原页岩身份与产物约束；IECR 2023，DOI 10.1021/acs.iecr.3c02759：酸洗干酪根四馏分先验，两样品不得混为同一实测。
- Stevenson 1994，DOI 10.1016/0378-3812(94)87016-0：dodecane/squalane代理体系，不转移为真实Heavy。
- Amani/Jia：仅高压bitumen背景与实现诊断，不作为主交付。

本次外部核对的是前两篇出版社摘要；未声称逐页重审全部文献。

## 3. 冻结物质身份，标明数据层级

|组分|馏程°C|四馏分产油质量占比%|MW g/mol|Tc K|Pc MPa|omega|
|---|---|---:|---:|---:|---:|---:|
|Gasoline|IBP–180|0.81|82.392|515.231|3.192179|0.269328|
|Diesel|180–350|23.73|215.583|736.423|1.751155|0.583473|
|Middle|350–500|34.11|387.349|870.437|1.138075|0.896554|
|Heavy|>500|41.35|660.132|982.864|0.808714|1.215676|

Heavy还有Tb=563.696°C、SG60/60=.942605、Vc=1793.890cm3/mol。馏分比例有实验来源，拟组分性质是provisional估计，Heavy状态为PROVISIONAL_HEAVY_TAIL_LOWER_BOUND_LIKE。

第一阶段烃基准重新归一化为100%Heavy，不能只填41.35%而把其余库存留空。初始纯Heavy有助于保证同温A/B/C库存一致；它是数值控制，不代表现场初始含水。

## 4. 数值主案

|项目|v1选择|
|---|---|
|名称|SCW_H01_HEAVY_DISPLACEMENT|
|空间|均质一维，200×1×1；网格复核100/400；xi=x/L∈[0,1]|
|温度|653.15K；全域和入口等温|
|初始/出口压力|28MPa；出口fixed BHP，入口压力是结果，不强称平均压力严格28MPa|
|初态|zH2O=0、zHeavy=1、SH=1；确认模型返回允许的稳定流体状态|
|注入|纯水；ReservoirTotalRate+maximumBhp|
|数值压力上限|30MPa仅为模型范围约束，不是硬件耐压许可|
|速率尺度|Delta p_ref=.10MPa；Qref=k A Delta p_ref/(mu_H_ref L)；Qtarget/Qref=1|
|时间|tau=t Qref/PV0；0–2实际PVI；1和2PVI主比较|
|步长|初始Delta tau=1e-4，最大1e-3；减半复核；建议单步Delta S<=.02|
|输出|每.02PVI；剖面0/.1/.25/.5/1/2PVI|
|重力/毛管/物理弥散|首轮0，均为显式机制假设|
|Corey基线|Swr=Sor=0，nw=no=2，端点1；不是实测SCW曲线|
|反应/吸附/焦炭|关闭|
|物性提供者|首审既有PR76；目标有效域内可追溯表格可替代；不做EOS排名|

无实测几何时使用孔隙体积归一化，每单元PV0/Nx。归一化Darcy式为q_alpha/Qref=−(mu_H_ref kr_alpha/mu_alpha) dp*/dxi，p*=(p−pout)/Delta p_ref。它不产生m、mD或分钟，实际L/A/k/PV和参考黏度到位后再换算Q及停留时间。Qref在A/B/C之间不重新选择。

入口切换BHP后记录实际注入量；不能设产出体积等于注入体积。整个25–30MPa旧0D范围只是数值覆盖，不是该范围内真实Heavy全部实测验证。

## 5. 对照定义

A_IMMISCIBLE：关闭交换，两纯流体；黏度为同T/P纯流体参考。它是数值反事实，不通过极端kij制造不混溶，不要求跨相逸度相等。

B_TRANSFER：开启局部相平衡，冻结组成黏度反馈；保持与A相同参考黏度。B−A包含相分率、密度等变化，不能命名为纯抽提贡献。

C_COUPLED：与B相同平衡模型，启用mu(T,P,x)。C−B为指定协议下的组成黏度反馈；不是独立线性贡献。C缺混合相黏度证据时只能给显式条件性闭合，不能人为设定有利降黏倍数作为结论。

在油富/水富两相区，B冻结为各纯相参考。进入混合单相区后，冻结黏度规则须另行预登记为单值光滑反事实函数并做敏感性，或只比较共同两相域；不能依相标签突然切黏度。物理单相使用kr=1，不强造油水两相。无法表示的拓扑/已知局部失稳返回MODEL_DOMAIN_UNSUPPORTED，不继续扩展通用CPA。

## 6. 方程与输入闭合

守恒组分i=W,H：d[phi sum_alpha S_alpha rho_alpha w_i,alpha]/dt + div[sum_alpha rho_alpha w_i,alpha u_alpha]=q_i；u_alpha=−k kr_alpha/mu_alpha grad p，sum S=1。本案无反应、无重力、共同压力、物理弥散0。B/C有效域内要求相间逸度相等。EOS摩尔组成转质量组成：w_i=x_i MW_i/sum(x_j MW_j)。

水复用现有IAPWS验证；不能把纯水闭合或CO2专属闭合直接套到Heavy富水混合相。LBC外推只算模型假设，不算实测降黏。

**缺口不是CPA全部开发，而是本案使用的平衡/密度与黏度两个输入文件。** 允许目标数据表或明确假设版本用于条件性机制研究；真实样品定量结论必须有同域证据。Heavy纯组分表可冻结做初始筛查，不能证明混合物正确。

pr_binary_matrix_screening.csv中Heavy的0.2398345519及温度系数是BENCHMARK_ONLY/NO_FIT的squalane诊断值，不自动赋给本案真实Heavy。默认真实Heavy kij为空；禁止根据增产挑参数。

A可先做声明纯物性假设的数值基线；B/C在指定提供者后开展条件性计算，不等同实验复现。真实L/A/k/PV/kr/Pc缺失只阻止实物映射，不阻止归一化协议设计。

最小PVT检查：653.15K，28/29/30MPa，zW=0/.01/.05/.1/.2/.35/.5/.65/.8/.9/.97/.995/1，共39点；实际轨迹中追加检查。输出相数/组成/相量/rho/mu/逸度/守恒/稳定性和外推标签；正反混合路径核对，39点不代表连续域认证。LTE物理对应需tau_equil<<tau_res<<tau_reaction；无时间证据时仅报告平衡交换极限。

## 7. 运行矩阵

M0：上述最小PVT与无交换守恒/单相、Buckley–Leverett型控制。
M1：A380/B380/C380，Nx200，0–2PVI；C物性未齐不编造结果。
M2：已启用组的Nx100/400、时间步减半，检查数值混合是否改变携带结论。
M3：porous_media/corey_sensitivity.csv的低/中/高滞留代理敏感性，不冒充一套SCW实测曲线。
M4：360°C/28MPa对照，先比同温A/B/C，再比较温差。温差不是临界点独立因果效应；固定初始饱和度和PV时不同温度的初始质量会改变，必须记录各自M_H0并归一化，不能同时强称M、S、PV、T、p全部相同。
M5：加入Diesel，H2O+Diesel+Heavy；烃质量比1:1作为明确设计控制，或恢复来源四馏分。此时才研究轻重选择性；随后二维60×20×1。

速率敏感性不在本轮被解释为传质动力学：LTE无毛管/弥散时相同PVI变速不自动验证停留时间机制。有限传质、弥散、毛管及物理时钟后置。

## 8. 观测量

PVI=integral Qinj,actual,res dt/PV0；RF_H=M_H,prod/M_H0。PV0固定为各次运行初始操作态可接触孔隙体积，不能用常温泵体积直接替代储层体积。

组分产出取真实边界/井源；按相分解m_dot_H,alpha=integral rho_alpha w_H,alpha u_alpha·n dA，仅计外流。分别积分水富相携带与油相直接产出并闭合总量；混合单相单列，不强归为水携带。未建模乳化/颗粒不归入溶解。

记录F_H,wr=M_H,wr/M_H,prod，分母0时NaN；Delta RF_BA、Delta RF_CB；p、S、x、rho、mu、kr/mu剖面；实际流率、BHP切换、压差。液力功integral Delta p Qactual dt不包含加热能耗。

阶段二用去水质量组成Y_iHC和质量基准E_LH=(m_dot_L/m_dot_H)/(M_L0/M_H0)，不得混摩尔/质量基准。Heavy近检出限时给NaN/检出限，不报无限富集。

必交图：RF–PVI、分相携带累计量、出口Heavy分相流率、压差/实际流率、固定PVI剖面、网格/物性/kr敏感性。单Heavy没有轻质化结论。

## 9. 预定验收与程序任务

数值目标：组分全程累计闭合<=1e-6；B/C活动相log逸度差<=1e-6；归一化误差<=1e-8；无负相量，rho/mu正且有限。零库存组分不强制逸度相等。A水相Heavy通量为0；C退回B黏度须重现B；只积分接受时间步。

网格/时间步复核：RF(1/2PVI)绝对差<=.005，携带量和压差有信号时相对差<=2%；机制差异需超过估计数值误差3倍，否则未能分辨。物性/kr假设内效应反号时报告条件性，不选择最好看的情况。上述是设计数值目标，不是文献物性误差。

建议新case/scw_heavy_mechanism_1d（尚未实现）。复用既有网格、传输、井控和accepted-step台账；新增A/B/C控制适配、黏度反馈开关、producer_phase_component.csv和物性域标签。已有EOS选择不等于这些开关已实现；不得报告本案已经编译或运行成功。毛管尚未在现有流动核完整接入；真实岩心需实现它，或用本试件的力尺度证明Pc=0可接受。

实体实验最少补：同批Heavy表征；目标T/P相组成、密度和mu；保温/接触时间空白验证非反应/LTE窗口；实际试件L/A/PV/k及泵/背压限制；关键驱替重复、取样降压二次分相与质量损失校正。装置安全由实验室审核，数值30MPa不构成安全许可。

## 10. 停止规则

本轮交付设计和调研，不生成结果曲线。下一步限定H01输入闭合、A/B最小驱替和有依据的C，不开发新CPA分支。条件性机制结果可先做；真实样品定量预测仍需适用证据。FIGURE7_REPRODUCTION=BLOCKED、FORMAL_FLOW=BLOCKED继续保留，不把算法CI绿色当成实验复现。
