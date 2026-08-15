# 公开人类椎间盘退变队列中髓核转录模块方向性与异质性的队列感知描述性分析

## 标题页

**文章类型：** 原创研究

**短题名：** IVDD 髓核转录模块的队列级异质性

**作者：** [待全体作者确认后补入]

**作者单位：** [待全体作者确认后补入]

**通讯作者：** [待补入姓名、通信地址、电子邮箱和电话]

**英文题名：** Cohort-Aware Descriptive Analysis of Directionality and Heterogeneity of Nucleus Pulposus Transcriptional Programs across Public Human Intervertebral Disc Degeneration Datasets

**图表：** 图形摘要；2 张主图；2 张主表；5 张补充图；9 组补充表

## 摘要

公开椎间盘退变（intervertebral disc degeneration，IVDD）单细胞研究常含有大量细胞，却只有较少的独立供体或样本观察键；不同矩阵处理流程和临床取材方式也使跨队列比较复杂。我们使用 NCBI Gene Expression Omnibus（GEO）的公共功能基因组学记录[1]，并按照单细胞研究中忽视生物学重复可能产生假阳性的警示[2]，将细胞视为嵌套观测而非独立重复。本研究围绕四个项目内预先锁定的髓核（nucleus pulposus，NP）转录模块开展证据受限、队列特异性的分析：细胞外基质和胶原重塑、炎症/NF-κB 反应、缺氧/氧化应激及椎间盘基质稳态。统计观察单位为供体，或在公共标识不完整时元数据中明确标记但患者层面未完全核实的推定供体/样本/文库键。细胞仅作为嵌套观测。我们将 GSE230809 父项目（GSE229711 与 GSE230808）作为单一探索性 AF/NP 来源项目[3-6]，并分别展示 GSE244889[7,8]、GSE153066[9] 和 GSE165722[10,11] 的 NP 分数。原始 10x 矩阵仅在完成归档、质量控制、注释和样本身份审计后汇总；GSE153066 使用贡献者保留的稠密计数矩阵，GSE165722 仅用于标准化计数的分数层面方向分析。冻结的默认汇总报告逐队列未加权组间差、Welch 95% 区间、供体/文库键自助法区间和逐一剔除观察键稳定性，不报告 p 值、正式荟萃分析或复制裁决。另行打包的 S7、S8 和 S9 为采用 REML 与 Knapp-Hartung 区间的异方差标准化均数差（SMDH）随机效应汇总，其 HKSJ 和 BH p 值仅作透明报告。

在默认 NP 汇总的四个队列中，缺氧/氧化应激模块在记录严重度较高组中的点估计均为正，但相应的四个 Welch 95% 区间均包含 0。细胞外基质、炎症/NF-κB 和椎间盘基质稳态模块的方向在队列间不一致。S8 事后六队列扩展中，缺氧模块的 SMDH 为 0.7694（95% CI 0.1706 至 1.3682）；然而，S9 以 GSE245147 替换 GSE167931 后，该值为 0.5746（95% CI -0.7231 至 1.8723）。GSE251686 的单独审计分析因一条记录未通过数据流完整性检查且可用比较不平衡而保持隔离。本研究提供了可审计、以队列为单位的描述，指出模块分数在何处呈现描述性方向对齐、在何处呈现异质性；结果并未建立普遍 IVDD 转录程序、机制、生物标志物或治疗靶点。

## 关键词

椎间盘退变；髓核；公共转录组数据；单细胞转录组；队列异质性；可重复性

## 图形摘要

公共人类 IVDD 队列以供体或推定样本/文库键作为观察单位进行审计，细胞保留为嵌套观测。四个模块在外部评分前锁定。四个默认 NP 对比中，缺氧/氧化应激的点估计均为正，但相应 Welch 区间均包含 0，其他模块方向不一致。GSE251686 的隔离敏感性分析未进入默认汇总。

[[FIGURE:graphical_abstract]]

## 引言

椎间盘退变是具有重要临床负担且生物学表现高度异质的疾病过程。公共单细胞转录组数据为比较人类样本中的分子程序提供了机会[3,7,10,12,13]，但细胞数不能替代独立的供体或样本键观察数[2,14]。直接跨队列合并细胞可能使方向估计看似精确，却掩盖供体嵌套、临床来源差异和数据处理边界。

近期人类 IVDD 研究已整合多个公共数据集，或聚焦于特定 NP 状态[12,13,15]。因此，本研究不声称是首次进行 IVDD 多队列整合。其更窄的目标是保存可获得的观察单位结构，并报告各队列特异性的 NP 转录模块方向在哪些地方一致、在哪些地方不同。

本研究的问题是，四个项目内预先锁定的 NP 表达模块在独立处理的公共数据集中是否表现出方向对齐或方向不一致。研究被设计为一项关注来源可追溯性的描述性证据审计，而不是用于估计年龄独立的疾病效应、建立细胞谱系或分子机制、开发临床分类器，或提名治疗靶点。

## 材料与方法

### 研究设计、队列和推断边界

表 1 和补充表 S2 描述队列角色、记录的分组结构、观察键和身份核对。本文中，
AF 指纤维环（annulus fibrosus），ECM 指细胞外基质（extracellular matrix），
NF-κB 指核因子 κB（nuclear factor-kappa B）。GSE229711 与 GSE230808 被视为单一 GSE230809 父项目，而非发现队列与验证队列。该项目每个部位仅有 3 名年轻低级别供体，并与年龄较大的晚期退变供体比较，因此年龄与记录疾病状态完全混杂。该父项目的 AF 和 NP 对比仅作为与晚期退变相关的探索性展示。

外部 NP 来源队列分别独立处理。GSE244889 提供 4 个轻度和 3 个重度、由标题推导的推定供体/文库键。GSE153066 提供 8 个相对正常和 8 个退变的样本前缀键，但其临床来源和年龄与疾病状态混杂。GSE165722 按来源论文的严重度分组提供 4 个轻度和 4 个重度推定样本键；GEO 将其整数样矩阵描述为标准化计数，因此该队列仅用于分数层面的方向和稳定性分析[10,11]。GSE251686 接受了独立审计的 2 对 3 探索性分数分析[16,17]。由于 `GSM7986002` 未通过数据流完整性审计，该记录被永久排除。该不平衡展示未纳入冻结的默认 20 效应汇总，因为其存在完整性失败以及样本身份和协变量信息的局限；它既不被用作验证，也不被称为复制。补充图 S2 展示默认分析流与隔离分析流的队列处置。

队列来源依据相应 GEO 系列和可获得的来源论文进行了核对，包括 GSE244889[7,8]、GSE153066[9] 和 GSE251686[16,17]。GSE153066 的 GEO 记录当前没有关联文献，故以数据库记录而非未经核实的文章进行引用[9]。

在默认分析冻结后，GSE186542 被审计为小样本 NP 分数层面对比，即早期 Pfirrmann I--III 与晚期 IV--V 之间的 3 对 3 比较[22]。其 GEO SOFT 元数据未列出 PubMed 关联文献，因此仅以登录号引用。GSE167931 提供一份正常与退变之间的 每百万比对读段每千碱基转录本片段数（FPKM）表示（4 对 5），其配对 每百万转录本（TPM）矩阵仅保留为同一样本的处理敏感性检查[23-25]。GSE245147 在排除传代和处理组后，提供原生 Degenerated 与 No-degenerated 的 每百万比对读段每千碱基读段数（RPKM）子集（3 对 3）[26,27]。由于公共元数据不能确定二者的数据来源家族关系和患者级独立性，GSE167931 与 GSE245147 从不合并。

### 锁定模块评分和效应展示

四个模块定义在外部评分前已在项目内锁定。这是带时间戳和版本信息的分析锁定，不是前瞻性研究注册。补充表 S3 记录锁定时间、来源标识、基因列表和 SHA-256 哈希。每个分数要求至少 80% 的锁定基因能够映射。模块组成参考 KEGG 与 Reactome 知识库[18,19]以及引用的 IVDD 来源研究[3,10]。较高分数仅表示列表基因表达较高，并不表示获益、有害、因果关系或治疗相关性。对于已核实的原始 10x 矩阵，按供体或文库在解剖来源受限、通过质量控制的细胞内进行汇总，每个模块分数为映射基因的平均 每百万计数（CPM）的 `log1p(CPM)` 值。对于外部提供的稠密或标准化矩阵，同一评分定义仅在其声明的矩阵处理边界内应用。原始归档、输入哈希和评分到台账身份键的对应关系保留在可重复性契约中（补充表 S6）。

对于每个队列、部位和模块，计算记录严重度较高组与较低组之间的未加权均值差。Welch 区间与组内供体/文库键自助法区间作为描述性不确定性展示。逐一剔除观察键（leave-one-key-out，LOKO）分析用于评估符号保留情况，完整默认分析结果见补充表 S4。冻结的默认跨队列方向对齐展示不包含合并效应、正式荟萃分析、p 值或复制裁决。

### 独立打包的探索性标准化汇总

S7--S9 为每个 NP 队列和模块使用一个“记录严重度较高组减较低组”的效应。主要效应为异方差标准化均数差（SMDH），采用限制性最大似然法（REML）随机效应模型和 Knapp-Hartung 区间拟合。以传统合并标准差计算的 Hedges *g* 和 Paule-Mandel tau 平方估计用于模型敏感性分析。三个包均记录 `metafor` REML 的 `maxiter = 10000`，其余 `metafor` 控制参数保持默认。四个模块的 Hartung--Knapp--Sidik--Jonkman（HKSJ）p 值及其 Benjamini--Hochberg（BH）调整仅用于透明报告，不作为确认或复制检验。

S7 对四个默认 NP 队列进行独立标准化（`k = 4`）。S8 为事后六队列分数层面扩展，加入 GSE186542 和 GSE167931 FPKM（`k = 6`）。S9 使用相同设计，但以 GSE245147 的原生子集替换 GSE167931（`k = 6`）。这些包不修改冻结的默认 20 效应汇总。较少的队列数、跨平台评分尺度、未核实的患者级独立性和来源家族不确定性均排除了确认性解释。

以供体为单位处理细胞符合多样本单细胞群体层面分析原则[14]。由于两组方差无需相等，未加权均值差采用 Welch 区间展示[20]；百分位数区间通过在每组内重抽样供体/文库键生成[21]。

### 计算可重复性

默认汇总使用 95% Welch 区间，并对每个对比进行 10,000 次独立的供体或文库键自助法重抽样。队列、部位和模块特异的确定性随机种子均由根种子 20260814 推导。映射基因比例低于 80% 的评分会连同审计记录一起排除。S7--S9 汇总使用 R 4.4.1、`metafor` 4.8.0 和 REML 控制参数 `maxiter = 10000`。输入、输出、环境和生成器哈希已在补充表 S6 与提交支持清单中索引。

## 结果

冻结的默认描述性汇总包含来自表 1 所列四个评分队列的 20 个队列/部位/模块效应，评分到台账的样本键 55/55 精确匹配（补充表 S2 和 S6）。所有四个队列均标记为 `confirmatory_eligible=false`。GSE230809 提供探索性 AF 和 NP 对比；GSE244889 提供 4 对 3 的 NP 方向支持对比；GSE153066 提供 8 对 8 的 NP 计数层面支持对比；GSE165722 提供 4 对 4 的 NP 标准化计数分数层面方向对比。GSE251686 不属于默认汇总、方向对齐或任一主图（补充图 S2）。

在 NP 中，缺氧/氧化应激模块在四个默认对比的记录严重度较高组中均出现正点估计：GSE230809 为 0.1776，GSE244889 为 0.2381，GSE153066 为 0.0882，GSE165722 为 0.4346（表 2 和图 1--2）。相应的描述性 Welch 95% 区间分别为 `[-0.0245, 0.3798]`、`[-0.1314, 0.6076]`、`[-0.2166, 0.3929]` 和 `[-0.0619, 0.9310]`，每个区间均包含 0。因此，四个正点估计构成描述性符号模式，并不是稳健的跨队列复制结果。细胞外基质/胶原重塑在 GSE230809、GSE244889 和 GSE165722 为正，而在 GSE153066 为负。炎症/NF-κB 在 GSE244889 和 GSE153066 为正，而在 GSE230809 和 GSE165722 为负。椎间盘基质稳态在 GSE230809、GSE244889 和 GSE165722 为正，而在 GSE153066 为负（表 2 和图 1--2）。

单独审计的 GSE251686 敏感性分析在永久排除 `GSM7986002` 后保留 5 个推定样本/文库键（轻度 n=2，重度 n=3）。其缺氧/氧化应激点估计为负（`-0.0778`，Welch 95% 区间 `[-1.7214, 1.5658]`），四个模块效应见补充表 S1 和补充图 S1。该隔离结果未被合并、未计入默认方向对齐，也未用于改变四队列默认结果。

独立的 S7 汇总中，细胞外基质/胶原重塑、炎症/NF-κB、缺氧/氧化应激和椎间盘基质稳态的 SMDH 分别为 1.0481（95% CI -0.9965 至 3.0926）、0.2809（-1.1432 至 1.7050）、0.8184（-0.1285 至 1.7654）和 0.1453（-0.3679 至 0.6585）（补充表 S7a--S7d 和补充图 S3）。这四个探索性区间均包含 0。

事后 S8 扩展中，细胞外基质、炎症、缺氧和基质稳态的 SMDH 分别为 0.7780（-0.3532 至 1.9093）、0.4032（-0.3991 至 1.2055）、0.7694（0.1706 至 1.3682）和 0.3762（-0.3017 至 1.0540）（补充表 S8a--S8d 和补充图 S4）。缺氧模块的 HKSJ p 值为 0.0214，四模块 BH 值为 0.0856，这些值仅作透明报告。在以 GSE245147 替换 GSE167931 的 S9 中，相应 SMDH 分别为 1.0931（-0.9186 至 3.1048）、0.4056（-1.0004 至 1.8115）、0.5746（-0.7231 至 1.8723）和 0.1046（-1.0200 至 1.2292）（补充表 S9a--S9d 和补充图 S5）。S9 的四个 BH 值分别为 0.6133、0.6556、0.6133 和 0.8205。替换后 S8 缺氧区间不再离开 0，是数据来源家族敏感性的结果，而不是确认。

AF 结果仅来自 GSE230809 父项目。由于低级别 AF 组仅有 3 名供体，且其年龄与晚期退变组完全混杂，这些估计保留为探索性背景，而不是外部证据。完整的默认逐一剔除观察键结果见补充表 S4。

锁定的来源项目保留细胞资格阈值敏感性分析没有改变 GSE230809 的输入集。全部 24 个供体/文库在 20、30 和 50 个细胞阈值下均通过，最小解剖来源受限且质量控制通过的细胞数为 471。96 个文库-模块分数和全部 8 个 AF/NP 描述性效应在这些阈值下完全相同（补充表 S5a 和 S5b）。该结果仅说明没有文库接近所选资格阈值，并不检验替代注释、细胞组成或随机细胞下采样下的稳健性。

## 讨论

这项以队列为单位的分析将一个具有审计支持但强度有限的观察与当前数据无法维持的强主张区分开来。四个默认 NP 对比中，缺氧/氧化应激在记录严重度较高组的点估计均为正，但每个相关 Welch 区间均包含 0。因此，这些结果不证明复制、机制、严重度生物标志物或治疗靶点。细胞外基质、炎症和基质稳态分数的方向不一致说明，不能凭借合并细胞或单一队列推断普遍的退变程序。

S7--S9 汇总提供了标准化效应量的描述，但没有消除这一限制。特别是，S8 中的缺氧区间未包含 0，而 S9 数据来源家族替换后的区间又包含 0。这种变化，加之事后队列选择、`k = 4` 或 `k = 6`、处理尺度差异和未知患者重叠，意味着其 p 值和区间都不能建立复制、生物标志物、机制或治疗效应。

本研究受到小样本量、若干资源中仅有推定而非完全核实的样本身份、GSE230809 中年龄与疾病状态完全混杂、GSE153066 中临床来源和年龄混杂、GSE165722 的标准化计数限制以及缺少独立 AF 严重度队列的限制。GSE251686 在评分后被透明地报告为隔离敏感性分析，而不是被选择性纳入默认汇总。四个分数是锁定的表达摘要而不是临床指数，且其绝对值不应跨平台合并。未来工作需要独立采集、样本量充分、样本嵌套关系已核实的人类队列，并结合正交生物学测量，才能提出因果、机制、预后或治疗结论。

## 结论

在纳入默认 NP 汇总的四个队列中，只有缺氧/氧化应激的点估计方向一致，且其四个描述性 Welch 95% 区间均包含 0。其余模块方向不一致。独立打包的 S7--S9 汇总同样为非确认性，且对数据来源家族替换敏感。这些结果支持透明的队列特异性报告，但不建立普遍的 IVDD 转录程序。

## 数据可用性

本研究分析的全部源数据可通过 NCBI GEO 获得，登录号包括 GSE230809、GSE229711、GSE230808、GSE244889、GSE153066、GSE165722、GSE251686、GSE186542、GSE167931 和 GSE245147。GSE56081 仅作为候选微阵列扩展接受审计，未用于模块效应分析。项目代码、配置、环境锁定、派生表和清单尚未归档至公共版本化仓库。投稿前，须将用于本文的确切版本存入公共仓库并使用可提供持久 DOI 的归档服务。请在此处补入最终仓库 URL 和 DOI：`[投稿前补入仓库 URL 和 DOI]`。

## 代码可用性

版本化的公开代码发布必须包含分析脚本、锁定配置文件、环境锁定、生成的结果表和补充表 S6 所引用的清单。本文不声称代码已经公开归档。

## 伦理声明

本研究重新分析已去标识化的公共数据，不涉及新的受试者招募、干预或标本采集。投稿前，通讯作者须确认适用的本地机构决定，并补入任何需要的审批、豁免或免除编号：`[投稿前补充伦理决定]`。

## 知情同意

本研究为已去标识化公共数据的二次分析，未新获取受试者知情同意。投稿前，通讯作者须确认是否需要披露任何数据集特异的知情同意信息。

## 经费

`[投稿前由作者补充经费来源、项目编号及资助方角色。]`

## 利益冲突

`[投稿前由全体作者确认并补充利益冲突声明。]`

## 作者贡献

`[投稿前由作者确认并补充作者姓名及 CRediT 贡献。]`

## 致谢

[待全体作者确认后补充。不得列入未经确认的贡献者或支持信息。]

## 参考文献

1. Barrett T, Wilhite SE, Ledoux P, Evangelista C, Kim IF, Tomashevsky M, Marshall KA, Phillippy KH, Sherman PM, Holko M, Yefanov A, Lee H, Zhang N, Robertson CL, Serova N, Davis S, Soboleva A. NCBI GEO: archive for functional genomics data sets-update. Nucleic Acids Research. 2013;41(Database issue):D991-D995. doi:10.1093/nar/gks1193.

2. Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R, Qaiser T, Matson KJE, Barraud Q, Levine AJ, La Manno G, Skinnider MA, Courtine G. Confronting false discoveries in single-cell differential expression. Nature Communications. 2021;12(1):5692. doi:10.1038/s41467-021-25960-2.

3. Swahn H, Mertens J, Olmer M, Myers K, Mondala TS, Natarajan P, Head SR, Alvarez-Garcia O, Lotz MK. Shared and Compartment-Specific Processes in Nucleus Pulposus and Annulus Fibrosus During Intervertebral Disc Degeneration. Advanced Science. 2024;11(17):e2309032. doi:10.1002/advs.202309032.

4. NCBI Gene Expression Omnibus. GSE230809: Shared and compartment-specific processes in nucleus pulposus and annulus fibrosus during intervertebral disc degeneration [Internet]. GEO Series accession GSE230809. 2024. SuperSeries containing GSE229711 and GSE230808; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230809.

5. NCBI Gene Expression Omnibus. GSE229711: The cellular landscape of the healthy human intervertebral disc [Internet]. GEO Series accession GSE229711. 2024. SubSeries of GSE230809; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229711.

6. NCBI Gene Expression Omnibus. GSE230808: Shared and compartment-specific processes in nucleus pulposus and annulus fibrosus during intervertebral disc degeneration [Internet]. GEO Series accession GSE230808. 2024. SubSeries of GSE230809; accessed 2026-08-14. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230808.

7. Chen F, Lei L, Chen S, Zhao Z, Huang Y, Jiang G, Guo X, Li Z, Zheng Z, Wang J. Serglycin secreted by late-stage nucleus pulposus cells is a biomarker of intervertebral disc degeneration. Nature Communications. 2024;15(1):47. doi:10.1038/s41467-023-44313-9.

8. NCBI Gene Expression Omnibus. GSE244889: Gene expression profile at single cell level of nucleus pulposus cells from mild and severe degenerative intervertebral discs [Internet]. GEO Series accession GSE244889. 2023. Accessed 2026-08-14; GEO record links PMID 38167807 and later related publications. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244889.

9. NCBI Gene Expression Omnibus. GSE153066: Single cell sequencing of human nucleus pulposus [Internet]. GEO Series accession GSE153066. 2023. Accessed 2026-08-14; GEO record explicitly lists citation as missing. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153066.

10. Tu J, Li W, Li W, Yang S, Yang P, Yan Q, Wang S, Lai K, Bai X, Wu C, Ding W, Cooper-White J, Diwan A, Yang C, Yang H, Zou J. Single-Cell Transcriptome Profiling Reveals Multicellular Ecosystem of Nucleus Pulposus during Degeneration Progression. Advanced Science. 2022;9(3):e2103631. doi:10.1002/advs.202103631.

11. NCBI Gene Expression Omnibus. GSE165722: Single-cell transcriptome profiling reveals nucleus pulposus heterogeneity and immunity during degeneration progression [Internet]. GEO Series accession GSE165722. 2021. Accessed 2026-08-14; GEO supplementary matrices are described as normalized counts. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165722.

12. Wang D, Li Z, Huang W, Cao S, Xie L, Chen Y, Li H, Wang L, Chen X, Yang JR. Single-cell transcriptomics reveals heterogeneity and intercellular crosstalk in human intervertebral disc degeneration. iScience. 2023;26(5):106692. doi:10.1016/j.isci.2023.106692.

13. Sun Y, Peng Y, Su Z, So KKH, Lu Q, Lyu M, Zuo J, Huang Y, Guan Z, Cheung KMC, Zheng Z, Zhang X, Leung VYL. Fibrocyte enrichment and myofibroblastic adaptation causes nucleus pulposus fibrosis and associates with disc degeneration severity. Bone Research. 2025;13(1):10. doi:10.1038/s41413-024-00372-2.

14. Crowell HL, Soneson C, Germain PL, Calini D, Collin L, Raposo C, Malhotra D, Robinson MD. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nature Communications. 2020;11(1):6077. doi:10.1038/s41467-020-19894-4.

15. Niu H, Qi H, Zhang P, Meng H, Liu N, Zhang D. Single-Cell Analysis Reveals Aspirin Restores Intervertebral Disc Integrity via Ferroptosis Regulation. Journal of Inflammation Research. 2025;18:6889-6905. doi:10.2147/JIR.S519218.

16. Jia S, Liu H, Yang T, Gao S, Li D, Zhang Z, Zhang Z, Gao X, Liang Y, Liang X, Wang Y, Meng C. Single-cell sequencing reveals cellular heterogeneity of nucleus pulposus in intervertebral disc degeneration. Scientific Reports. 2024;14(1):27245. doi:10.1038/s41598-024-78675-x.

17. NCBI Gene Expression Omnibus. GSE251686: Single-cell sequencing reveals cellular heterogeneity of nucleus pulposus in intervertebral disc degeneration [Internet]. GEO Series accession GSE251686. 2024. Accessed 2026-08-14; GSM7986002 is excluded in the present audit for failed stream integrity. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251686.

18. Kanehisa M, Furumichi M, Sato Y, Kawashima M, Ishiguro-Watanabe M. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research. 2023;51(D1):D587-D592. doi:10.1093/nar/gkac963.

19. Jassal B, Matthews L, Viteri G, Gong C, Lorente P, Fabregat A, Sidiropoulos K, Cook J, Gillespie M, Haw R, Loney F, May B, Milacic M, Rothfels K, Sevilla C, Shamovsky V, Shorser S, Varusai T, Weiser J, Wu G, Stein L, Hermjakob H, D'Eustachio P. The Reactome pathway knowledgebase. Nucleic Acids Research. 2020;48(D1):D498-D503. doi:10.1093/nar/gkz1031.

20. Welch BL. The generalization of Student's problem when several different population variances are involved. Biometrika. 1947;34(1-2):28-35. doi:10.1093/biomet/34.1-2.28.

21. Efron B. Bootstrap methods: another look at the jackknife. The Annals of Statistics. 1979;7(1):1-26. doi:10.1214/aos/1176344552.

22. NCBI Gene Expression Omnibus. GSE186542: Nucleus pulposus related lncRNA and mRNA expression profiles in intervertebral disc degeneration [Internet]. GEO Series accession GSE186542. 2021. Submitted 2021; accessed 2026-08-14; GEO SOFT lists no linked PubMed citation. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186542.

23. NCBI Gene Expression Omnibus. GSE167931: Next Generation Sequencing analysis at single-cell level of normal and degenerated nucleus pulposus cells transcriptomes [Internet]. GEO Series accession GSE167931. 2021. Accessed 2026-08-14; GEO SOFT lists PMIDs 35304463 and 35340126. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167931.

24. Li G, Ma L, He S, Luo R, Wang B, Zhang W, Song Y, Liao Z, Ke W, Xiang Q, Feng X, Wu X, Zhang Y, Wang K, Yang C. WTAP-mediated m6A modification of lncRNA NORAD promotes intervertebral disc degeneration. Nature Communications. 2022;13(1):1469. doi:10.1038/s41467-022-28990-6.

25. Li G, Luo R, Zhang W, He S, Wang B, Liang H, Song Y, Ke W, Shi Y, Feng X, Zhao K, Wu X, Zhang Y, Wang K, Yang C. m6A hypomethylation of DNMT3B regulated by ALKBH5 promotes intervertebral disc degeneration via E4F1 deficiency. Clinical and Translational Medicine. 2022;12(3):e765. doi:10.1002/ctm2.765.

26. NCBI Gene Expression Omnibus. GSE245147: CytoDNA triggered NP cell inflammatory senescence via cGAS-STING axis sensing but not AIM2 inflammasome activation [Internet]. GEO Series accession GSE245147. 2024. Accessed 2026-08-14; GEO SOFT lists PMID 38488012. Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245147.

27. Zhang W, Li G, Zhou X, Liang H, Tong B, Wu D, Yang K, Song Y, Wang B, Liao Z, Ma L, Ke W, Zhang X, Lei J, Lei C, Feng X, Wang K, Zhao K, Yang C. Disassembly of the TRIM56-ATR complex promotes cytoDNA/cGAS/STING axis-dependent intervertebral disc inflammatory degeneration. The Journal of Clinical Investigation. 2024;134(6):e165140. doi:10.1172/JCI165140.

## 补充表与图题注

**补充表 S1. GSE251686 隔离探索性敏感性分析的效应。** 本表报告永久排除 `GSM7986002` 后轻度 n=2 与重度 n=3 比较的四个模块效应。该队列不属于默认汇总、方向对齐或主图，且不作为验证或复制证据。

**补充表 S2. 队列处置和身份审计。** 本表列出默认和隔离分析流、观察键、分组结构、身份或资格核对、默认汇总纳入状态和解释边界。

**补充表 S3. 锁定的转录模块定义。** 本表给出来源类别、标识符、完整锁定基因列表、时间戳、评分方向约定和每个模块的 SHA-256 哈希。

**补充表 S4. 默认汇总的完整逐一剔除观察键结果。** 本表记录默认 20 效应汇总中的每次逐一剔除观察键效应计算，包括来自 GSE230809 父项目的 AF 和 NP 展示。

**补充表 S5a. 来源项目保留细胞阈值汇总。** 本表记录 20、30 和 50 个细胞资格运行、通过的文库数以及最小解剖来源受限且质量控制通过的细胞数。

**补充表 S5b. 来源项目保留细胞阈值下的效应稳定性。** 本表将三个阈值下的全部 8 个 GSE230809 AF/NP 描述性效应与 30 细胞参考运行进行比较。

**补充表 S6. 可重复性契约。** 本表以 SHA-256 哈希索引默认汇总清单、身份交叉表、锁定程序台账、阈值敏感性输入、隔离 GSE251686 清单、S7--S9 清单和生成器、产物生成器以及项目本地 Python 环境锁定。

**补充表 S7a--S7d. 四队列探索性 NP 随机效应效应量汇总。** SMDH 指异方差标准化均数差，REML 指限制性最大似然法。这些表给出 S7 的研究层效应、主要 SMDH/REML/Knapp-Hartung 结果、Hedges *g* 和 Paule-Mandel 敏感性以及逐一剔除队列结果。该分析为单独打包的非确认性探索性汇总，不替代默认描述性分析。

**补充表 S8a--S8d. 事后六队列 NP 扩展。** SMDH 指异方差标准化均数差，REML 指限制性最大似然法。这些表在加入 GSE186542 和 GSE167931 FPKM 后给出单独打包的 S8 研究层效应、主要结果、模型敏感性和逐一剔除队列结果。HKSJ 和 BH p 值仅作透明报告，不构成患者级验证或确认。

**补充表 S9a--S9d. NP 数据来源家族替换敏感性分析。** SMDH 指异方差标准化均数差，REML 指限制性最大似然法。这些表给出以 GSE245147 的原生对比替代 GSE167931 后的 S9 结果。两个数据来源家族队列从不合并；S9 是非确认性敏感性分析，而不是第七个独立队列。

**补充图 S1. GSE251686 隔离探索性敏感性展示。** 本图展示 GSE251686 的队列特异性效应点估计、描述性 Welch 95% 区间和逐一剔除观察键方向保留。该展示保持隔离，不进入默认汇总。

**补充图 S2. 队列处置和分析边界。** 流程图区分四个默认评分队列与单独评分的 GSE251686 敏感性包，并说明分析的推断边界。

**补充图 S3. 四队列探索性 NP 随机效应效应量汇总。** 森林图给出 S7 SMDH 估计和 Knapp-Hartung 区间；SMDH 指异方差标准化均数差。该图为单独的探索性定量摘要，不是确认或复制展示。

**补充图 S4. 事后六队列 NP 扩展。** 森林图给出加入 GSE186542 和 GSE167931 FPKM 后的 S8 结果。图示 SMDH 和区间均为非确认性；HKSJ 与 BH p 值仅在补充表 S8b 中报告。

**补充图 S5. NP 数据来源家族替换敏感性分析。** 森林图给出以 GSE245147 原生子集替代 GSE167931 后的 S9 结果。该图显示数据来源家族敏感性，不增加独立验证队列。不同面板的横轴范围可能不同，跨处理尺度的效应量不应直接比较。

## 主表

**表 1. 默认描述性汇总中各队列的角色与推断边界。**

[[TABLE:1]]

**表 2. 默认描述性汇总中各队列特异性的 NP 模块效应。** 默认层共包含 20 个效应（16 个 NP 效应和 4 个探索性 AF 效应）；本表展示其中 16 个 NP 效应。

[[TABLE:2]]

## 主图

**图 1. 各队列的 NP 模块分数差及描述性 Welch 95% 区间。** 颜色仅用于区分队列，不表示效应大小。

[[FIGURE:figure_1]]

**图 2. NP 队列特异性方向和描述性方向对齐。** 蓝色和橙色仅表示正向和负向，不表示效应大小。

[[FIGURE:figure_2]]

## 补充材料说明

补充表 S1--S4、S5a--S5b、S6、S7a--S7d、S8a--S8d、S9a--S9d 以及补充图 S1--S5 应作为独立的投稿附件提供。S7--S9 的每个文件均应保留其单独打包、探索性和非确认性的标签。图形摘要可按目标期刊要求单独上传。
