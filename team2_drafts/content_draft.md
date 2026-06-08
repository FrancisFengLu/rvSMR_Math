# Content Draft — Sub-worker 2A

Source-of-truth: `main.tex` §0–§14, `briefing_for_wei.md` §§1–4, `steps_5_to_9_logic.md`, `citation_audit_2026-05-27.md`.

This is the prose draft per step. Each step follows: 直觉 → 数学 → 用了哪篇 paper → 替代方案 → 代码位置 → 陷阱.

---

## Step 0 · Assemble inputs (fix notation)

### 直觉

在我们开始任何因果推断之前, 必须把每一项数字的"身份证"固定下来. rvSMR 的一行 input 不是一个 SNP 的 effect, 而是 "(基因 g, 变体掩码 k, 细胞类型 c) 这个三元组的 burden effect" — 也就是把 mask k 内所有 rare variant 用同一个加权和压成单一信号后, 它对 X (表达) 和对 Y (结局) 的回归斜率. 同一个三元组在 exposure 侧 (单细胞 eQTL) 和 outcome 侧 (人群-scale exome GWAS) 必须用 *完全相同* 的 mask 定义和 *完全相同* 的权重, 否则 $\hat b_x$ 和 $\hat b_y$ 指的就不是一个东西, 后面再做的所有 ratio 都没有意义.

我们想知道什么: "对于基因 g, 在细胞类型 c 中, 它的 mask-k burden 表达水平如果被人为干预上调 1 个单位, 结局 Y 会变多少?" 这是个 *局部因果斜率* (per-allele causal effect on Y per unit of X).

所以 Step 0 要做的就是 *固定记号*, 让所有后续公式都说同一种语言.

### 数学

对每个 $(g, k, c)$:

- 暴露侧: $(\hat b_{x,k,c},\ \mathrm{SE}_{x,k,c},\ n_x)$ — SAIGE-QTL Step 2 group-burden score test 在 TenK10K Phase 1 单细胞 RNA-seq 上跑出来.
- 结局侧: $(\hat b_{y,k},\ \mathrm{SE}_{y,k},\ n_y)$ — Genebass 或 RGC-ME 的 gene-mask burden test.
- 可选 pQTL anchor: $(\hat b_{\mathrm{burden}\to\mathrm{protein},k},\ \mathrm{SE}_{\mathrm{burden}\to\mathrm{protein},k})$ — UKB-PPP rare-variant 或 deCODE.
- Sample overlap: $\rho$ (K=1 scalar) 或 $R_{xy}\in\mathbb{R}^{K\times K}$ block.

在固定 $(g, c)$ 时把 K 个 mask 堆成向量:
$$
\hat{\boldsymbol b}_x = (\hat b_{x,1},\dots,\hat b_{x,K})^\top,\quad \hat{\boldsymbol b}_y = (\hat b_{y,1},\dots,\hat b_{y,K})^\top,\quad D_x = \mathrm{diag}(\mathrm{SE}_x),\ D_y = \mathrm{diag}(\mathrm{SE}_y).
$$

### 这一步用了哪篇 paper 的什么方法

- **两样本 IV 设定**: Pierce & Burgess 2013 *AJE* 178(7):1177 — efficient design for two-sample MR.
- **Summary-stat MR 框架**: Burgess-Butterworth-Thompson 2013 *Genetic Epidemiology* 37(7):658. (注意期刊是 *Genet Epidemiol*, 不是 Stat Med.)
- **样本重叠为零的默认假设**: Burgess et al. 2016 *Genet Epidemiol* 40(7):597 — bias due to participant overlap in two-sample MR.

### 为什么不用替代方案

不用 *individual-level* IV (e.g. 把 SAIGE-QTL 和 Genebass 的原始基因型重新合并跑一次 2SLS): rare variant data sharing 障碍极大 (HIPAA/GDPR), summary-stat 路径是唯一可扩展的方案. 不用 *single-variant* MR 框架直接对 rare variant 跑: MAF $\lesssim 10^{-3}$ 时第一阶段 F → 0, Wald 比退化 (Zhu et al. 2016 SMR 在 rare-variant 上失败的核心原因).

### 代码位置

`validate_summary_input()` in `rvMR/R/utils.R:32` — checks dimensions, finiteness, no-NA.

### 陷阱

- 把 RareEffect Step 4 BLUP / PEV (Bayesian shrinkage) 替换 SAIGE-QTL Step 2 score-test $\hat b$: BLUP 不是无偏估计, 下游 $\chi^2$ 参考分布失效.
- 把 RGC-ME 引为 "Sun BB 2024" — 错的, 真正第一作者是 **Sun KY** (Kathie Y. Sun). Sun BB 是 UKB-PPP common-variant proteomics 的不同人.
- $n_x$ 用全 donor 数而非 SAIGE-QTL pseudobulk 后的 cell-type-stratified effective N.
- 一样本分析时漏掉 $\rho \neq 0$, Step 4 / Step 6 的 denominator 会少 cross-term.

---

## Step 1 · Construct the linear burden instrument

### 直觉

Rare variant (MAF < 1‰) 单个跑没意义 — 第一阶段 F 太小, eQTL 信号没法被稳定地检测出来. 解决方案: **把 mask 内所有 rare variant 用权重相加**, 让单一变体的微弱信号在 burden 这层级聚合成一个能跑的信号. 这就是 burden test 的根本逻辑: 用"许多很罕见的变体合起来出现的次数"作为 IV 而不是"某一个特定变体".

为什么是 *linear* 加权和? 因为只有 linear 形式才能在 IV 框架里写出 Wald 比 $\hat b_y / \hat b_x$. SKAT 是二次型, ACAT 是 Cauchy-combine, 它们都不提供 *signed* 单一斜率 — 你没法对一个 p-value 做 ratio.

权重怎么选? 默认是 STAAR 推广的 Beta(MAF; 1, 25) ≈ $(1-p)^{24}$, 给 *更罕见* 的变体更高权重 (因为越罕见的越可能是 deleterious, 信噪比越高).

### 数学

$$
Z_k \;=\; \sum_{j \in \mathcal{M}_k} w_j\, G_j, \qquad w_j \propto \mathrm{Beta}(p_j; 1, 25) = (1 - p_j)^{24}.
$$

Madsen-Browning 原始权重: $w_j = 1/\sqrt{p_j(1-p_j)}$.

### 这一步用了哪篇 paper 的什么方法

- **Burden 加权和 $B = \sum_j w_j G_j$ 的开山之作**: Madsen BE, Browning SR 2009 *PLOS Genetics* 5(2):e1000384. "A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic."
- **STAAR Beta(1,25) MAF 权重**: Li X et al. 2020 *Nature Genetics* 52:969 (注意: **Li X** 2020, 不是 Li Z 2022 — Li Z 2022 是 STAARpipeline 的 *Nature Methods* 实现论文, 不是权重方法的原始 paper).
- **Beta(1,25) MAF 权重的更早出处**: Wu MC et al. 2011 *AJHG* 89(1):82 — SKAT 论文, Beta(1,25) 是它定义的 default weight, STAAR 沿用并 generalized 给 burden + SKAT + ACAT 三家.

### 为什么不用替代方案

- **CAST** (Morgenthaler-Thilly 2007 *Mutation Research* 615:28): 把"携带至少 1 个 rare allele"二值化, 损失剂量信息.
- **CMC** (Li B-Leal S 2008 *AJHG* 83:311): combined multivariate and collapsing — 类似 CAST 的指示函数, 同样损失剂量.
- **SKAT** (Wu 2011 *AJHG* 89:82): 二次型 $\sum w_j^2 (G_j^\top y)^2$. **没有 signed 斜率**, 不能跑 Wald 比.
- **ACAT** (Liu Y 2019 *AJHG* 104:410): Cauchy-combined p-values. 也无 signed 斜率.
- **STAAR-O**: STAAR 的 omnibus 组合 (Cauchy-combine burden + SKAT + ACAT). p-value 而非 effect.

rvSMR 选择 *linear weighted burden* 是因为它是 Wald-ratio 的*唯一*可识别 aggregator. STAAR-B (B = burden sub-statistic) 给我们的就是这个 linear form.

### 代码位置

`rvMR` 包消费 *预先计算好的* $\hat b_x, \hat b_y$ — burden 的实际拼装在上游 SAIGE-QTL / STAAR 完成. Downstream entry: `wald_burden()` in `rvMR/R/wald_burden.R:77`, roxygen 记 unit convention.

### 陷阱

- 直接把 STAAR-O p-value 当 burden 用 — STAAR-O 没有 signed effect.
- 用 SKAT / ACAT 当 IV — 无法做 Wald 比.
- *数据驱动* 重 scale 权重 $w_j$ (e.g. 根据样本 MAF 重新归一化): 破坏 IV 的 exogeneity.
- Mask 跨 gene 用 STAARpipeline 的默认 "putative regulatory" mask 集合而不重定义: STAARpipeline mask 可能与 cell-type 和 outcome 端的 mask 不匹配.

---

## Step 2 · Acquire per-mask burden summary statistics

### 直觉

Step 1 是 "怎么定义 burden"; Step 2 是 "去哪里把 burden 的回归斜率拿到". 两端必须 *同一组 variant + 同一组 weight* 跑出两个独立的 score test:

- 暴露端: 在 TenK10K 单细胞数据上跑 SAIGE-QTL Step 2 burden score test, 输出 $\hat b_x = \partial \log E[\mathrm{expression} | Z_k] / \partial Z_k$.
- 结局端: 在 Genebass 或 RGC-ME UKB exomes 上跑 SAIGE-GENE+ burden score test, 输出 $\hat b_y = \partial \log E[\mathrm{outcome} | Z_k] / \partial Z_k$.

两端必须用同一种 mask 和权重, 否则 ratio 无意义.

### 数学

$$
(\hat b_{x,k,c},\, \mathrm{SE}_{x,k,c}) \leftarrow \text{SAIGE-QTL Step 2 burden score test on } (g, k, c).
$$
$$
(\hat b_{y,k},\, \mathrm{SE}_{y,k}) \leftarrow \text{Genebass / RGC-ME burden test on } (g, k).
$$

### 这一步用了哪篇 paper 的什么方法

- **SAIGE-QTL Step 2 burden score test (单细胞 eQTL 版)**: Zhou W, Cuomo ASE et al. 2024 medRxiv 2024.05.15.24307317 — "SAIGE-QTL maps genetic regulation of gene expression in single cells with multi-ancestry sample sizes."
- **SAIGE-GENE+ 的 burden score test (人群 exome 版, 给 outcome 用)**: Zhou W et al. 2022 *Nature Genetics* 54(10):1466 — "SAIGE-GENE+ improves the efficiency and accuracy of set-based rare variant association tests."
- **Genebass UKB exomes substrate**: Karczewski KJ et al. 2022 *Cell Genomics* 2:100168 — 394,841 UKB exomes 上所有 phenotype 的预算 burden p-values + effect sizes.
- **RGC-ME alternative substrate**: Sun KY et al. 2024 *Nature* 631:583 — Regeneron Genetics Center, 983,578 individuals.

### 为什么不用替代方案

- **RareEffect Step 4 BLUP**: Zhou et al. 后续工作给出 *posterior mean* 的 random-effect estimator. BLUP 是 Bayesian shrinkage, 不是 无偏 estimator, 不再有 $\chi^2_1$ asymptotics — 下游 AR 的 reference distribution 失效. **关键的方法 lineage 区分**: rvSMR 只接 Step 2 score test, 不接 Step 4 BLUP.
- **STAARpipeline 输出**: Li Z et al. 2022 *Nat Methods* — STAARpipeline 是 pipeline 实现, 输出常是 STAAR-O p-value (Cauchy-combined), 不带 signed burden $\hat b$. 需要 STAAR-B 子统计量.
- **Burden-RVAT 类 collapsing 方法**: 输出 binary 携带者状态, 损失剂量.

### 代码位置

被 `wald_burden()` (`wald_burden.R:77`) 消费, `validate_summary_input()` (`utils.R:32`) 验证. **rvMR 不跑 SAIGE-QTL / Genebass — 它消费它们的输出.**

### 陷阱

- 用 RareEffect Step 4 BLUP 取代 Step 2 score-test $\hat b$: $\chi^2$ asymptotics 失效.
- 两端 mask 定义不一致: $\hat b_x$ 和 $\hat b_y$ 指代不同 IV, ratio 没意义.
- 二值结局 (Genebass binary) 直接给 log-OR 而 SAIGE-QTL 给连续表达 — 注意 scale 转换.
- $n_x$ 误用全 donor 数 (而非 cell-type-stratified pseudobulk).

---

## Step 3 · Per-mask Wald ratio and first-stage F

### 直觉

一旦两端的 $\hat b_x$ 和 $\hat b_y$ 在同一个 burden $Z_k$ 上, 经典 MR 答案就是 Wald 比 $\hat b_y / \hat b_x$. 直觉上: "burden 每多 1 单位拉动表达 $\hat b_x$, 同时拉动结局 $\hat b_y$, 那么 *表达每 1 单位拉动结局* 应该是 $\hat b_y / \hat b_x$". 这是 IV 思想的最朴素形式.

第一阶段 F = $(\hat b_x / \mathrm{SE}_x)^2$ 衡量 "burden 移动 X 有多稳". $F \gg 1$ 表示 burden 是个强 IV; $F \le 1$ 表示弱 IV.

**关键: rvSMR 不在 F = 10 处截断**. 因为 AR (Step 4) 在 F < 1 时仍保持 $1 - \alpha$ 覆盖率. 这是 rvSMR 区别于传统 MR 的根本设计选择.

### 数学

$$
\hat b_{xy} \;=\; \frac{\hat b_y}{\hat b_x}, \qquad F \;=\; \left(\frac{\hat b_x}{\mathrm{SE}_x}\right)^2.
$$

Delta-method SE (含 sample overlap):
$$
\widehat{\mathrm{Var}}(\hat b_{xy}) \;\approx\; \frac{\mathrm{SE}_y^2}{\hat b_x^2} \;+\; \frac{\hat b_y^2\, \mathrm{SE}_x^2}{\hat b_x^4} \;-\; 2\,\frac{\hat b_y}{\hat b_x^3}\,\rho\, \mathrm{SE}_x\, \mathrm{SE}_y.
$$

### 这一步用了哪篇 paper 的什么方法

- **Wald ratio MR 解释**: Burgess S, Labrecque JA 2018 *European Journal of Epidemiology* 33(10):947 — interpretation and presentation of causal estimates.
- **IV in epidemiology / Wald 比的因果识别**: Didelez V, Sheehan N 2007 *Statistical Methods in Medical Research* 16(4):309.
- **F > 10 obsolete 的现代替代**: Lee DS, McCrary J, Moreira MJ, Porter J 2022 *American Economic Review* 112(10):3260 — "Valid t-ratio inference for IV". 给出 tF critical values, **替代** Stock & Yogo 2005 的 F > 10 阈值.
- **关键设计选择**: rvSMR *不过滤 F > 10*. AR (Step 4) 在 F → 0 时仍 size-correct.

### 为什么不用替代方案

- **dIVW (debiased IVW)** (Ye T et al. 2021 *Ann Stat*): point-estimate debias, 但 CI 仍依赖 F > $\chi^2_{1, 1-\alpha} \approx 3.84$. F < 3.84 时退化.
- **MR-RAPS** (Zhao Q et al. 2020 *Ann Stat*): profile-score, robust to invalid IV; F 极弱时仍 unstable.
- **Stock-Yogo F > 10**: 自 LMMP 2022 起被 superseded. F > 10 在 just-identified case 不能保证 size correctness.
- **rvSMR 的 AR (Step 4)**: uniform coverage over IV strength, 包括 F < 1.

### 代码位置

`wald_burden()` in `rvMR/R/wald_burden.R:77` 算三个量. Delta SE helper: `delta_method_ratio_se()` in `utils.R:89`. F helper: `f_statistic()` in `utils.R:133`. **已实现, 跑过 testthat.**

### 陷阱

- 把 delta SE 直接用作 CI 的 SE 在弱 IV regime 下: Gaussian 形状 undercover.
- 加 F > 10 过滤: 把 AR 的鲁棒性优势丢光.
- 漏 cross-term $-2(\hat b_y/\hat b_x^3) \rho \mathrm{SE}_x \mathrm{SE}_y$ 在 sample overlap 情况.
- 引 Bowden-Vansteelandt 2011 *Stat Med* 作为 delta SE source: 那篇是 case-control SMM, **不是** ratio SE 公式.

---

## Step 4 · K=1 closed-form Anderson-Rubin CI (核心难点)

### 直觉 (慢, 多段)

**第一段 — 问题的根**. Wald CI ($\hat\beta \pm 1.96 \cdot \mathrm{SE}$) 的 SE 公式分母里是 $\hat b_x^2$ 或 $\hat b_x^4$. 当 $\hat b_x$ 很小 (rare-variant 典型情况), 这个分母接近 0, SE 爆炸, CI 宽到没意义. 在 $\hat b_x$ 恰好换号附近, CI 还会出现 "中间挖空" 的几何 — 但 Wald 框架根本无法表达这种 CI, 它强行给你一个 misleading 对称区间.

**第二段 — Anderson-Rubin 的转换**. AR (Anderson & Rubin 1949, *Annals of Mathematical Statistics* 20:46) 换了一个角度: *不要估算 $\beta$ 然后建 CI*, 而是 *对每个候选 $\beta_0$ 测试 "$\beta = \beta_0$ 是否兼容数据"*. 在候选 $\beta_0$ 下定义残差 $\hat b_y - \beta_0 \hat b_x$. 如果 $\beta_0$ 等于真值, 这个残差 *在 IV exclusion 假设下* 期望为零 — 不管 IV 多弱. 把残差跟它的方差比, 在 $H_0$ 下服从 $\chi^2_1$. 接受集 $\{\beta_0 : AR(\beta_0) \le 3.84\}$ 就是 95% CI.

**第三段 — 为什么这是弱 IV 鲁棒的**. AR 公式的 *分母* 是 $\mathrm{SE}_y^2 + \beta_0^2 \mathrm{SE}_x^2 - 2\beta_0 \rho \mathrm{SE}_x \mathrm{SE}_y$. 注意 *从未 invert $\hat b_x$*. 不管 $\hat b_x \to 0$, 分母被 $\mathrm{SE}_y^2 > 0$ 兜底, AR(\beta_0) well-defined. 所以 coverage 在 F → 0 时 *仍是* $1 - \alpha$. 这是 "pivotal" 性质.

**第四段 — 四种 CI 形状 (rvSMR 的视觉特色)**. 把 AR($\beta_0$) ≤ c 化简为 $\beta_0$ 的二次不等式 $A \beta_0^2 + B \beta_0 + C \le 0$. 二次开口方向 = sign(A) = sign(F - c); 判别式 $\Delta = B^2 - 4AC$ 决定根存在与否. 四个组合 → 四种 CI 几何:

1. **A > 0, Δ ≥ 0**: 强 IV, 有解, **bounded interval** [a, b]. 跟传统 CI 长得一样.
2. **A > 0, Δ < 0**: 强 IV, 无解, **empty** (CI 为空; 在 $H_0$ 下概率 ≤ α, 是数值 artifact).
3. **A < 0, Δ ≥ 0**: **弱 IV, disconnected union** $(-\infty, a] \cup [b, \infty)$ — 数据排除了中间一段, 但允许 $|\beta|$ 任意大.
4. **A < 0, Δ < 0**: **弱 IV, 无信息**, **whole line** $(-\infty, \infty)$ — 数据对 $\beta$ 没限制.

后两种几何形状是 Wald 框架无法表达的; rvSMR 诚实地报告它们.

### 数学

$$
AR(\beta_0) \;=\; \frac{(\hat b_y - \beta_0 \hat b_x)^2}{\mathrm{SE}_y^2 + \beta_0^2 \mathrm{SE}_x^2 - 2 \beta_0 \rho\, \mathrm{SE}_x \mathrm{SE}_y} \;\xrightarrow{H_0: \beta = \beta_0}\; \chi^2_1.
$$

$$
\mathcal{C}_{1-\alpha} \;=\; \{\beta_0 : AR(\beta_0) \le c\},\quad c = \chi^2_{1, 1-\alpha} = 3.841 \text{ at } \alpha = 0.05.
$$

二次不等式形式:
$$
\underbrace{(\hat b_x^2 - c\,\mathrm{SE}_x^2)}_{A} \beta_0^2 + \underbrace{[-2(\hat b_x \hat b_y - c \rho \mathrm{SE}_x \mathrm{SE}_y)]}_{B} \beta_0 + \underbrace{(\hat b_y^2 - c\,\mathrm{SE}_y^2)}_{C} \le 0,
$$
其中 $A = \mathrm{SE}_x^2 (F - c)$, 所以 sign(A) = sign(F − c).

### 这一步用了哪篇 paper 的什么方法

- **AR 检验的起源**: Anderson TW, Rubin H 1949 *Annals of Mathematical Statistics* 20(1):46. §3, Eq 3 — 在 IV regression 中检验结构参数的检验.
- **AR-K=1 与 Fieller 比 CI 的代数等价**: Fieller EC 1954 *JRSSB* 16(2):175 — Inverted Wald-ratio CI for a ratio of normal means.
- **AR 在 MR 中的现代用法**: Wang S, Kang H 2022 *Biometrics* 78(4):1699 — "Weak-instrument-robust tests in two-sample summary-data Mendelian randomization". 单 IV / 多 IV 的 AR variant 在 MR 设定的系统化.
- **MVMR-AR 延伸 (Step 6 用)**: Patel A, Lane J, Burgess S 2024 arXiv:2408.09868 — 实际标题 "Weak instruments in multivariable Mendelian randomization: methods and practice", *不是* "AR tests for MR" (HANDOVER 早期版本误引).
- **F > 10 obsolete 的现代理论支撑**: Lee-McCrary-Moreira-Porter 2022 *AER* 112(10):3260.

### 为什么不用替代方案

- **Wald CI $\hat\beta_{xy} \pm 1.96 \cdot \mathrm{SE}_\delta$**: 弱 IV 时 SE 爆炸或 CI 几何根本不存在 (e.g. disconnected). Wald 框架强行给一个 misleading 对称区间.
- **Fieller 1954 inverted Wald**: 数学上跟 K=1 AR 等价. 但 Fieller 没有 over-id 推广 — K ≥ 2 时它无能为力, AR 直接推到 $\chi^2_K$ 多 IV 形式.
- **dIVW** (Ye 2021): debias 点估计, 但 F < 3.84 时 CI 仍崩溃.
- **MR-RAPS** (Zhao 2020): profile-score, 在 horizontal pleiotropy 下 robust, 但弱 IV 下 SE 估算不稳.
- **Stock-Yogo F > 10 过滤**: 抛弃弱 IV 数据 — 但 rare-variant 几乎全是弱 IV, 等于抛弃整个数据集.

### 代码位置

`mrAR()` in `rvMR/R/mrAR.R:88`. 二次系数 A, B, C 在 `mrAR.R:115-117`; 几何分类基于 (sign A, sign Δ) 在 `mrAR.R:135-191`. **已实现, 跑过 testthat (34 个断言).**

### 陷阱

- 把 AR CI 和 F > 10 过滤一起用: 把 AR 的鲁棒性优势打掉.
- 把 empty CI 当 "不可能" — 在 $H_0$ 下以 ≤ α 概率发生, 是数值 artifact, 应当 flag 而非压制.
- 漏掉 cross-term $-2 \beta_0 \rho \mathrm{SE}_x \mathrm{SE}_y$ 在 $\rho \neq 0$ 情况.
- 把 disconnected CI 当 bug 强行 re-bound 成 bounded interval: 这等于声称数据支撑了一个数据并不支撑的覆盖率.

---

## Step 5 · Stack K ≥ 3 masks per gene

### 直觉

Step 4 是单 mask 单 IV. Step 5 开始走 *over-identification* 路线: 同一基因的三个 mask (pLoF, missense:LC, regulatory) 提供 *三个独立的 reduced-form moment*. 如果三个 mask 都真正在 instrument 同一个 $X \to Y$ 因果效应, 它们应该 *彼此同意*. 这种 "同意" 就是 Step 8 Sargan-J 检验的来源.

为什么至少 K ≥ 3 而不是 K ≥ 2? 数学上, Sargan-J df = K - 1. K=2 时 df = 1, 一个 outlier mask 就能让 J 通过或拒绝, 几乎无诊断力. K=3 时 df = 2, J 才有 power 区分 "两 mask 同意 vs 三 mask 同意", 能定位异常 mask. 所以 *K ≥ 3 是结构性要求*, 不是 power 优化.

### 数学

$$
\hat{\boldsymbol b}_x = (\hat b_{x,1}, \dots, \hat b_{x,K})^\top, \quad \hat{\boldsymbol b}_y = (\hat b_{y,1}, \dots, \hat b_{y,K})^\top \in \mathbb{R}^K,
$$
$$
D_x = \mathrm{diag}(\mathrm{SE}_{x,1}, \dots, \mathrm{SE}_{x,K}), \quad D_y = \mathrm{diag}(\mathrm{SE}_{y,1}, \dots, \mathrm{SE}_{y,K}),
$$

加 K×K 相关矩阵 $R_{xx}, R_{yy}, R_{xy}$. 默认 $R_{xx} = R_{yy} = I_K$ (mask 间 LD ≈ 0), $R_{xy} = \mathbf{0}$ (两样本无重叠).

### 这一步用了哪篇 paper 的什么方法

- **K ≥ 3 over-id 设计 commitment**: 来自 rvSMR HANDOVER §2 设计选择, 数学根据是 Sargan 1958 *Econometrica* 26:393 的 over-identifying restrictions test.
- **多 IV MR 的现代 AR 框架**: Patel-Lane-Burgess 2024 arXiv:2408.09868 — MVMR setting 下的 AR test, 包含 over-id 的处理.
- **summary-stat MR 多 IV 的 reduced-form moment 解释**: Burgess-Butterworth-Thompson 2013 *Genet Epidemiol* 37(7):658.

### 为什么不用替代方案

- **K=2**: df=1, J test 几乎无 power.
- **K=1 跨细胞类型 stratification 来凑 over-id**: 也可以, 但混杂了"细胞类型是否一致"和"mask 是否一致"两种异质性 — 把 K = mask 数和 C = 细胞类型数 分开是 rvSMR 的设计 (Step 12 处理 cell-type 那一轴).
- **Cauchy-combine 三 mask 之后再跑 MR**: ACAT 把三个 p-value 压成 1 个, 完全摧毁 over-id 结构.

### 代码位置

`mrAR_multi()` in `rvMR/R/mrAR_multi.R:114`; 对角 SE 构造 在 `mrAR_multi.R:153-154`. **已实现, 跑过 testthat.**

### 陷阱

- mask 间共享变体: $R_{xx}, R_{yy}$ 不再对角.
- $K=2$ 时把 J 当真过-id 检验解读.
- 把"mask"理解成 STAARpipeline 默认 mask 集合而非 annotation class.
- Cauchy-combine masks 把 over-id 杀掉.

---

## Step 6 · K ≥ 2 AR statistic and confidence set

### 直觉

把 Step 4 的 scalar AR 推广到向量. 残差 $m(\beta_0) = \hat{\boldsymbol b}_y - \beta_0 \hat{\boldsymbol b}_x \in \mathbb{R}^K$. 它的方差矩阵 $V(\beta_0)$ 是 $K \times K$, 含三项 — outcome 端 SE, exposure 端 SE 按 $\beta_0^2$ 加权, cross-sample 修正项. AR 统计量是 quadratic form $m^\top V^{-1} m$, 在 $H_0: \beta = \beta_0$ 下服从 $\chi^2_K$.

关键直觉: *从未 invert $\hat{\boldsymbol b}_x$*. AR 依赖 *reduced-form* moment, 不依赖 *first-stage* estimator. 所以弱 IV 鲁棒性从 K=1 直接继承到 K≥2.

每个 mask 自己的第一阶段 F 都比 *pooled* burden 弱, 这看起来是坏事 — 但 AR 把这 K 个弱信号合到一个 $\chi^2_K$, 既保留 uniform coverage, 又拿到 K-1 个 over-id 自由度. 这是 stacking 的 *双赢*.

### 数学

$$
m(\beta_0) \;=\; \hat{\boldsymbol b}_y - \beta_0 \hat{\boldsymbol b}_x \in \mathbb{R}^K
$$
$$
V(\beta_0) \;=\; D_y R_{yy} D_y + \beta_0^2 D_x R_{xx} D_x - 2 \beta_0 D_y R_{xy} D_x \in \mathbb{R}^{K \times K}
$$
$$
\boxed{AR(\beta_0) \;=\; m(\beta_0)^\top V(\beta_0)^{-1} m(\beta_0) \;\xrightarrow{H_0}\; \chi^2_K}
$$
$$
\mathcal{C}_{1-\alpha} = \{\beta_0 : AR(\beta_0) \le \chi^2_{K, 1-\alpha}\}.
$$

### 这一步用了哪篇 paper 的什么方法

- **多 IV AR 的统计构造**: Anderson-Rubin 1949 *AMS* 20:46 是起源; 多 IV 的 quadratic form 形式来自 Wang-Kang 2022 *Biometrics* 78:1699 §2.2.
- **MVMR setting 的 AR 推广**: Patel-Lane-Burgess 2024 arXiv:2408.09868 — "Weak instruments in multivariable Mendelian randomization: methods and practice".

### 为什么不用替代方案

- **直接对每个 mask 跑 K=1 AR 然后 meta-analyze**: 损失 over-id 信息 (meta-analysis 只能 pool point estimates, 不能 jointly test 同质性).
- **IVW (Inverse-Variance Weighted)** (Burgess-Butterworth-Thompson 2013): point + Wald SE; 弱 IV 下崩溃.
- **2SLS**: 需要 individual-level data; 即便用 summary-stat 近似, 弱 IV 下 size 偏差.
- **MR-Egger** (Bowden 2015 *Int J Epidemiol* 44:512): pleiotropy intercept; 但弱 IV 下 SE 不稳.

### 代码位置

`mrAR_multi()` in `rvMR/R/mrAR_multi.R:114`. $V_{yy}, V_{xx}, V_{xy}$ 在 `mrAR_multi.R:155-157`; AR closure $\beta_0 \mapsto m^\top V^{-1} m$ 在 `mrAR_multi.R:162-171`. **已实现.**

### 陷阱

- 漏 cross-term $-2 \beta_0 D_y R_{xy} D_x$ 在 sample overlap.
- 顺序写反: 应是 $D_y R_{xy} D_x$ (即 (i,j) entry 是 $\mathrm{SE}_{y,i} R_{xy,ij} \mathrm{SE}_{x,j}$), *不是* $D_x R_{xy} D_y$.
- 用 $\chi^2_{K-1}$ 作 AR 的参考分布 (那是 J 的, Step 8). AR 自己是 $\chi^2_K$.
- $V(\beta_0)$ 可能在大 $|\beta_0|$ + 高 overlap 下变 indefinite, AR 在该点未定义, 需 flag.

---

## Step 7 · Grid + uniroot inversion of the AR confidence set

### 直觉

K=1 时, $AR(\beta_0) \le c$ 是 $\beta_0$ 的二次不等式 — 有闭式根. K ≥ 2 时, $V(\beta_0)$ 是 $\beta_0$ 的矩阵值二次函数, 它的逆是 *有理函数*, AR 不再有闭式根. 我们必须 *数值* 求 level set.

策略很标准:
1. 选一个初始扫描范围 $[\beta_\min, \beta_\max]$ (取 per-IV Wald 比的 min/max + padding).
2. 在 4000 个 grid 点上算 $AR(\beta_0) - c$.
3. 找 sign change (相邻两 grid 点异号).
4. 每个 sign change 用 `uniroot()` 精修到 1e-8.
5. 用精修 roots 把实轴切段, 段内取中点判 accept / reject.

输出: CI shape label ∈ {bounded_interval, disconnected_union, whole_line, empty} + endpoints list. K=1 的四分类几何在这里 *自然继承*.

### 数学

不是新公式, 而是算法:

$$
\text{Find } \{\beta_0 : AR(\beta_0) - c = 0\} \text{ via grid + uniroot};\quad \mathcal{C}_{1-\alpha} = \bigcup_{(l_i, u_i) \text{ accepted}} [l_i, u_i].
$$

### 这一步用了哪篇 paper 的什么方法

- **数值 level-set 求解 (grid + uniroot)**: 标准统计计算实践. R 的 `stats::uniroot()` 实现 Brent 1973 的 root-finding algorithm.
- **AR CI 数值反演的统计基础**: Anderson-Rubin 1949 给出统计量定义; 数值反演 (CI 几何分类) 是工程实现, 不需要额外 paper anchor.
- **CI 几何分类来源**: Fieller 1954 *JRSSB* 16:175 在 K=1 下给出四种几何, rvSMR 把它推到 K≥2 数值地继承.

### 为什么不用替代方案

- **Closed-form (K=1)**: K ≥ 2 不存在.
- **Profile likelihood maximization**: 给的是 argmin (Step 8 用), 不是 level set 边界.
- **Bootstrap CI**: rare-variant + 弱 IV 下 bootstrap 也崩溃; 重采样不能修复 first-stage 接近 0 的问题.

### 代码位置

`mrAR_multi()` in `rvMR/R/mrAR_multi.R:114`. Envelope 在 `mrAR_multi.R:175-189`; grid extension 在 `mrAR_multi.R:201-220`; sign-change 在 `mrAR_multi.R:222-225`; uniroot 在 `mrAR_multi.R:228-242`; 分类 在 `mrAR_multi.R:246-319`. **已实现.**

### 陷阱

- Grid 太粗: 两个挨近 root 漏掉, disconnected_union 被误判 bounded_interval.
- Envelope 太窄: 真无界 ray 被截断成有限 endpoint; 必须靠 extension loop 救.
- 把 uniroot 失败 (返 NA) 当 root 接受: 实现里 `mrAR_multi.R:243` 显式过滤非有限 root.
- $V(\beta_0)$ singular 子区间不 inspect: 算法 mark +∞, 但用户应看 grid trace.

---

## Step 8 · Sargan-J over-identification at the AR argmin

### 直觉

如果三个 mask 都真正在 instrument 同一个 $\beta$, 它们应该 "在某个 $\hat\beta_{\rm AR}$ 处同时让残差变小". 残差 *都*小, AR 最小值 $\min_{\beta_0} AR(\beta_0)$ 就小 — 在 $H_0$ 下服从 $\chi^2_{K-1}$. 大的 J 说明 "没有一个 $\beta$ 能让三 mask 都 happy" → 至少 1 个 mask 在违反 IV 假设.

为什么 df = K-1 而不是 K? 直观: $AR(\beta_0)$ 在固定 $\beta_0$ 下是 K 维 $\chi^2$. 取 min 等于估了 1 个参数 ($\hat\beta_{\rm AR}$), 类比 OLS 的 $\chi^2_{n-p}$, df 减 1. 所以 J $\sim \chi^2_{K-1}$.

J 是 *联合* 检验: "至少一个 mask 错了". 不告诉你 *哪个* mask 错了. 定位需要 Step 10 (HEIDI-rv, mask 内) 和 Step 11 (annotation Q, mask 间).

### 数学

$$
\hat\beta_{\rm AR} \;=\; \arg\min_{\beta_0}\, AR(\beta_0)
$$
$$
J \;=\; AR(\hat\beta_{\rm AR}) \;\xrightarrow{\text{joint IV validity}}\; \chi^2_{K-1}
$$

P-value: $1 - F_{\chi^2_{K-1}}(J)$. Reject homogeneity if $J > \chi^2_{K-1, 1-\alpha}$.

### 这一步用了哪篇 paper 的什么方法

- **Sargan 1958 *Econometrica* 26(3):393** — over-identifying restrictions test 的起源 ("The Estimation of Economic Relationships Using Instrumental Variables").
- **Hansen 1982 *Econometrica* 50(4):1029** — GMM J-statistic 的现代 (generalized) 推广. rvSMR 的 J 是这种 GMM 形式在 two-sample summary MR 设定下的特例.
- **MR setting 的 AR-J**: Patel-Lane-Burgess 2024 arXiv:2408.09868. 给出 MR 多 IV 下的 J.
- **df = K-1 specifically**: 直接来自 Sargan 1958 推导 (over-id 自由度 = #IVs − #endogenous parameters; 这里 endogenous parameter 只有 $\beta$, 所以 = K - 1).

### 为什么不用替代方案

- **Cochran's Q on Wald ratios** (Bowden 2017): 直接对 K 个 mask 的 Wald 比做 inverse-variance Q. 但 Q 依赖 Wald 比的 delta SE, 弱 IV 下崩溃.
- **MR-Egger intercept test** (Bowden 2015 *Int J Epidemiol* 44:512): 测试 pleiotropy 的均值偏 (directional pleiotropy). 但 MR-Egger 假设 InSIDE (instrument strength independent of direct effect), 罕见变体下不稳.
- **MR-PRESSO** (Verbanck 2018 *Nat Genet*): outlier removal — 但 K=3 时移除 1 个 = 33% 数据, 太激进.

### 代码位置

`mrAR_multi()` in `rvMR/R/mrAR_multi.R:114`. `stats::optimize()` argmin + grid fallback 在 `mrAR_multi.R:321-335`; J p-value (df = K-1) 在 `mrAR_multi.R:337-341`. **已实现.**

### 陷阱

- 报 J $\sim \chi^2_K$: 自由度错, 反保守, 漏检 pleiotropy.
- 把小 J 解读成 "没 pleiotropy": 实际是 "没 *跨 mask* 异质性"; coherent pleiotropy (Step 14) 无法被 J 检测.
- 用别的 argmin (e.g. delta-method midpoint) 算 J: J 特指 min AR.
- K=1 时强行报 J: 不存在 over-id (df = 0), 函数返 NA.

---

## Step 9 · Sample-overlap correction via the AR cross-term

### 直觉

如果 exposure 端和 outcome 端有 *样本重叠* (e.g. 都用 UKB), 两端的 $\hat b_x$ 和 $\hat b_y$ 不再独立, 它们的 sampling error 有 correlation $\rho \neq 0$. AR 公式分母里要加 cross-term 来吃掉这层相关.

关键洞察: cross-term 跟 $\beta_0$ 成线性, **而且符号会翻**. 当 $\beta_0$ 跨过 0, cross-term 反号 — 同一个 $\rho$ 在 $\beta_0 > 0$ 端膨胀 CI, 在 $\beta_0 < 0$ 端收缩 CI. 这让 AR CI 关于 0 *不对称*, 即使 $\rho$ 是常数.

默认两样本 ($\rho = 0$ / $R_{xy} = \mathbf{0}$) 时该项消失.

### 数学

K=1:
$$
\mathrm{denom}(\beta_0) = \mathrm{SE}_y^2 + \beta_0^2 \mathrm{SE}_x^2 \boxed{- 2\beta_0 \rho \mathrm{SE}_x \mathrm{SE}_y}.
$$

K ≥ 2:
$$
V(\beta_0) = D_y R_{yy} D_y + \beta_0^2 D_x R_{xx} D_x \boxed{- 2 \beta_0 D_y R_{xy} D_x}.
$$

### 这一步用了哪篇 paper 的什么方法

- **Sample overlap 修正在两样本 MR 中的处理**: Burgess S, Davies NM, Thompson SG 2016 *Genetic Epidemiology* 40(7):597 — "Bias due to participant overlap in two-sample Mendelian randomization". **注意: 是 Burgess et al. 2016, *不是* Burgess-Butterworth-Thompson 2013** (后者是 summary-stat MR 框架的 base paper).
- **Cross-term 推导**: 标准 IV 变量代数; 在 AR 框架下来自 Anderson-Rubin 1949 + Wang-Kang 2022 §2.3.

### 为什么不用替代方案

- **LDSC intercept correction** (Bulik-Sullivan 2015 *Nat Genet*): 用 cross-trait LD score regression intercept 估 sample overlap. 但 intercept 把 *overlap* 和 *cross-trait genetic covariance* 混淆, 在 rare-variant 下尤其不可靠.
- **Re-running on disjoint subsamples**: 解决 overlap 但损失 50% 样本量.
- **Ignoring overlap (设 $\rho = 0$)**: silent type-I error inflation 当 $\rho > 0$ 且 $\beta_0$ 与 $\hat\beta$ 同号.

### 代码位置

K=1: `mrAR.R:116` 的 B 系数 $-2(\hat b_x \hat b_y - c \rho \mathrm{SE}_x \mathrm{SE}_y)$. K≥2: `mrAR_multi.R:157` 的 $V_{xy} = D_y R_{xy} D_x$, 在 `mrAR_multi.R:164` propagate 到 closure. **已实现.**

### 陷阱

- 一样本分析默认 $\rho = 0$: 常见 silent bug.
- 顺序写反 $D_x R_{xy} D_y$ 而非 $D_y R_{xy} D_x$.
- LDSC intercept 估 $R_{xy}$ 不查重叠比例.
- 忘记 cross-term sign-flip 让 CI 关于 0 不对称.

---

## Step 10 · HEIDI-rv within-burden heterogeneity

### 直觉

Sargan-J (Step 8) 检验 *mask 间* 同质性; HEIDI-rv 检验 *mask 内* 同质性. 也就是: 在同一 mask 里的 m 个 variant, 它们的 *per-variant Wald 比* $\hat\beta_{xy}^{(j)} = \hat b_{y,j} / \hat b_{x,j}$ 应该彼此一致 — 如果 mask 的 IV 假设成立, 每个 variant 估的都是同一个 $\beta$.

经典 HEIDI (Zhu 2016 *Nat Genet* 48:481, common-variant SMR) 是给同一基因附近 LD-相关的 SNP 设计的. Rare-variant 没有 LD partner — HEIDI-rv 必须重写为 *contrast test*: 选一个零均 contrast 矩阵 C, 算 $\delta = C \hat{\boldsymbol\beta}_{xy}^{\rm per-var}$, 它的方差 $V_\delta = C \Sigma_{\hat\beta} C^\top$. 检验 $\delta = 0$ — 即 m 个 per-variant Wald 比互相相等.

关键统计细节: $V_\delta$ 通常 *rank-deficient* (有 m-1 行但秩 ≤ m-1). 用 Moore-Penrose 伪逆 $V_\delta^+$ 写 $T = \delta^\top V_\delta^+ \delta$. 在 $H_0$ 下 T 服从 *广义 $\chi^2$* — 独立 $\chi^2_1$ 的加权和, 权重是 $V_\delta$ *本身* 的非零特征值. **不是** $V_\delta^+ V_\delta$ 的特征值 (那是 projector, 特征值 0/1, 会把广义 $\chi^2$ 退化成 plain $\chi^2$, 反保守).

尾概率用 Davies 1980 算法 (精确, 用 CompQuadForm::davies); 在尾部可以用 Kuonen 1999 鞍点近似.

诚实地说: HEIDI-rv 在 mask 大小 m 大时 *power 差*. 单个 pleiotropy variant 在 burden 内只贡献 $\mathcal{O}(1/m)$ 信号, m=50 时几乎检测不到. 所以 rvSMR 把它定位为 *sanity check*, 不是 primary pleiotropy test.

### 数学

$$
\delta \;=\; C \hat{\boldsymbol\beta}_{xy}^{\rm per-var} \in \mathbb{R}^{m-1}, \qquad V_\delta = C \Sigma_{\hat\beta} C^\top \in \mathbb{R}^{(m-1) \times (m-1)}.
$$
$$
T \;=\; \delta^\top V_\delta^+ \delta \;\sim\; \sum_{i=1}^{\mathrm{rank}(V_\delta)} \lambda_i \chi^2_{1,i}, \quad \boxed{\lambda_i = \text{非零特征值 of } V_\delta \text{ 本身}}
$$
$$
\Pr(T \ge t \mid H_0) \;=\; \texttt{CompQuadForm::davies}(q = t, \lambda = \mathrm{eig}_{\ne 0}(V_\delta)).
$$

### 这一步用了哪篇 paper 的什么方法

- **HEIDI 起源 (common-variant SMR)**: Zhu Z et al. 2016 *Nature Genetics* 48(5):481 — "Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets". HEIDI = HEterogeneity In Dependent Instruments.
- **广义 $\chi^2$ 尾概率算法**: Davies RB 1980 *Applied Statistics (JRSS-C)* 29(3):323 — "Algorithm AS 155: the distribution of a linear combination of $\chi^2$ random variables".
- **Saddlepoint approximation 备用**: Kuonen D 1999 *Biometrika* 86(4):929 — "Saddlepoint approximations for distributions of quadratic forms in normal variables".
- **rvSMR 把 HEIDI 适配到 rare-variant**: 新方法 (within-burden contrast + Moore-Penrose), Step 10 描述的就是 rvSMR 的 contribution.

### 为什么不用替代方案

- **LD-based HEIDI 原版**: rare variant 没有 LD partner, 不可用.
- **MR-PRESSO outlier removal**: 适合 common-variant 多 SNP MR, rare-variant burden 内不便.
- **Single-variant leave-one-out**: 等价于 HEIDI-rv 的特殊 contrast (每次删 1 个 variant), 但漏掉对全 contrast space 的检验.

### 代码位置

`heidi_rv()` in `rvMR/R/heidi_rv.R:102` — 🔴 **stub** (`stop("not implemented")` at line 103, 等 Team 1 实现). 应当传 `lambda = eig_nz(V_delta)` 给 `CompQuadForm::davies`.

### 陷阱

- **headline trap**: 把 $V_\delta^+ V_\delta$ 的特征值传给 Davies. $V_\delta^+ V_\delta$ 是 projector, 特征值 {1, ..., 1, 0, ..., 0}, 会把 T 的 reference 塌成 plain $\chi^2_{\mathrm{rank}(V_\delta)} = \chi^2_{m-1}$ — 反保守 null, 多检测 pleiotropy.
- Zero eigenvalues 不丢: monomorphic-in-sample variants 让 $V_\delta$ rank-deficient, 0-weight $\chi^2_1$ component 数值上仍贡献 df.
- 当 primary pleiotropy test 用: $\mathcal{O}(1/m)$ 不敏感性.
- 期望从 STAARpipeline aggregate output 跑 HEIDI-rv: 需要 SAIGE-QTL 和 Genebass single-variant mode 重跑 + LD 矩阵.

---

## Step 11 · Annotation-class concordance via pQTL-anchor normalization

### 直觉

Step 8 的 Sargan-J 假设三 mask 估的都是同一个 $\beta$. 但 *raw* burden Wald 比在三 mask 间天然不可比: pLoF 把蛋白完全打掉, missense:LC 部分损失, regulatory 改变表达水平 — 三类 mask 把 *burden 映射到 protein abundance* 的 effect 完全不同. 同一个下游 $X \to Y$ 因果效应, 在三 mask 上会被放大不同的 mediator scale.

解决: 用 pQTL anchor (UKB-PPP rare-variant burden-on-protein effect $\hat\beta^{(k)}_{\rm burden \to protein}$) 做 mediator-scale normalization. 把每个 class 的 burden-on-outcome 除以它的 burden-on-protein, 得到 *per-protein-unit* slope $\tilde\beta^{(k)}_{xy}$. 现在三 class 在同一尺度上, 可以做 Cochran's Q 检验同质性.

大 Q → mechanism 异质: 可能是 dominant-negative (pLoF vs missense 效应反号), gain-of-function (missense 异常激活), regulatory-only effect — *flag for mechanism investigation*, 不是自动 reject.

### 数学

$$
\tilde\beta^{(k)}_{xy} \;=\; \frac{\hat\beta^{(k)}_{\rm burden \to y}}{\hat\beta^{(k)}_{\rm burden \to protein}}
$$
$$
Q \;=\; \sum_{k=1}^K w_k (\tilde\beta^{(k)}_{xy} - \bar{\tilde\beta})^2 \;\xrightarrow{H_0}\; \chi^2_{K-1},
\quad w_k = 1/\mathrm{Var}(\tilde\beta^{(k)}_{xy}),\ \bar{\tilde\beta} = \frac{\sum_k w_k \tilde\beta^{(k)}_{xy}}{\sum_k w_k}.
$$

### 这一步用了哪篇 paper 的什么方法

- **Cochran's Q 同质性检验**: Cochran WG 1954 *Biometrics* 10(1):101 — "The combination of estimates from different experiments".
- **UKB-PPP rare-variant pQTL anchor**: Dhindsa RS et al. 2023 *Nature* 622:339 — "Rare variant associations with plasma protein levels in the UK Biobank".
- **deCODE pQTL alternative**: Ferkingstad E et al. 2021 *Nature Genetics* 53:1712 — "Large-scale integration of the plasma proteome with genetics and disease".
- **不要引: UKB-PPP common-variant Sun BB 2023 *Nature* 622:329** (这是同期但不同 paper, common-variant proteomics; rvSMR Step 11 需要 *rare-variant* burden-on-protein).

### 为什么不用替代方案

- **Skip normalization** (直接对 raw $\hat\beta^{(k)}_{\rm burden, y}$ 做 Q): raw 三 class 因 mediator scale 差距大, 总是 disagree. Q 永远 reject.
- **Sun BB 2018 INTERVAL pQTL** (Sun BB 2018 *Nature* 558:73) 或 **Sun BB 2023 UKB-PPP common-variant pQTL**: 是 *common-variant* pQTL, scale 转换不对 — 应当用同样 mask 的 *rare-variant burden-on-protein*.
- **Han-Eskin 2011 *AJHG* 88:586** (sign concordance): 那篇是 random-effects meta-analysis, *不是* sign-concordance test. 不要引.

### 代码位置

`annotation_concord()` in `rvMR/R/annotation_concord.R:85` — 🔴 **stub** (`stop("not implemented")` at line 86). Roxygen 记 `estimates_list` schema 含 `b_burden_to_protein`, `se_burden_to_protein`.

### 陷阱

- Skip pQTL normalization.
- 用 common-variant pQTL 而不是 rare-variant burden-on-protein.
- 引 Han-Eskin 2011.
- 把小 Q 读成 "no pleiotropy" — coherent pleiotropy 通过同上游 confounder 同时影响三 class, Q 不动.

---

## Step 12 · Cell-type concordance

### 直觉

第三个 over-id 轴: cell type. 同一个 (gene, mask) 在 C=28 个免疫 cell type 里跑出 C 个 burden Wald 比. 如果 $\beta$ 是基因层面性质, 它们应该都同意 — Cochran Q 检验.

异质性可能是 *feature* (e.g. hepatocyte-specific PCSK9 → LDL, 是真生物学) 也可能是 *bug* (cell-type-specific pleiotropy). Q 本身不区分, 需要 mechanism follow-up.

注意: TenK10K Phase 1 只有 PBMC (28 个免疫亚型). Hepatocyte / adipocyte 要等 Phase 2.

### 数学

$$
Q_{\rm cell} \;=\; \sum_{c=1}^C w_c (\hat\beta^{(c)}_{xy} - \bar\beta)^2 \;\xrightarrow{H_0}\; \chi^2_{C-1}, \quad w_c = 1/\mathrm{Var}(\hat\beta^{(c)}_{xy}).
$$

或完整 stratified AR: 把 K mask × C cell type 一起 stack 成 $K\cdot C$ 个 IV, 跑 Step 6 的 multi-IV AR.

### 这一步用了哪篇 paper 的什么方法

- **TenK10K Phase 1 substrate (28 PBMC immune subsets)**: Cuomo ASE et al. 2025 medRxiv 2025.03.20.25324352 — "TenK10K Phase 1: cell-type-specific eQTLs in 28 PBMC immune cell types".
- **Cell-type cis-MR comparator (first published method)**: Ray D et al. 2025 *American Journal of Human Genetics* 112(7):1597 — single-cell cis-MR across 14 immune cell types (common-variant). **不是 Ge 2025** (那是 confabulation).
- **Cochran Q (同 Step 11)**: Cochran 1954 *Biometrics* 10:101.

### 为什么不用替代方案

- **OneK1K** (Yazar 2022 *Science* 376:eabf3041): 14 immune subsets, 较旧, 样本量小. TenK10K 是 OneK1K 升级版.
- **Bulk eQTL (GTEx)**: 没有 cell-type resolution, 跟 single-cell substrate 的精度不匹配.
- **Cell-type pseudobulk + bulk MR**: 损失 cell-type stratification 的 over-id 信息.

### 代码位置

不是单独函数; cell-type 轴是高层 wrapper, loops `wald_burden()` (`wald_burden.R:77`) + `mrAR()` (`mrAR.R:88`) over cell types, 再用 Cochran-Q (类似 `annotation_concord()` 的机器). **wrapper 待实现 (Team 1 路线图).**

### 陷阱

- 充 C 太多: cell type 里只有几个 donor 携带 mask variant, SE 不稳, Q 被噪声主导.
- 声称 hepatocyte 分辨率从 TenK10K Phase 1 跑出来: Phase 1 是 PBMC, Phase 2 才有 hepatocyte.
- 引 "Ge 2025" 作 sc-cis-MR comparator: 错的, 应当是 **Ray 2025 *AJHG* 112(7):1597**.
- 把 cell-type Q 和 annotation Q 用 Fisher 合并而不考虑 shared (g) axis: 不独立.

---

## Step 13 · Sensitivity scalars (Cinelli-Hazlett and E-value)

### 直觉

前 12 步都假设 IV exclusion restriction 成立 (no confounder 联通 IV → Y 绕过 X). 但 *coherent pleiotropy* — 三 mask 都通过同一上游 confounder 影响 Y — Sargan-J + HEIDI-rv + class Q + cell-type Q 都看不见.

最后的防线是 *sensitivity analysis*: 计算 "要破坏当前结论, unobserved confounder 至少要多强". 两个 metric:

1. **Robustness value RV** (Cinelli-Hazlett): 最小的 partial $R^2$, 让一个 unobserved 违反 exclusion 的 confounder 同时跟 IV 和 outcome 关联, 把效应推到 null.
2. **E-value** (VanderWeele-Ding, Swanson-VanderWeele MR adaptation): 在 risk-ratio 尺度上, 最小的 confounder-arm association 强度可以把 $\hat\beta_{xy}$ 解释掉.

高 RV + 高 E → effect 难以被 confounder 解释; 低 → fragile.

诚实表态: 这是 *honest limit, not a fix*. Sensitivity 量化弱点, 不修复弱点 — Wang-Tchetgen 2018 *JRSS-B* 80:531 的 partial-identification bounds 是更强但更保守的工具.

### 数学

IV partial $R^2$ + RV:
$$
R^2_{Z\to X} \;\approx\; \frac{t^2}{t^2 + n_x - 2}, \quad t = \frac{\hat b_x}{\mathrm{SE}_x}
$$
$$
\boxed{RV \;=\; \frac{\sqrt{t^2 + 4} - t}{2}}
$$

E-value (continuous outcome via Swanson-VanderWeele's RR approximation):
$$
\mathrm{RR} \;\approx\; \exp(0.91 \cdot \beta_{\rm std}), \quad E \;=\; \mathrm{RR} + \sqrt{\mathrm{RR}(\mathrm{RR} - 1)}.
$$

报 $E_{\rm point}$ (at point estimate) 和 $E_{\rm CI}$ (at CI bound nearest null) — 后者更保守, 是标准报告值.

### 这一步用了哪篇 paper 的什么方法

- **PRIMARY (IV setting)**: Cinelli C, Hazlett C 2025 *Biometrika* asaf004 — "An omitted variable bias framework for sensitivity analysis of instrumental variables". 这是 IV-specific 的 sensitivity 框架.
- **SECONDARY (OLS origin)**: Cinelli C, Hazlett C 2020 *JRSS-B* 82(1):39 — "Making Sense of Sensitivity: Extending Omitted Variable Bias". OLS 框架的 partial-$R^2$ / RV 起源.
- **MR E-value adaptation**: Swanson SA, VanderWeele TJ 2020 *Epidemiology* 31(3):e23 — "E-values for Mendelian randomization".
- **E-value 起源 (OLS observational)**: VanderWeele TJ, Ding P 2017 *Annals of Internal Medicine* 167(4):268 — "Sensitivity analysis in observational research: introducing the E-value".

### 为什么不用替代方案

- **MR-Egger intercept** (Bowden 2015): 测 directional pleiotropy, 但需要 ≥ 3 IV 且 InSIDE 假设.
- **Mendelian randomization with bidirectional effects** (Hemani 2017): direction test, 不是 sensitivity bound.
- **Partial-identification bounds** (Wang-Tchetgen 2018 *JRSS-B* 80:531): 更保守, 给 [lower, upper] 区间. 是 Step 14 的辅助工具, Step 13 用 RV/E 是更简洁的 sensitivity scalar.

### 代码位置

`iv_partial_r2()` in `rvMR/R/sensitivity.R:48` — 🔴 **stub** (`stop("not implemented")` at line 49). `e_value()` in `rvMR/R/sensitivity.R:99` — 🔴 **stub** (`stop("not implemented")` at line 100).

### 陷阱

- **Sign error**: RV 应是 $(\sqrt{t^2+4} - t)/2$ (positive). 反号 $(t - \sqrt{t^2+4})/2$ 是负数, 无意义.
- 只引 2020 *JRSS-B* OLS 而不引 2025 *Biometrika* IV: 2025 是 PRIMARY IV 参考, 2020 是 OLS 起源.
- RR 近似系数错: Swanson-VanderWeele 2020 fix 0.91, 不要随便改.
- 只报 $E_{\rm point}$: 标准是 $E_{\rm CI}$.

---

## Step 14 · Gene-level decision rule

### 直觉

把六根诊断 (AR CI, J, HEIDI-rv, class Q, cell Q, RV/E) 合到一个 gene-level 标签. 决策树是 *合取式* (所有 over-id 轴必须通过) + *保守* (coherent pleiotropy 下默认 inconclusive).

标签三类:
- **causal**: AR CI 排除 0, J 不显著, 每 mask 的 HEIDI-rv 不显著, class Q 不显著, cell Q 不显著, $RV, E_{\rm CI}$ 高于阈值.
- **inconclusive**: AR CI 包含 0 (whole_line / 包含 0 的 bounded) 或 sensitivity 低于阈值.
- **pleiotropy-suspect**: J / HEIDI-rv / class Q 任一显著.

诚实表态: 即便所有轴通过, **coherent pleiotropy 仍然能逃过**. 三 mask 都通过同一 confounder → 每 mask 内同质 (HEIDI-rv 通过), mask 间一致 (Sargan-J + class Q 通过), 每 cell type 一致 (cell Q 通过). 唯一防线是 RV/E 量化 confounder 强度 + Wang-Tchetgen 2018 的 partial-identification bounds.

### 数学

决策树 (D1 - D6):

1. **D1 — AR CI 排除 0?** 若 bounded / disconnected 且不含 0, → 非空信号; whole_line → inconclusive; empty → numerical artifact, refine grid.
2. **D2 — J 不显著?** 若 $J > \chi^2_{K-1, 1-\alpha}$ → pleiotropy-suspect.
3. **D3 — HEIDI-rv 不显著?** 任一 mask reject → 降级 pleiotropy-suspect.
4. **D4 — Class Q 不显著?** Reject → pleiotropy-suspect / mechanism-investigate.
5. **D5 — Cell Q 不显著?** Reject → cell-type-specific note, judge biology.
6. **D6 — Sensitivity floor?** $RV < 0.05$ 或 $E_{\rm CI} < 1.5$ → inconclusive (fragile).

只有 D1-D6 全过, 才报 **causal**.

### 这一步用了哪篇 paper 的什么方法

- **Decision-rule 合成**: 是 rvSMR 的设计选择, 综合前 13 步. 没单一 paper 出处.
- **Coherent pleiotropy 警告 + partial-identification fallback**: Wang L, Tchetgen Tchetgen EJ 2018 *JRSS-B* 80(3):531 — "Bounded, Efficient and Multiply Robust Estimation of Average Treatment Effects Using Instrumental Variables". 在 IV invalid 下给出 [lower, upper] 区间. 是 Step 14 的诚实退路.
- **Pleiotropy 的传统对照**: Bowden J et al. 2015 *International Journal of Epidemiology* 44(2):512 — "Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression". MR-Egger 是 directional pleiotropy 的传统检测; rvSMR 的多轴 over-id + sensitivity 是更新方法.
- **Honest-limit framing**: Cinelli-Hazlett 2025 *Biometrika* asaf004 的 framing — sensitivity quantifies, doesn't fix.

### 为什么不用替代方案

- **Single-test paradigm** (e.g. 只看 AR CI): 忽视 over-id 信息, 当 J 反 alarm 时仍 report.
- **Voting across methods** (e.g. IVW + MR-Egger + Median-IVW vote majority): 每个 method 的 weak-IV failure mode 不同, vote 不解决 coherent pleiotropy.
- **Bayesian model average** (RARE 用类似思路): credible interval 不给 frequentist coverage 保证, 监管口径下不被接受.

### 代码位置

决策树不在单一函数里, 是 calling convention: 组合 `wald_burden()` (`wald_burden.R:77`), `mrAR()` (`mrAR.R:88`), `mrAR_multi()` (`mrAR_multi.R:114`), `heidi_rv()` (`heidi_rv.R:102` 🔴 stub), `annotation_concord()` (`annotation_concord.R:85` 🔴 stub), `iv_partial_r2()` (`sensitivity.R:48` 🔴 stub), `e_value()` (`sensitivity.R:99` 🔴 stub).

### 陷阱

- 单一轴作 sufficient: 干净 AR CI + 失败 J 不是 causal call.
- Non-rejection 当 proof of validity: coherent pleiotropy 看不见.
- Stock-Yogo F > 10 pre-filter: 把 AR 整套设计的优势废掉.
- 报 $\hat\beta_{\rm AR}$ 不带 CI shape label: disconnected/whole_line 携带的信息被埋掉.

---

# Intro draft

rvSMR is a rare-variant Summary-data Mendelian Randomization framework. It extends classical SMR (Zhu et al. 2016) by replacing single-SNP Wald inference with Anderson-Rubin (AR) weak-IV-robust confidence sets on linear weighted burden instruments, with three over-identification axes (mask, within-burden, cell-type) plus sensitivity scalars (RV, E-value).

This walkthrough takes you slowly through all 14 algorithmic steps, explaining intuition before formalism, citing the paper that originated each method and the paper that adapted it to MR / rare variants, and listing both why we chose what we chose and why we did not choose the obvious alternatives. Every step has a code reference into the `rvMR` R package.

---

# Summary draft

The pipeline composes:
- Steps 0-2: assemble + construct burden instruments (Madsen-Browning, STAAR, SAIGE-QTL/SAIGE-GENE+).
- Steps 3-4: per-mask Wald ratio + closed-form AR (single IV).
- Steps 5-8: multi-IV AR + Sargan-J over-identification (K≥3 commitment).
- Step 9: sample-overlap cross-term.
- Steps 10-12: three over-id axes (HEIDI-rv within-burden, annotation Q across class, cell-type Q across c).
- Step 13: sensitivity (Cinelli-Hazlett RV, Swanson-VanderWeele E).
- Step 14: gene-level decision tree.

Honest limits: AR has uniform coverage but the CI can be unbounded; J / Q tests do not see coherent pleiotropy; sensitivity quantifies but does not fix; partial-identification bounds (Wang-Tchetgen 2018) are the most conservative fallback.

End of content draft.
