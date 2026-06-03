# Steps 5–9 · 逻辑梳理 (rvSMR AR inference 核心)

*Companion to `main.tex` §§5–9. Source artifacts: `main.tex:190-336`, `mrAR_multi.R` (360 lines), `mrAR.R` (200 lines), `briefing_for_wei.md` §4, `HANDOVER_2026-05-27.md` §§2,10.*

---

## 总览: 这 5 步在做什么 (one-pager)

Steps 5–9 是 rvSMR 把 per-mask burden summary statistics 转成 **weak-IV-robust 置信集 + over-identification 诊断** 的核心机器. 在 Step 4 (K=1 闭式) 之后, Step 5 把信息从 scalar 推到 vector; Step 6 写下多 IV AR 统计量; Step 7 用 grid + uniroot 数值反演置信集 (K=1 的四分类几何在这里自然继承); Step 8 在 AR argmin 处给 Sargan-J 检验, df = K−1 (over-id 的 1 个自由度被 β̂ 吃掉); Step 9 是个"横向"补丁 — 把 sample overlap 通过 cross-term $-2\beta_0 D_y R_{xy} D_x$ 接进 $V(\beta_0)$, 两样本默认 $R_{xy}=0$.

```
Step 3 (per-mask Wald: b_x,SE_x,b_y,SE_y)
                 │
                 ▼
    Step 5: stack K masks into vectors
       (b_x, b_y ∈ R^K; D_x, D_y diag; R_xx, R_yy = I_K; R_xy = 0)
                 │
                 ▼
    Step 6: AR(β₀) = m(β₀)ᵀ V(β₀)⁻¹ m(β₀) ~ χ²_K
       │                                          │
       ▼                                          ▼
    Step 7: invert {β₀ : AR ≤ c_crit}     Step 8: J = min AR ~ χ²_{K−1}
       (grid + uniroot, 4 CI shapes)        (Sargan over-id test)
       │                                          │
       └─── Step 9 enters both via V(β₀): ────────┘
            R_xy ≠ 0  ⇒  − 2 β₀ D_y R_xy D_x cross-term
```

三个关键转变:
- **(a) Scalar → vector**: $\hat b_x,\hat b_y$ 从标量变成 $\mathbb{R}^K$ 向量, $\mathrm{SE}^2$ 从标量变成 $K\times K$ 对角阵.
- **(b) 闭式 → 数值**: Step 4 的二次不等式有解析根; Step 6 的 $V(\beta_0)$ 是 $\beta_0$ 的二次函数, $AR$ 是有理函数, level set 没有闭式 (除 K=1), 必须 grid + uniroot.
- **(c) 单一假设检验 → over-id**: K=1 时 $AR(\beta_0)$ 只能逐点检验固定 $\beta_0$; K≥2 时 $\min_{\beta_0} AR(\beta_0)$ 本身是 Sargan-J $\sim \chi^2_{K-1}$, 能检验 K 个 mask 互不矛盾的联合假设.

---

## Step 5 · K ≥ 3 mask burden 堆向量

### 目标
把 Steps 1–3 输出的每个 (gene $g$, mask $k$, cell type $c$) 的标量四元组 $(\hat b_{x,k}, \mathrm{SE}_{x,k}, \hat b_{y,k}, \mathrm{SE}_{y,k})$ 沿 **mask 轴** 堆叠成长度-$K$ 向量, 给 Step 6 准备多 IV 输入.

### 数据结构
$$
\hat{\boldsymbol b}_x = \begin{pmatrix}\hat b_{x,1}\\ \vdots\\ \hat b_{x,K}\end{pmatrix},\quad
\hat{\boldsymbol b}_y = \begin{pmatrix}\hat b_{y,1}\\ \vdots\\ \hat b_{y,K}\end{pmatrix} \in \mathbb{R}^K,\quad
D_x = \mathrm{diag}(\mathrm{SE}_{x,1},\dots,\mathrm{SE}_{x,K}),\quad D_y = \mathrm{diag}(\mathrm{SE}_{y,1},\dots,\mathrm{SE}_{y,K}).
$$

加上三个 $K\times K$ 相关矩阵:
- $R_{xx}$ — 暴露 estimators 跨 mask 的相关 (默认 $I_K$).
- $R_{yy}$ — 结果 estimators 跨 mask 的相关 (默认 $I_K$).
- $R_{xy}$ — 跨样本 across-sample correlation block (默认 $\mathbf{0}_{K\times K}$, 见 Step 9).

### 默认值的合法性
- **$R_{xx} \approx R_{yy} \approx I_K$**: pLoF / missense:LC / regulatory 是三个 *不相交的* 变体集 (按 VEP annotation class 切分). 同一 gene 内, mask A 的变体跟 mask B 的变体之间的 LD 几乎为 0 (rare-variant LD blocks 比 mask 内的变体距离短得多). 因此 burden estimators 跨 mask 几乎不相关.
- **$R_{xy} = \mathbf{0}$**: 两样本 MR (eQTL 来自 OneK1K/TenK10K, GWAS 来自 Genebass) 默认无样本重叠. 一样本 / 重叠样本时改 $R_{xy}\neq 0$ (Step 9 详).

### 为什么 K ≥ 3 而非 K ≥ 2
HANDOVER §2 明文 commit: pLoF / missense:LC / regulatory, $K\ge 3$. 数学理由:
- $K=2$ ⇒ Sargan-J df = $K-1 = 1$. 1 个 over-id 约束, 一个 outlier mask 就能让 J 通过或拒绝, **诊断力很脆弱**.
- $K=3$ ⇒ df = 2, 真正能区分 "两 mask 同意 vs 三 mask 同意", J 检验有 power 定位异常 mask.
- $K\ge 3$ 是 **结构性要求**, 不是 power 优化 (`main.tex:198`).

### 代码定位 (`mrAR_multi.R`)
- 行 114–121: 函数签名 + 默认值 (`R_xx = diag(K)`, `R_yy = diag(K)`, `R_xy = matrix(0, K, K)`).
- 行 124–135: validate K 一致, R 矩阵维度 + 有限性.
- 行 153–154: 构造 $D_x = \mathrm{diag}(\mathrm{se\_x})$, $D_y = \mathrm{diag}(\mathrm{se\_y})$.

### 常见混淆点
- **mask 内的变体 ≠ Step 5 堆叠的 K**: mask 内的变体 $j\in\mathrm{mask}_k$ 是 *Step 1 burden 聚合的对象*; "K" 是 *mask 跨 annotation class 的个数* (pLoF, mis:LC, reg). 两者用同一字母 $j$ vs $k$ 在 `briefing_for_wei.md` §4.1 公式里区分.
- **K 不等于 HEIDI-rv 的 m**: HEIDI-rv (Step 10) 的 $m$ 是 mask 内 *变体数*, 用 leave-one-out 检验 burden 内异质性. Step 5 的 K 是 mask 数, 用 Sargan-J 检验 mask 间异质性. 这是 rvSMR 的 *两个独立 over-id 轴*.

---

## Step 6 · 多 IV AR 统计量 + 置信集

### 核心公式 (`main.tex:215-228`)
$$
m(\beta_0) \;=\; \hat{\boldsymbol b}_y \;-\; \beta_0\, \hat{\boldsymbol b}_x \;\in\; \mathbb{R}^K \tag{moment}
$$
$$
V(\beta_0) \;=\; D_y\, R_{yy}\, D_y \;+\; \beta_0^{\,2}\, D_x\, R_{xx}\, D_x \;-\; 2\,\beta_0\, D_y\, R_{xy}\, D_x \;\in\; \mathbb{R}^{K\times K} \tag{V}
$$
$$
\mathrm{AR}(\beta_0) \;=\; m(\beta_0)^\top\, V(\beta_0)^{-1}\, m(\beta_0) \;\xrightarrow{H_0:\,\beta=\beta_0}\; \chi^2_K \tag{AR}
$$
$$
\mathcal{C}_{1-\alpha} \;=\; \bigl\{\beta_0 : \mathrm{AR}(\beta_0)\le \chi^2_{K,\,1-\alpha}\bigr\}.
$$

### 解读每一项

**$m(\beta_0)$ — reduced-form 残差.** 在候选因果效应 $\beta_0$ 下, $\hat{\boldsymbol b}_y$ 减去 IV 暗示的 indirect path. 真值 $\beta_0=\beta$ 处 $\mathbb{E}[m]=0$ (依 exclusion). 注意 *从未 invert $\hat{\boldsymbol b}_x$* — 这就是 weak-IV 鲁棒性的根源.

**$V(\beta_0)$ 三项**:
1. $D_y R_{yy} D_y$ — 结果 SE 贡献, 跟 $\beta_0$ 无关.
2. $\beta_0^2 D_x R_{xx} D_x$ — 暴露 SE 贡献, 按 $\beta_0^2$ 缩放 (因为 $\beta_0 \hat{\boldsymbol b}_x$ 在 $m$ 里乘了 $\beta_0$).
3. $-2\beta_0 D_y R_{xy} D_x$ — 跨样本 cross-term (Step 9 详). 顺序是 $D_y R_{xy} D_x$, 不是 $D_x R_{xy} D_y$ (见 pitfall 3).

**为什么 $V$ 不能简化成 $V_{yy}$ 而已**: $m = \hat{\boldsymbol b}_y - \beta_0 \hat{\boldsymbol b}_x$ 同时包含 $\hat{\boldsymbol b}_y$ 和 $\beta_0 \hat{\boldsymbol b}_x$ 的随机性. 即使在 $H_0$ 下 $\beta_0$ 是固定常数, $\hat{\boldsymbol b}_x$ 仍是随机的; 把它的 SE 漏掉相当于假设 $\hat{\boldsymbol b}_x$ 无误差, 直接退化成 OLS reduced-form 检验, 不再 pivotal 在弱 IV 下.

### R 代码对应 (`mrAR_multi.R`)
- 行 153–157: 预计算 $V_{yy} = D_y R_{yy} D_y$, $V_{xx} = D_x R_{xx} D_x$, $V_{xy} = D_y R_{xy} D_x$. 这三个矩阵跟 $\beta_0$ 无关, 提到 closure 外避免每次评估重算.
- 行 159: `c_crit <- stats::qchisq(1 - alpha, df = K)` — **df = K, 不是 K−1** (Sargan-J 才是 K−1, 见 pitfall 1).
- 行 162–171: `ar_fun` closure, 每次输入 $\beta_0$ 标量, 返回 $AR(\beta_0)$:
  - 行 163: `m <- b_y - b0 * b_x`
  - 行 164: `Vb <- V_yy + (b0 * b0) * V_xx - (2 * b0) * V_xy`
  - 行 166–168: `solve(Vb, m)` 再 `crossprod(m, sol)` 得 $m^\top V^{-1} m$. 失败 (singular $V$) 返 NA, 后续算法在 grid 上把 NA 视为 +∞.

### 关键设计选择 — 为什么要数值反演
$V(\beta_0)$ 是 $\beta_0$ 的 *矩阵值二次函数*, 它的逆是有理函数; $m^\top V^{-1} m$ 在 $\beta_0$ 维上是 *有理函数*. 在 K=1 时分子分母都是关于 $\beta_0$ 的二次, 不等式 $AR \le c$ 化简为标量二次不等式有闭式根. K≥2 时分子是 $\beta_0$ 的二次型的内积, 分母是 $K\times K$ 矩阵的行列式 (degree $\le 2K$), 没有解析根. 所以必须数值 (Step 7).

### 弱 IV 鲁棒性的来源
$V(\beta_0)$ 的下界由 $V_{yy}$ 给出 (正定, 跟 $\hat{\boldsymbol b}_x$ 大小无关). 因此 $AR(\beta_0)$ 在 first-stage $F\to 0$ 时仍 well-defined; coverage 不依赖 IV 强度 (`main.tex:234` 强调 "uniform in IV strength" 的关键就是 $\hat{\boldsymbol b}_x$ 从未被 invert). 这是 §4.2 of `briefing_for_wei.md` 在 K=1 时讨论的 pivotal 性质, K≥2 直接继承.

---

## Step 7 · 数值反演: grid + uniroot

### 算法逻辑 (5 步)

1. **初始 envelope**: per-IV Wald 比 $\{\hat b_{y,k}/\hat b_{x,k}\}_{k=1}^K$ 的 min/max, 加 padding `grid_pad_mult * max|Wald_k|` (默认乘数 3). $\hat b_{x,k}$ 太接近 0 的 Wald 比 NA 过滤掉.
2. **Grid 评估**: 在 envelope 内 `n_grid = 4000` 等距点上算 $AR(\beta_0) - c_{\rm crit}$; $V(\beta_0)$ singular 处 mark 为 +∞ (≥ c_crit, 算作拒绝).
3. **Envelope 延伸**: 若 grid 两端点的 $AR$ 仍 $\le c_{\rm crit}$, 说明 level set 漏出 envelope, half-width 翻倍重 grid, 最多 `grid_extend_max = 3` 次. 之后仍漏就分类为 `whole_line`.
4. **Sign-change 检测**: 找 $AR - c_{\rm crit}$ 在相邻 grid 点之间换号的 index.
5. **uniroot 精修**: 每个 sign-change 区间调 `stats::uniroot()`, `tol = 1e-8`, 最多 200 iter. 失败返 NA, 后续过滤掉非有限 root, 排序去重.
6. **分类**: 用精修 roots 把实轴切成段, 段内取中点 (端点 ±∞ 取 ±1 偏移) 评 $AR$ 判段是否属于置信集. 分类:
   - `bounded_interval`: 1 段 accept, 两端有限 (或一端 ±∞ 且 grid 已穷尽延伸).
   - `disconnected_union`: ≥ 2 段 accept.
   - `whole_line`: 0 roots 且中心点 accept.
   - `empty`: 没有 accept 段.

### CI shape 跨 K 的传承
K=1 闭式四分类 (`briefing_for_wei.md` §4.3 table) — Strong+identified→bounded, Strong+empty (≤α prob), Weak+disconnected_union, Weak+whole_line — 在 K≥2 由 grid 几何 *自然* 出现. 没有 K 维专属的新 CI shape; 只是判定路径从"求二次不等式的根"换成"扫 level set".

### 代码定位 (`mrAR_multi.R`)
- 行 175–189: 初始 envelope. `wald_per <- ifelse(abs(b_x) > eps, b_y/b_x, NA)`, 过滤有限值, 设 `pad = grid_pad_mult * max(|wald_per|, 1)`, `grid_lo/hi` 为 padded 极值.
- 行 191–220: grid build + extension `repeat` loop. 每轮算 `ar_grid`, singular 处 → +Inf, 检查两端是否 leak, leak 且未到 max 就翻倍 span.
- 行 222–225: `sc_idx <- which(diffs[-length(diffs)] * diffs[-1L] < 0)` — 经典 sign-change.
- 行 228–242: `refine()` 用 `stats::uniroot()` + `tryCatch`, 应用到每个 `sc_idx`.
- 行 243–244: `roots <- roots[is.finite(roots)]`, `sort(unique(roots))` — 过滤 uniroot 失败.
- 行 246–319: CI 分类. 注意 行 271 端点序列 `c(-Inf, roots, Inf)`, 中点 logic 在 行 276–284 (handle ±∞ 端).

### 数值陷阱
- **`grid_extend_max` 默认 3**: 防止 whole_line 时无限延伸. 3 次翻倍 = 原 envelope 的 $2^3 = 8$ 倍宽度, 通常够.
- **$V(\beta_0)$ singular**: 在某些 $\beta_0$ (e.g. 大 $|\beta_0|$ + 重叠样本) 可让 $V$ indefinite. `tryCatch(solve(Vb, m))` 失败 → NA → 算法 mark +Inf (行 207–208), 不破坏 sign-change 逻辑.
- **`tol = 1e-8`**: uniroot 边界点精度 ≤ 1e-8, publication-grade.
- **`n_grid = 4000`**: 太疏可能漏掉两个挨得很近的 root, 把 `disconnected_union` 误判成 `bounded_interval` (pitfall in `main.tex:274`).

### 为什么不用 optim
`optim` (或 `optimize`) 给的是 **argmin**, 那是 Step 8 的事. Step 7 要的是 **level set 边界** ($\{\beta_0 : AR(\beta_0) = c_{\rm crit}\}$). 数值反演 level set 的标准做法就是 sign-change 扫 + uniroot 精修; `optim` 不直接给边界.

---

## Step 8 · Sargan-J 过识别检验

### 核心公式 (`main.tex:284-291`)
$$
\hat\beta_{\mathrm{AR}} \;=\; \arg\min_{\beta_0}\, \mathrm{AR}(\beta_0)
$$
$$
J \;=\; \mathrm{AR}\bigl(\hat\beta_{\mathrm{AR}}\bigr) \;\xrightarrow{\text{joint IV validity}}\; \chi^2_{K-1}
$$
$$
\text{p-value} \;=\; 1 - F_{\chi^2_{K-1}}(J).
$$

拒绝同质性 (即声明 horizontal pleiotropy across $K$ masks) 当且仅当 $J > \chi^2_{K-1, 1-\alpha}$.

### 为什么 df = K−1 而非 K
- 在 *固定* $\beta_0$ 下 $AR(\beta_0)\sim \chi^2_K$ (这是 Step 6 的 reference).
- 取 $\min_{\beta_0}$ 等于估了 1 个参数 ($\hat\beta_{\rm AR}$). 类比 OLS: 残差平方和 $\sum(y_i - \hat y_i)^2 / \sigma^2 \sim \chi^2_{n-p}$ — 拟合 $p$ 参数后 df 减 $p$.
- 这里 $p=1$ (单一 causal effect), 所以 J ~ $\chi^2_{K-1}$.
- 这是 Hansen 1982 GMM J-statistic 在两样本 summary MR 设定下的实例 (`HANDOVER §2` 引 Patel–Lane–Burgess 2024).

### R 代码 (`mrAR_multi.R`)
- 行 322–326: `stats::optimize(ar_fun, lower=grid_lo, upper=grid_hi, tol=1e-8)`. 在 grid envelope 内做 1D 凹面搜索. `tryCatch` 包错防 NA-边界.
- 行 327–331: fallback: `which.min(ar_grid_finite)`, 用 grid 上的 argmin (粗但稳).
- 行 332–335: 成功时取 `opt$minimum` (β̂) 和 `opt$objective` (J).
- 行 337–341: `J_pvalue <- pchisq(J_stat, df = K-1, lower.tail = FALSE)`, K==1 时返 `NA_real_` (没有 over-id).

### 诊断意义
- **J 不显著 (大 p)**: K 个 mask 在 $\hat\beta$ 处 reduced-form 残差互不矛盾 ⇒ 支持联合 IV 有效 ⇒ 报 $\hat\beta_{\rm AR}$ 作为点估计 + AR CI 作为区间.
- **J 显著 (小 p)**: 至少 1 个 mask 违反 exclusion / monotonicity / linearity. **但 J 是 joint test, 不告诉你是哪个 mask**. 定位需要:
  - HEIDI-rv (Step 10): mask 内 leave-one-out.
  - Annotation-class concordance (Step 11): mask 间 Cochran-Q.
  - Cell-type concordance: cell type 间一致性.
- **Coherent pleiotropy 警告** (`briefing_for_wei.md` §4.4): 如果三 mask 都被同一上游 confounder 污染, J 仍可能不显著 (因为 mask 间一致地被偏). 三轴 over-id 都无法捕捉, 需要 Cinelli–Hazlett 2025 sensitivity (Step 12).

---

## Step 9 · 样本重叠修正

### K=1 形式 (`mrAR.R:116`, `main.tex:312-313`)
$$
\mathrm{Var}\!\bigl(\hat\beta_y - \beta_0 \hat\beta_x\bigr) \;=\; \mathrm{SE}_y^2 \;+\; \beta_0^2 \mathrm{SE}_x^2 \;-\; 2\beta_0\,\rho\, \mathrm{SE}_x \mathrm{SE}_y.
$$
进入二次系数 $B = -2(\hat b_x \hat b_y - c\, \rho\, \mathrm{SE}_x \mathrm{SE}_y)$ (`mrAR.R:116`).

### K≥2 形式 (`main.tex:316-319`)
$V(\beta_0)$ 含矩阵 cross-term
$$
-\,2\,\beta_0\, D_y\, R_{xy}\, D_x.
$$
进入 `mrAR_multi.R:157` 的 `V_xy <- Dy %*% R_xy %*% Dx`, 再在 closure 行 164 `Vb <- V_yy + b0² V_xx − 2 b0 V_xy` 合成.

### 核心 geometric property — sign-flip with β₀
Cross-term 系数是 $-2\beta_0$. 当 $\beta_0$ 跨过 0, cross-term 反号 ⇒
- 在 $\beta_0 = 0$ 处 sample overlap **不影响 $AR(0)$**.
- 远离 0 时, $\beta_0$ 与 $\rho$ 同号时 cross-term 减小 $V$ ⇒ $AR$ 增大 ⇒ CI 在那侧 *变窄*.
- 反号时 $V$ 变大 ⇒ $AR$ 变小 ⇒ CI 在那侧 *变宽*.
- 结果: $\rho\neq 0$ 让 CI **关于 0 不对称** (`main.tex:325`).

**Type-I 含义**: 一样本分析里默认 $\rho=0$ 是 silent type-I error inflation, 尤其当 $\beta_0$ 和 $\hat\beta$ 同号时 ($V$ 被低估 ⇒ AR 被高估 ⇒ 拒绝更激进).

### 典型场景
- **两样本不重叠**: $\rho = 0$, $R_{xy} = \mathbf{0}$ (rvSMR v1 的默认 — eQTL OneK1K, GWAS Genebass).
- **一样本**: $\rho \approx 1 - \mathrm{noise\ overlap}$, 两 estimator 共享个体.
- **部分重叠**: $\rho \approx N_{\rm overlap} / \sqrt{N_x N_y}$ × (phenotypic correlation).
- 真值估计: LD-Score regression cross-trait intercept, 但要小心 — LDSC intercept conflates overlap + cross-trait genetic covariance (`main.tex:333`).

### R 代码默认与 plug-in
- `mrAR.R:88` 默认 `cor_xy = 0`.
- `mrAR_multi.R:117` 默认 `R_xy = matrix(0, length(b_x), length(b_y))` — 即 $K\times K$ 零矩阵 (注意这里用 `length(b_y)` 维度, 但 K=length(b_x)=length(b_y) 已 validate, OK).
- 非零时调用方传 `R_xy = ρ * J_K` (常数 block) 或 mask-specific 矩阵.

### 代码定位 (`mrAR_multi.R`)
- 行 157: `V_xy <- Dy %*% R_xy %*% Dx` — **顺序是 $D_y R_{xy} D_x$, 不是 $D_x R_{xy} D_y$**. $(i,j)$ entry 是 $\mathrm{SE}_{y,i} R_{xy,ij} \mathrm{SE}_{x,j}$ (`main.tex:241`).
- 行 164: `Vb <- V_yy + (b0 * b0) * V_xx - (2 * b0) * V_xy` — closure 里 V_xy 进入.

---

## 5 步之间的依赖图

```
                          Step 5 (stack)
                            │
                            ▼
              ┌──── Step 9 (R_xy default = 0) ────┐
              │                                   │
              │  ┌── Step 6 (define AR statistic) ┤
              │  │                                │
              ▼  ▼                                ▼
          Step 7 (grid + uniroot)            Step 8 (J statistic)
              │                                    │
              ▼                                    ▼
          CI shape (4 types)                   pleiotropy decision
                                            (homogeneous / reject)
```

Step 9 不是顺序后置, 而是 **Step 6 的参数注入**: $R_{xy}$ 默认 0 时 Step 6 仍能算; 非零时仅改 $V(\beta_0)$ 的 cross-term, 之后 Step 7/8 完全继承.

---

## 数据流: 从 Step 3 输出到 Step 8 终点 (具体例子)

设 PCSK9 gene 在 1 个 cell type, K=3 masks. Step 3 的 per-mask 输出 (借用 `mrAR.R:81-83` 的 K=1 示例数值扩成 K=3):

```
  pLoF    :  b_x = 0.18,  SE_x = 0.04   b_y = -0.11,  SE_y = 0.025
  mis:LC  :  b_x = 0.15,  SE_x = 0.05   b_y = -0.10,  SE_y = 0.030
  reg     :  b_x = 0.10,  SE_x = 0.07   b_y = -0.07,  SE_y = 0.035
```

### Step 5 输出
$$
\hat{\boldsymbol b}_x = (0.18, 0.15, 0.10)^\top, \quad \hat{\boldsymbol b}_y = (-0.11, -0.10, -0.07)^\top,
$$
$$
D_x = \mathrm{diag}(0.04, 0.05, 0.07), \quad D_y = \mathrm{diag}(0.025, 0.030, 0.035), \quad R_{xx} = R_{yy} = I_3, \quad R_{xy} = \mathbf{0}_3.
$$

每 mask Wald 比: $-0.11/0.18 \approx -0.611$, $-0.10/0.15 \approx -0.667$, $-0.07/0.10 = -0.700$. 三 mask 大致同意 $\beta \approx -0.65$ — 这是 IV 联合有效的 *预期* 签名.

### Step 6 评估 — 取 $\beta_0 = -0.6$
$$
m(-0.6) = \hat{\boldsymbol b}_y - (-0.6)\hat{\boldsymbol b}_x = (-0.11 + 0.108,\; -0.10 + 0.090,\; -0.07 + 0.060) = (-0.002,\; -0.010,\; -0.010).
$$
$$
V(-0.6) = D_y^2 + 0.36\, D_x^2 = \mathrm{diag}(0.025^2 + 0.36\cdot 0.04^2,\; 0.030^2 + 0.36\cdot 0.05^2,\; 0.035^2 + 0.36\cdot 0.07^2)
$$
$$
= \mathrm{diag}(1.201\times 10^{-3},\; 1.800\times 10^{-3},\; 2.989\times 10^{-3}).
$$
$$
AR(-0.6) = \sum_k m_k^2 / V_{kk} = \frac{(0.002)^2}{1.201\!\times\!10^{-3}} + \frac{(0.010)^2}{1.800\!\times\!10^{-3}} + \frac{(0.010)^2}{2.989\!\times\!10^{-3}}
$$
$$
\approx 0.0033 + 0.0556 + 0.0335 \approx 0.092.
$$
$c_{\rm crit} = \chi^2_{3, 0.95} \approx 7.815$. $AR(-0.6) = 0.092 \ll 7.815$ ⇒ $\beta_0 = -0.6$ **属于 CI**.

### Step 7 反演
扫 grid: per-IV Wald 比的 envelope 大致是 $[-0.70, -0.61]$, padding 后比如 $[-2.8, +1.4]$. 在 envelope 上 $AR$ 应该在 $\beta_0\approx -0.65$ 附近 minimum (约 0 量级), 两侧 $AR$ 单调上升, 两 root 大约对称地框出 bounded interval (例如 $\approx [-0.92, -0.41]$ 这种宽度, 实际值取决于精确 grid 评估).

### Step 8 J 检验
$\hat\beta_{\rm AR} = \arg\min AR \approx -0.65$, $J = AR(\hat\beta_{\rm AR})$ 很小 (因为三 Wald 比几乎一致, residual $m(\hat\beta) \approx 0$). $J \sim \chi^2_{K-1} = \chi^2_2$, $\chi^2_{2, 0.95} = 5.99$. $J \ll 5.99$ ⇒ p-value 接近 1 ⇒ **不拒绝 mask 间 homogeneity** ⇒ 报 $\hat\beta \approx -0.65$ + CI 区间作为最终结果.

(若数据里三个 Wald 比明显分歧 — e.g. pLoF 给 -1.5, reg 给 -0.3 — $\min AR$ 不会接近 0, $J$ 可能 > 5.99, 触发 pleiotropy alarm.)

---

## 五步在论文里要分别说什么 (writing guidance)

| Step | 写法建议 |
|---|---|
| **5** | 一段散文 + 一行公式 (stacked vectors + $D_x, D_y$) + 一句 "$K\ge 3$ commitment" 引 `HANDOVER §2`. 强调 $R_{xx}=R_{yy}=I_K$ 的合法性来自 mask 不相交. |
| **6** | 两段 + 三行公式 ($m$, $V$, $AR$) + 一段 "为什么 $V$ 不简化为 $V_{yy}$" + 一句 pivotal/weak-IV-robust. |
| **7** | 一段 + algorithm float (5-6 行 pseudocode: envelope, grid, sign-change, uniroot, classify) + 一句 "为什么 grid+uniroot 不是 optim". |
| **8** | 一段 + 一行公式 + 一段 "df = $K-1$ 的直觉" (OLS 类比) + 一段 J 不告诉你是哪个 mask. |
| **9** | 一段 + 强调 sign-flip with $\beta_0$ + 一句"两样本默认 0". |

---

## 容易踩的坑 (避雷)

1. **写成 $AR(\beta_0)\sim \chi^2_{K-1}$**: 错. 在 $H_0:\beta=\beta_0$ *固定* 下是 $\chi^2_K$; 只有 $J = \min AR = AR(\hat\beta_{\rm AR})$ 才是 $\chi^2_{K-1}$. `main.tex:242` 显式警告.
2. **写 $V(\beta_0) = D_y R_{yy} D_y + \beta_0^2 D_x R_{xx} D_x$ 丢掉 cross-term**: 当 $R_{xy}\neq 0$ 时错; 即使默认 0 也要保留在公式里, 否则 sample-overlap 通道消失.
3. **写 $V_{xy} = D_x R_{xy} D_y$ 顺序反**: 应是 $D_y R_{xy} D_x$. 跟 `mrAR_multi.R:157` 和 `main.tex:241` 一致. $(i,j)$ entry 是 $\mathrm{SE}_{y,i} \cdot R_{xy,ij} \cdot \mathrm{SE}_{x,j}$, 不是反过来.
4. **把 Sargan-J 跟 Hansen-J 写成不同检验**: 在 over-identified 单参数 IV 下两者等价 (Hansen 1982 把 Sargan 推广到 GMM). `main.tex:291` 用 "Sargan / GMM over-identification" 双名.
5. **用 Stock-Yogo $F > 10$ filter 然后说 AR 鲁棒**: 自相矛盾. AR 的全部意义就是 *不需要* filter (`briefing_for_wei.md` §4.2). `mrAR.R:26-27` 显式写 "Do not combine an AR CI with $F > 10$".
6. **CCT across the K masks within a gene**: `HANDOVER §2` 明令禁止 — 会把 $K$ 个 p-value collapse 成一个, 消掉 over-id 信息. CCT 只用于 across-gene / across-cell-type / across-MAF-weight ensembles.
7. **报 $\hat\beta_{\rm AR}$ 不报 CI**: 在弱 IV 下 argmin 可能远离 accepted set 主体. CI 的 type (bounded / disconnected / whole_line / empty) 是 *第一类* 输出, $\hat\beta$ 是 *附加* 信息.

---

## 跟 K=1 closed-form (Step 4) 的差异表

| 维度 | Step 4 (K=1, `mrAR.R`) | Steps 5–8 (K≥2, `mrAR_multi.R`) |
|---|---|---|
| Statistic | scalar quadratic $AR(\beta_0) = (b_y-\beta_0 b_x)^2 / \mathrm{den}$ | matrix quadratic form $m^\top V^{-1} m$ |
| Reference null | $\chi^2_1$ | $\chi^2_K$ |
| Solver | 闭式 quadratic roots ($A\beta_0^2 + B\beta_0 + C \le 0$) | grid + uniroot, level-set inversion |
| Sign-of-$A$ criterion | $A = b_x^2 - c\,\mathrm{SE}_x^2$, $\mathrm{sign}(A) = \mathrm{sign}(F-c)$ | implicit in CI shape from grid geometry |
| Over-id test | none ($\equiv$ Fieller 1954) | Sargan-$J \sim \chi^2_{K-1}$ at argmin |
| CI shapes | 4 (bounded / disconnected / whole_line / empty) | same 4 (numerically detected) |
| Sample-overlap entry | $\rho$ scalar in $B$ coefficient (`mrAR.R:116`) | $R_{xy}$ $K\times K$ in $V_{xy}$ (`mrAR_multi.R:157`) |
| Computational cost | $O(1)$ | $O(n_{\rm grid} \times K^3)$ per CI |
| R impl | ~200 lines (`mrAR.R`) | ~360 lines (`mrAR_multi.R`) |
| Test coverage | 27 testthat assertions | 34 + K=1 cross-check (~1e-10 agreement) |

K=1 cross-check: `mrAR_multi(K=1)` 必须跟 `mrAR()` 在所有四个 CI shape 上数值同步到 ~1e-10. 这是 `HANDOVER §4` 报告的 regression test, 也是 `mrAR_multi.R` 头部 docstring 行 7 "the closed-form K=1 sibling" 的来源 — multi 版本通过 K=1 cross-check 自验证.
