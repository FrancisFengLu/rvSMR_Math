# Pedagogy Review — Sub-worker 2C

Reads `content_draft.md` as a SMART BUT UNTRAINED bioinformatician would. Flags (a) intuition-before-formula gaps, (b) method introductions without paper attribution, (c) awkward Chinese, (d) missing SVG / diagrams that would help, (e) missing "why not alternatives".

---

## Step-by-step pedagogical review

### Step 0 (assemble inputs)

- (a) Intuition: ✅ Good — opens with "what we're trying to know" and works up to notation.
- (b) Paper attribution: ✅ All three (Pierce-Burgess 2013, BBT 2013, Burgess 2016) tied to method.
- (c) Chinese: mostly fine. One awkward phrase: "把每一项数字的'身份证'固定下来" — keep, it's vivid.
- (d) SVG needed: **YES — a small DAG showing (g, k, c) tuple data flow from SAIGE-QTL + Genebass → rvMR input table.** Add.
- (e) Why not alternatives: ✅ Has "no individual-level IV, no single-variant MR".

### Step 1 (burden construction)

- (a) Intuition: ✅ Clear — "many rare variants combined as one signal".
- (b) Paper attribution: ✅ Madsen-Browning + STAAR + SKAT.
- (c) Chinese: clean.
- (d) SVG needed: **YES — a small illustration showing 5 rare variants in a mask, each with MAF + weight, summing to one burden Z_k.** This is the most-confusing concept for a newcomer.
- (e) Why not alternatives: ✅ CAST, CMC, SKAT, ACAT, STAAR-O all named.

### Step 2 (acquire summary stats)

- (a) Intuition: ✅ "Step 1 is HOW; Step 2 is WHERE TO GET".
- (b) Paper attribution: ✅ SAIGE-QTL, SAIGE-GENE+, Genebass, RGC-ME.
- (c) Chinese: clean.
- (d) SVG: optional — could show two-arm flow (TenK10K → b_x; UKB exomes → b_y). Skip if tight on space.
- (e) Alternatives: ✅ RareEffect BLUP, STAARpipeline output, RVAT all named with reason.

### Step 3 (Wald ratio + F)

- (a) Intuition: ✅ Builds Wald ratio from IV first principles before formula.
- (b) Paper attribution: ✅ Burgess-Labrecque, Didelez-Sheehan, Lee-McCrary-Moreira-Porter.
- (c) Chinese: clean.
- (d) SVG needed: **YES — small DAG: Z → X → Y, with side arrow Z → Y blocked by exclusion**. This is the core IV identification picture; critical to show explicitly.
- (e) Alternatives: ✅ dIVW, MR-RAPS, Stock-Yogo named.

### Step 4 (K=1 AR — the HARDEST step)

- (a) Intuition: ✅✅ Multi-paragraph build (4 paragraphs) — this is the load-bearing pedagogy. Reasonable.
- (b) Paper attribution: ✅ Anderson-Rubin 1949, Fieller 1954, Wang-Kang 2022, Patel-Lane-Burgess 2024.
- (c) Chinese: clean.
- (d) SVG needed: **CRITICAL — must include a parabola / phase-diagram showing the 4 CI shapes** (bounded interval, empty, disconnected union, whole line) as a function of (sign A, sign Δ). This is the visual key to Step 4.
- (e) Alternatives: ✅ Wald, Fieller, dIVW, MR-RAPS, Stock-Yogo all named.
- **KEY PEDAGOGICAL MOVE (per user spec)**: The 4-paragraph build is the move — paragraph 1 = "Wald fails", paragraph 2 = "AR flips the question", paragraph 3 = "why pivotal", paragraph 4 = "4 shapes". Then a parabola SVG cements the geometry. The reader can SEE why disconnected/whole_line are honest answers.

### Step 5 (stack K masks)

- (a) Intuition: ✅ Sets up over-id motivation.
- (b) Paper attribution: ✅ Sargan 1958 + Patel-Lane-Burgess 2024.
- (c) Chinese: clean.
- (d) SVG: optional — could show 3 masks (pLoF, mis:LC, reg) stacked as K=3 IV vector. Skip if tight.
- (e) Alternatives: ✅ K=2, K=1 stratified by cell type, ACAT collapse all named.

### Step 6 (multi-IV AR)

- (a) Intuition: ✅ Generalizes Step 4. "Never inverts b_x" — core property emphasized.
- (b) Paper attribution: ✅ A-R + Wang-Kang + Patel-Lane-Burgess.
- (c) Chinese: clean.
- (d) SVG: optional — could show V(β₀) decomposition into 3 terms.
- (e) Alternatives: ✅ Meta-analyze K=1 ARs, IVW, 2SLS, MR-Egger all named.

### Step 7 (grid + uniroot)

- (a) Intuition: ✅ "No closed form for K≥2, so we go numerical".
- (b) Paper attribution: ⚠️ Light — only Brent 1973 implicit for uniroot. Says "engineering implementation, not new paper". OK.
- (c) Chinese: clean.
- (d) SVG needed: **YES — a small inversion picture showing AR(β₀) function with horizontal line at c_crit, sign changes marked, roots boxed, intervals colored accept/reject.** This is the algorithm picture.
- (e) Alternatives: ✅ Closed-form, profile likelihood, bootstrap all named.

### Step 8 (Sargan-J)

- (a) Intuition: ✅✅ Build is clear — "if K mask agree, min AR is small; df = K-1 because we fit 1 parameter".
- (b) Paper attribution: ✅ Sargan 1958 + Hansen 1982 + Patel-Lane-Burgess 2024.
- (c) Chinese: clean.
- (d) SVG: optional — could illustrate K=3 case with 3 AR curves crossing.
- (e) Alternatives: ✅ Cochran Q on Wald, MR-Egger, MR-PRESSO all named.

### Step 9 (sample overlap)

- (a) Intuition: ✅ Sign-flip insight explicit.
- (b) Paper attribution: ✅ Burgess-Davies-Thompson 2016 + A-R 1949 + Wang-Kang 2022.
- (c) Chinese: clean.
- (d) SVG: optional — could show cross-term parabola showing asymmetry.
- (e) Alternatives: ✅ LDSC intercept, disjoint subsamples, ignore named.

### Step 10 (HEIDI-rv)

- (a) Intuition: ✅✅ Long build — "Step 8 = mask-between, Step 10 = mask-within"; explains contrast matrix.
- (b) Paper attribution: ✅ Zhu 2016 + Davies 1980 + Kuonen 1999.
- (c) Chinese: clean.
- (d) SVG: optional — could show m per-variant Wald ratios as scatter points + their mean line, with one outlier flagged.
- (e) Alternatives: ✅ LD-based HEIDI, MR-PRESSO, single-variant LOO named.
- **CRITICAL trap callout**: V_δ^+ V_δ vs V_δ headline trap is explicit. ✅

### Step 11 (annotation Q)

- (a) Intuition: ✅ "raw is not comparable" → "pQTL anchor normalizes" → "Q on normalized".
- (b) Paper attribution: ✅ Cochran 1954 + Dhindsa 2023 + Ferkingstad 2021.
- (c) Chinese: clean.
- (d) SVG: optional — could show 3 class slopes before & after normalization.
- (e) Alternatives: ✅ Skip normalization, common-variant anchor, Han-Eskin all named as "do not".

### Step 12 (cell-type Q)

- (a) Intuition: ✅ "third over-id axis".
- (b) Paper attribution: ✅ Cuomo 2025 + Ray 2025 + Cochran 1954.
- (c) Chinese: clean.
- (d) SVG: optional — bar chart of 28 cell-type Wald ratios with bar heights.
- (e) Alternatives: ✅ OneK1K, bulk GTEx, pseudobulk named.

### Step 13 (sensitivity)

- (a) Intuition: ✅ "Coherent pleiotropy bypasses all over-id; sensitivity quantifies the bypass strength".
- (b) Paper attribution: ✅ Cinelli-Hazlett dual + Swanson-VanderWeele + VanderWeele-Ding.
- (c) Chinese: clean.
- (d) SVG: optional — could plot RV as fn of t.
- (e) Alternatives: ✅ MR-Egger, bidirectional MR, partial-id bounds named.

### Step 14 (decision rule)

- (a) Intuition: ✅ "conjunctive (all axes must pass), conservative (default inconclusive under coherent pleiotropy)".
- (b) Paper attribution: ✅ Wang-Tchetgen 2018 + Bowden 2015 + Cinelli-Hazlett 2025.
- (c) Chinese: clean.
- (d) SVG needed: **YES — final decision flowchart D1 → D2 → ... → D6 → label**. This is the closing summary visual.
- (e) Alternatives: ✅ Single-test, voting, Bayesian average named.

---

## Cross-cutting suggestions

1. **Top-of-doc TOC + pipeline overview SVG**: One global SVG at top showing all 14 steps as boxes with arrows. Anchor for navigation + mental model.

2. **Color coding consistency**: Use distinct accent colors per step group:
   - Steps 0-2 (input prep): teal `#2dd4bf`
   - Steps 3-4 (K=1 inference): blue `#58a6ff`
   - Steps 5-9 (multi-IV core): purple `#a78bfa`
   - Steps 10-12 (over-id axes): amber `#d29922`
   - Steps 13-14 (sensitivity + decision): green `#3fb950`

3. **Citation cards**: Each "用了哪篇 paper" subsection should be a distinct card with amber `#d29922` left border. User explicitly requested this style.

4. **MathJax equations**: All formulas via MathJax v3 (per user spec).

5. **Mobile**: Cards stack vertically; sticky TOC collapses to hamburger.

6. **Pitfall callouts**: Use a small "陷阱" badge with red `#f47067` background to make traps visually distinct.

---

## Final recommendation: 6 SVGs to inline

1. **TOC pipeline overview** (top of doc): all 14 steps as boxes.
2. **Step 0 — (g, k, c) tuple flow**: SAIGE-QTL + Genebass → rvMR.
3. **Step 1 — burden aggregation**: 5 variants → 1 burden.
4. **Step 3 — IV DAG**: Z → X → Y with exclusion.
5. **Step 4 — 4 CI shapes**: 2x2 grid (sign A) × (sign Δ).
6. **Step 7 — AR inversion**: AR function with c_crit line + sign-change roots.
7. **Step 14 — Decision flowchart**: D1 → D6 → label.

These 7 SVGs cover the key conceptual leaps. Other steps can rely on text + math.

---

## Verdict

Content draft is solid. Add SVGs (especially Steps 4, 7, 14), apply color coding, render in HTML with citation cards. No content rewrites needed.

End of pedagogy review.
