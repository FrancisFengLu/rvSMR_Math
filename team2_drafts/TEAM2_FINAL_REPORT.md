# Team 2 Final Report

**Deliverable:** `/home/francisfenglu4/projects/rvSMR_Math/algorithm_paper_walkthrough.html`

## Headline numbers

- HTML line count: **1900 lines**
- Sections: **14 step sections + intro + pipeline overview + summary + references** (17 logical sections, 14 step-sections per spec)
- Papers cited (unique, in references list): **39**
- Citation cards (with amber border per spec): **19**
- Inline SVGs: **6** (pipeline overview, burden aggregation, IV DAG, 4 CI shapes, AR inversion, decision flowchart)

## Sub-worker outputs

All three sub-worker drafts saved to `/home/francisfenglu4/projects/rvSMR_Math/team2_drafts/`:
- `content_draft.md` — sub-worker 2A
- `citation_audit.md` — sub-worker 2B (16 trap items verified; no errors in draft)
- `pedagogy_review.md` — sub-worker 2C (recommended 7 SVGs; 6 inlined)

## Key pedagogical move for Step 4 (the conceptually-hardest step)

**4-paragraph build before the formula**:
1. "Wald fails" — establish problem from $\hat b_x^2 / \hat b_x^4$ in delta-SE denominator.
2. "AR flips the question" — test candidate $\beta_0$ rather than estimate $\beta$.
3. "Why pivotal" — denominator bounded below by $\mathrm{SE}_y^2$; never invert $\hat b_x$.
4. "Four shapes" — sign(A) × sign(Δ) → 2×2 geometric typology.

Then a 4-panel 2×2 SVG showing parabolas + critical line + accepted regions for each shape. The reader can *see* why disconnected and whole-line are honest answers, not pathologies. This is the move that the previous explainer "jumped over" — here it has 4 paragraphs of intuition before the first formula appears.

## Citation issues flagged or corrected

Two minor categories worth recording:

1. **Step 10 math erratum (CAUGHT BY TEAM 1, integrated by Team 2)**: literal main.tex eq (32-34) pairs Mahalanobis form $T = \delta^\top V_\delta^+ \delta$ with Davies weights = nonzero eigenvalues of $V_\delta$. These are inconsistent — Mahalanobis form follows plain $\chi^2_{\mathrm{rank}(V_\delta)}$, not generalized $\chi^2$. Team 1 simulated and verified literal pairing → 99% Type-I error at m=4. rvMR implements corrected $T = \delta^\top \delta$ + Davies weights = eig($V_\delta$). HTML Step 10 documents this erratum and uses the corrected math.

2. **Comparator references not in main.tex bibliography**: Bowden 2015 MR-Egger, Verbanck 2018 MR-PRESSO, Liu 2019 ACAT, Hemani 2017 bidirectional MR, Ye 2021 dIVW, Zhao 2020 MR-RAPS. All cited as comparator background in the "为什么不用替代方案" subsections. Where each is needed for academic accuracy (Bowden 2015 specifically called out by user spec for Step 14), it's now in the references list.

## Cross-team coordination outcome

Team 1 implemented all 5 stubs by 2026-06-08 MS1 (160 tests PASS / 0 FAIL). HTML's code-location subsections at Steps 10, 11, 12, 13, 14 all flipped from 🔴 stub to ✅ implemented with file references to the actual implementation.

## Style requirement compliance

- [x] Single self-contained HTML (no external CSS or JS except MathJax CDN).
- [x] Dark mode default with light toggle via shared `rvsmr-theme` localStorage key.
- [x] MathJax v3 via cdn.jsdelivr.net for formulas.
- [x] Inline SVG for conceptual diagrams.
- [x] PingFang SC → Microsoft YaHei → Noto Sans CJK SC → sans font stack.
- [x] Color palette: bg `#0d1117`, fg `#e6edf3`, accent `#58a6ff`; unique accent per step group.
- [x] Each step is `<details open>` collapsible.
- [x] Sticky TOC on left, collapses to hamburger button under 760px.
- [x] Citation cards with amber `#d29922` left border.
- [x] Mobile-friendly: container `max-width: 1200px`, grid collapses to single column under 760px.

## Hard constraint compliance

- [x] rvMR R package not modified.
- [x] STATUS_TEAM1.md not touched (read-only).
- [x] test_run_*/ directories not touched.
- [x] All citations verified against `citation_audit_2026-05-27.md`.

End of report.
