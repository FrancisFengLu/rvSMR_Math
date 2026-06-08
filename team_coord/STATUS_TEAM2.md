# TEAM 2 STATUS — pedagogical HTML walkthrough

**Lead deliverable:** `/home/francisfenglu4/projects/rvSMR_Math/algorithm_paper_walkthrough.html`

## Milestone 1 / 100% — HTML delivered + sub-worker drafts + erratum integrated (2026-06-08)

### Deliverable

- `algorithm_paper_walkthrough.html` — 1900 lines, single self-contained file.
- Structure: intro + pipeline overview SVG + 14 step sections (each w/ 6 sub-subsections) + summary + references.
- Style: dark mode default with light toggle (shared `rvsmr-theme` localStorage), MathJax v3 CDN, PingFang SC font stack, sticky TOC (collapses on mobile <760px).
- 6 inline SVGs: pipeline overview, burden aggregation (Step 1), IV DAG (Step 3), 4 CI shapes (Step 4), AR inversion (Step 7), decision flowchart (Step 14).
- 19 citation cards (amber `#d29922` border per spec).
- Color-coded section banners: Part I input (teal), Part II K=1 (blue), Part III multi-IV (purple), Part IV over-id (amber), Part V decision (green).

### Sub-worker drafts in team2_drafts/

- `content_draft.md` — sub-worker 2A (writer)
- `citation_audit.md` — sub-worker 2B (verifier; all 16 trap items checked, no errors found)
- `pedagogy_review.md` — sub-worker 2C (reviewer; 7 SVGs recommended, 6 implemented in final)

### Cross-team integration

- Read Team 1's status at MS1 / 50%: **all 5 stubs implemented; 160 tests PASS, 0 FAIL.**
- Updated HTML's "代码位置" subsections at Steps 10, 11, 12, 13, 14: all switched from 🔴 stub to ✅ implemented.

### CRITICAL math erratum integrated (from Team 1's flag)

Team 1 caught a mathematical inconsistency in main.tex Step 10 eq (32-34):
- Literal main.tex says: $T = \delta^\top V_\delta^+ \delta$ with Davies weights = nonzero eigenvalues of $V_\delta$.
- These two pairings are inconsistent. Mahalanobis form $\delta^\top V_\delta^+ \delta$ follows plain $\chi^2_{\mathrm{rank}(V_\delta)}$ (NOT generalized $\chi^2$); the un-normalized $T = \delta^\top \delta$ is what has the generalized $\chi^2$ distribution with weights = eigenvalues of $V_\delta$.
- Team 1 empirically verified: literal main.tex pairing gives ~99% Type-I error at nominal 5% in m=4 simulation. Corrected pairing ($T = \delta^\top \delta$ + Davies weights = eig($V_\delta$)) gives 4.8% Type-I error and uniform null p-values.
- rvMR's `heidi_rv()` implements the corrected form ($T = \delta^\top \delta$) and reports Mahalanobis form as a diagnostic sister ($T_{\rm mahalanobis}$, $p_{\rm mahalanobis}$).
- HTML Step 10 has been updated: math block now boxes the correct form, intuition paragraph and erratum callout box explain the inconsistency in main.tex, pitfall section rewritten to flag the corrected headline trap.

### Questions for Team 1 — none

Citation chain locked: Davies 1980, Kuonen 1999, Zhu 2016 all consistent with the corrected math.

### Hard-constraint checklist

- [x] Did NOT modify rvMR R package (Team 1's territory).
- [x] Did NOT touch STATUS_TEAM1.md.
- [x] Did NOT touch test_run_*/ directories.
- [x] Cited only verified papers (per `citation_audit_2026-05-27.md`).
- [x] Mobile-friendly grid + hamburger TOC.
- [x] All 7 trap items in citation_audit covered correctly.

---
2026-06-08 04:30 — DELIVERED. All 14 steps render; HTML opens at file:///home/francisfenglu4/projects/rvSMR_Math/algorithm_paper_walkthrough.html.
