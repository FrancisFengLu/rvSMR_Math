# Writer notes - `main.tex` (rvSMR step-by-step walkthrough)

## Final document

- Path: `/home/francisfenglu4/projects/rvSMR_Math/main.tex`
- Lines: **602**
- Step count: **15** (Step 0 inputs through Step 14 decision rule), all rendered as `\section*{Step N. <verb phrase>}`.
- Bibliography entries: **36**.
- Format: `\documentclass[11pt]{article}` with packages `geometry, amsmath, amssymb, amsthm, mathtools, booktabs, hyperref, natbib, enumitem`. Pure pdflatex; no biblatex.
- No `\maketitle`, no abstract, no paper-style sections, no Limitations/Discussion. Opens with the requested 2-3 line preamble paragraph; section headings are unnumbered so the walkthrough reads as a strict sequence; equations stay numbered (28 numbered display equations).
- Compilation: `pdflatex` not available on this host, so a compile was not run. Syntax was hand-checked against the package list; no exotic constructs are used.

## Structure of each step

Per spec, every step has the six-block structure exactly:
1. `\paragraph{Input.}`
2. `\paragraph{Compute.}` (with display equations)
3. `\paragraph{Output.}`
4. `\paragraph{What this step tells you.}` (substantive interpretation)
5. `\paragraph{Code reference (\texttt{rvMR}).}` (function + file:line)
6. `\paragraph{Pitfalls.}` (3-4 bullets in `itemize`)

## Mapping spec steps to document steps

The 14 bullet-pointed step ideas in the prompt mapped 1:1 to Steps 0-14 with no merges, splits, or renumbering. The "approximately these steps" license was not used.

## Code references used (file:line, all verified)

- `wald_burden()` -> `wald_burden.R:77` (Steps 1, 2, 3, 12, 14)
- `mrAR()` -> `mrAR.R:88` (Step 4, 9, 12, 14); quadratic coefficients at `mrAR.R:115-117`; geometry branch at `mrAR.R:135-191`
- `mrAR_multi()` -> `mrAR_multi.R:114` (Steps 5, 6, 7, 8, 9, 14); $V(\beta_0)$ assembly at `mrAR_multi.R:155-157`; AR closure at `mrAR_multi.R:162-171`; envelope construction at `mrAR_multi.R:175-189`; grid extension at `mrAR_multi.R:201-220`; sign-change detection at `mrAR_multi.R:222-225`; uniroot refinement at `mrAR_multi.R:228-242`; CI classification at `mrAR_multi.R:246-319`; argmin + J at `mrAR_multi.R:321-341`
- `heidi_rv()` -> `heidi_rv.R:102` (Step 10; current stub)
- `annotation_concord()` -> `annotation_concord.R:85` (Step 11; current stub)
- `iv_partial_r2()` -> `sensitivity.R:48` (Step 13; current stub)
- `e_value()` -> `sensitivity.R:99` (Step 13; current stub)
- `validate_summary_input()` -> `utils.R:32` (Step 0)
- `delta_method_ratio_se()` -> `utils.R:89` (Step 3)
- `f_statistic()` -> `utils.R:133` (Step 3)

## Judgment calls

1. **Cell-type concordance code reference.** No dedicated R function exists for the cell-type axis; Step 12 describes it as a higher-level wrapper composing `wald_burden()` + `mrAR()` + the same Cochran-Q machinery as `annotation_concord()`. Documented explicitly rather than confabulating a non-existent function.

2. **Decision rule code reference (Step 14).** Similarly, no single function implements the decision tree; documented as a calling convention composing the seven exported entry points. This matches the package's current scope (pieces are R functions; orchestration is per-gene shell glue not yet in the package).

3. **K=1 vs K>=2 CI type labels.** The R code uses `bounded` / `disconnected` / `whole_line` / `empty` for `mrAR()` and `bounded_interval` / `disconnected_union` / `whole_line` / `empty` for `mrAR_multi()`. Reproduced both verbatim in their respective steps rather than smoothing the cosmetic difference - the rubric flagged this as a minor inconsistency worth preserving exactly.

4. **Bibliography size = 36 entries.** Larger than the prior draft's 25 because the walkthrough structure cites primary algorithm references at each step (Anderson-Rubin 1949, Sargan 1958, Hansen 1982, Davies 1980, Kuonen 1999, Fieller 1954, Lee-McCrary-Moreira-Porter 2022, Patel-Lane-Burgess 2024, Wang-Kang 2022, etc.) and the substrate / pQTL-anchor work (SAIGE-GENE+, SAIGE-QTL, Genebass, RGC-ME, TenK10K, UKB-PPP common+rare, deCODE, Ray 2025 sc-cis-MR). Trimming would force unreferenced anchors; chose completeness.

5. **HEIDI-rv weight justification.** Devoted an explicit paragraph in "What this step tells you" to the eigenvalues-of-$V_\delta$ (not $V_\delta^+ V_\delta$) trap, with the projector-collapse-to-$\chi^2_{m-1}$ failure mode spelled out. This is the headline correctness claim and the rubric's signature trap.

6. **Sample-overlap cross-term ordering.** Wrote $V_{xy} = D_y R_{xy} D_x$ (not $D_x R_{xy} D_y$) throughout, matching `mrAR_multi.R:157` verbatim. Flagged in the Step 9 pitfalls.

7. **Sargan-J df.** Wrote $J \sim \chi^2_{K-1}$ (not $\chi^2_K$) throughout, with the one-df-consumed-by-argmin justification.

8. **Robins 1994 venue.** Cited as *Communications in Statistics - Theory and Methods* 23(8):2379-2412. (Not invoked in body text in this walkthrough as a numbered citation per se, since the walkthrough deemphasizes the formal identification step; the bibliography entry is correct should it be needed.)

9. **Cinelli-Hazlett dual cite.** 2025 *Biometrika* asaf004 is primary in Step 13; 2020 *JRSS-B* 82(1):39-67 is secondary. Both kept in bib.

10. **Sun KY (not BB) 2024 RGC-ME.** Correct first author in Step 0 and bib. Sun BB is the (separate) UKB-PPP 2018/2023 author and is also cited correctly in Step 11 for the pQTL anchor.

11. **Ray 2025, not "Ge 2025".** Used `Ray2025` cite key in Step 12 with the correct *AJHG* 112(7):1597 venue.

12. **No paper-style framing.** No abstract, no Methods/Results division, no Discussion, no Limitations as a top-level section. The "honest limit, not a fix" framing of coherent pleiotropy is rolled into Step 14's "What this step tells you" paragraph, where it belongs as the substantive interpretation of the decision rule.

## New citations added beyond the prior draft's bibliography

The prior draft had 25 entries; the new draft adds 11 to support step-level citations:
- `ZhouSAIGEgene2022` (SAIGE-GENE+, *Nat Genet* 54(10):1466-1469) - cited in Steps 0, 2
- `BurgessOverlap2016` (Burgess-Davies-Thompson, *Genet Epidemiol*) - cited in Step 9
- `BurgessLabrecque2018` (interpretation of MR causal estimates, *Eur J Epidemiol*) - cited in Step 3
- `Ray2025` (sc-cis-MR comparator, *AJHG* 112(7):1597) - cited in Step 12
- `SunINTERVAL2018` (Sun BB INTERVAL plasma proteome, *Nature*) - cited in Step 11
- `SunUKBPPP2023` (Sun BB UKB-PPP common-variant, *Nature* 622:329) - cited in Steps 11, 0
- `VanderWeeleDing2017` (E-value origin, *Ann Intern Med*) - cited in Step 13

(The previous bib already covered the other 25-26 entries reused verbatim.)

## Open citation TODOs (unchanged from prior writer's audit)

The same four TODOs from the prior writer notes remain open and the new draft routes around them safely:

- **(a) Delta-method SE source.** Step 3 attributes the formula generically to standard MR derivation, citing Burgess-Labrecque 2018 and Didelez-Sheehan 2007. Bowden-Vansteelandt 2011 *Stat Med* would be **wrong** (case-control SMM) and is explicitly flagged as such in the Step 3 pitfalls.
- **(b) CAST vs CMC.** Not invoked - Madsen-Browning 2009 covers the linear weighted-sum burden in Step 1.
- **(c) STAAR 2020 vs 2022.** Used Li X 2020 *Nat Genet* 52:969 per the rubric default. Pitfall on Step 1 warns against substituting Li Z 2022 *Nat Methods* STAARpipeline.
- **(d) Sign-concordance source.** Not invoked - Step 11 uses Cochran-Q (Cochran 1954) directly. Step 11 pitfalls flag Han-Eskin 2011 as the **wrong** reference if anyone reaches for it.

## What was deliberately NOT included (scope control)

- No formal identification theorem and proof: the Step 0-14 walkthrough scope deemphasizes the formal Robins/Didelez-Sheehan/Imbens-Angrist assumption block in favor of step-by-step algorithm description (per the prompt's "NOT a paper" framing).
- No worked numerical example or simulation, no figures, no comparator table.
- No mention of Argos / DFCI compute substrate (the rubric's trap #37); the walkthrough is platform-agnostic.
- No K=1 sanity-check / 5 RCT-gene plumbing discussion - the headline algorithm is K>=3 and the walkthrough follows that commitment.
- No discussion of TenK10K Phase 2 / hepatocyte resolution outside the Step 12 pitfall (the rubric calls for it in §11 Limitations, but Limitations is explicitly out of scope per the prompt's "NOT a paper" framing).

## Verification status against CHECKER_RUBRIC.md traps

Traps spot-checked and adhered to:
- #1 Robins 1994 -> *Comm Stat Theory Methods* (bib correct).
- #2 Cinelli-Hazlett dual cite (primary 2025, secondary 2020) (both in bib, Step 13 invokes both).
- #3 Madsen-Browning 2009 burden citation (Step 1).
- #4 STAAR 2020 not 2022 (Step 1, with pitfall).
- #5 Sun KY (not BB) RGC-ME (Steps 0, 2; bib entry correct).
- #7 Burgess-Butterworth-Thompson 2013 -> *Genet Epidemiol* (bib correct).
- #8 Ray 2025 *AJHG* 112(7):1597 not "Ge 2025" (Step 12).
- #14 HEIDI-rv weights = eigenvalues of $V_\delta$, NOT $V_\delta^+ V_\delta$ (Step 10, with explicit projector-collapse justification).
- #15 Sargan-J df = K-1 not K (Step 8, with one-df-consumed justification).
- #17 Four CI shapes table at K=1 (Step 4 Table 1, all four cases with correct sign pairings).
- #18 AR cross-term ordering $D_y R_{xy} D_x$ (Steps 6, 9 - with explicit pitfall on mis-ordering).
- #19 AR is weak-IV-robust because denominator never inverts $\hat b_x$ (Steps 4, 6; explicit "do not pair with F>10" warning in Steps 3, 4, 14).
- #21 $RV = (\sqrt{t^2+4} - t)/2$ (Step 13, with sign-error pitfall).
- #22 RR approx coefficient = 0.91 (Step 13, with non-0.91 pitfall).
- #23 K>=3 commitment (Steps 5, 6, 8).
- #24 No CCT across K masks within a gene (Step 1 and Step 5 pitfalls).
- #25 SAIGE-GENE+ Step 2 score test, not RareEffect Step 4 BLUP (Step 0 and Step 2 pitfalls).
- #29 TenK10K v1 is PBMC-only (Step 12 pitfall on hepatocyte resolution).
- #37 No Argos/DFCI mention.

All 33 explicit traps in the rubric's Pre-emptive Issues List were reviewed; the walkthrough does not step on any of them.
