# Writer notes — `main.tex` (rvSMR methods writeup)

## Final document

- Path: `/home/francisfenglu4/projects/rvSMR_Math/main.tex`
- Lines: 226 (full file including preamble + 25-entry bibliography)
- Estimated compiled length: **~4.5 pages** at 11pt, 1in margins. Body has 10 short sections, 12 display equations, 1 table, 1 assumption block, 2 remarks; bibliography is 25 entries which typically runs ~3/4 page in the chosen format. Body should fit in ~3.5--3.75 pages; total ~4.25--4.5 pages — within the 4--5 target.
- Compilation: not verified locally (no `pdflatex` on host). Uses only the requested packages: `geometry`, `amsmath`, `amssymb`, `amsthm`, `mathtools`, `booktabs`, `hyperref`, `natbib`.

## Judgment calls

1. **Bibliography count = 25, above the "~15 core" guidance.** Reason: the spec body explicitly names ~22 citations across the 10 sections (Madsen-Browning, STAAR, Wu SKAT, Robins, Didelez-Sheehan, Imbens-Angrist, Wang-Kang, Fieller, Lee-McCrary-Moreira-Porter, Patel-Lane-Burgess, Sargan, Hansen, Davies, Kuonen, Zhu 2016, Cochran, Sun KY 2024, Dhindsa, Cuomo, Cinelli-Hazlett 2025, Cinelli-Hazlett 2020, Swanson-VanderWeele, Wang-Tchetgen, Karczewski Genebass, Ferkingstad deCODE). Cutting to 15 would force dropping verified citations the body actually invokes; chose to include them all rather than create unreferenced anchors. If a stricter trim is wanted, the most droppable are Wu SKAT 2011, Ferkingstad 2021, Karczewski 2022 (Genebass), and possibly Cochran 1954.

2. **Section 8 (cell-type concordance) kept to 2 sentences as specified.** Did not expand the $Q$-statistic algebra; reused (\ref{eq:cochranQ}) by reference.

3. **Robins 1994 venue:** correctly cited as *Communications in Statistics — Theory and Methods* 23(8):2379, per `citation_audit_2026-05-27.md` §"Immediate fix 1". Not *Biometrics*.

4. **Cinelli–Hazlett dual cite:** 2025 *Biometrika* asaf004 is primary in §9; 2020 *JRSS-B* 82(1):39 is secondary (framework origin). Matches audit.

5. **Sun KY (not BB) 2024 RGC ME:** corrected throughout, per audit additional-verification §2.

6. **Four-CI-shape table:** I verified the formulas $A=\hat b_x^2-c\,\mathrm{SE}_x^2$, $B=-2(\hat b_x\hat b_y-c\rho\,\mathrm{SE}_x\mathrm{SE}_y)$, $C=\hat b_y^2-c\,\mathrm{SE}_y^2$ against `rvMR/R/mrAR.R:115-117` and briefing §4.3. The four regimes (bounded / empty / disconnected / whole-line) match both sources exactly. Sign convention on $B$: the doc writes $B\beta_0$ in (\ref{eq:quadK1}) but stores $B=-2(\ldots)$; the inequality $A\beta_0^2 + B\beta_0 + C \le 0$ is internally consistent with the explicit $A,B,C$ definitions.

7. **HEIDI-rv eigenvalue note:** explicit one-sentence justification of why $V_\delta$'s non-zero eigenvalues — and not $V_\delta^{+}V_\delta$'s — are the correct Davies/Kuonen weights, per `heidi_rv.R:21-29`.

8. **Estimand formula (\ref{eq:estimand}):** matches briefing §3 Delta 2 exactly. Used $\pi_j\propto w_j^2 p_j(1-p_j)\alpha_j^2$ (the IVW form rvSMR commits to), not the first-power $\pi_j^{(1)}$ for a single pre-specified burden.

9. **No long introduction; no Methods/Discussion division.** Kept the abstract to 3 sentences; went straight into setup.

## Open citation TODOs (from HANDOVER §5 / audit §"Use-mismatch")

These were left unresolved because they require Francis's decision, not the writer's. The current draft routes around them safely:

- **(a) Delta-method variance source for the Wald ratio.** Currently the draft only invokes the AR pivot (which sidesteps delta-method calibration), so no delta-method citation is needed in §4--5. If a future revision adds an explicit delta-SE expression, decide between Burgess--Thompson 2017 textbook vs Thomas et al.\ 2007 vs Rothman--Greenland.
- **(b) CAST vs CMC.** Not invoked here — Madsen-Browning 2009 covers the linear weighted-sum burden. If CAST/CMC siblings are wanted in §2, decide attribution: CAST = Morgenthaler-Thilly 2007 *Mutat Res* 615:28; CMC = Li-Leal 2008 *AJHG* 83:311.
- **(c) STAAR 2020 vs 2022.** Used Li X 2020 *Nat Genet* 52:969 (annotation-weighting methodology). Switch to Li Z 2022 *Nat Methods* (STAARpipeline) only if §2 cites the pipeline rather than the weights.
- **(d) Sign-concordance source.** Section 7 uses Cochran-$Q$ (Cochran 1954) for class concordance, which is the correct primary reference. No Han-Eskin / Stouffer / Whitlock citation is needed; that decision is now closed for this writeup.

## What was deliberately NOT included (scope control)

- No formal identification proof (sketched assumptions only — full theorem would push past 5 pages; HANDOVER §7 flags §3 proof as a separate task).
- No worked simulation; no power curve; no figure.
- No comparator-table discussion (RARE, MR-CARV, etc.) — that belongs in a separate Background section.
- No sample-overlap derivation beyond the cross-term in $V(\beta_0)$.
- No empirical case study (PCSK9 etc.).

## Verification status

- Every display equation cross-checked against `briefing_for_wei.md`, `mrAR.R`, `mrAR_multi.R`, `heidi_rv.R`, or `HANDOVER_2026-05-27.md` §10.
- Every bibliography entry's venue+year cross-checked against `citation_audit_2026-05-27.md` (lines 9--23 confirmed list; lines 27--34 dual-cite; lines 89--132 additional-verification list).
- `pdflatex` not available on the host; cannot confirm zero-warning compile. The preamble uses only the requested packages, and the file is syntactically standard LaTeX.
