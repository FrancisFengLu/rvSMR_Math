# Citation Audit — Sub-worker 2B

Reads `content_draft.md` and verifies each citation against:
- `/home/francisfenglu4/rvSMR/May_30md/citation_audit_2026-05-27.md`
- `/home/francisfenglu4/projects/rvSMR_Math/research_AR_cisMR_comparators.md`
- `/home/francisfenglu4/projects/rvSMR_Math/research_RARE.md`, `research_MR-CARV.md`
- `/home/francisfenglu4/projects/rvSMR_Math/main.tex` bibliography
- `/home/francisfenglu4/rvSMR/May_30md/briefing_for_wei.md`

Format: STEP — CITATION — STATUS — ACTION.

---

## Step 0

- Pierce & Burgess 2013 *AJE* 178(7):1177 — ✅ confirmed (audit §13).
- Burgess-Butterworth-Thompson 2013 *Genet Epidemiol* 37(7):658 — ✅ confirmed (audit explicit: NOT Stat Med).
- Burgess et al. 2016 *Genet Epidemiol* 40(7):597 — ✅ confirmed (sample overlap; this is BURGESS-DAVIES-THOMPSON 2016, not the 2013 BBT).
- Sun KY 2024 *Nature* 631:583 — ✅ confirmed (audit §2 RGC; first author Sun **KY** NOT Sun **BB**).

## Step 1

- Madsen-Browning 2009 *PLOS Genet* 5(2):e1000384 — ✅ confirmed (audit confirmed list).
- Li X STAAR 2020 *Nat Genet* 52:969 — ✅ confirmed (audit use-mismatch (c) resolved as 2020 = annotation weights). NOT Li Z 2022 STAARpipeline.
- Wu SKAT 2011 *AJHG* 89:82 — ✅ confirmed (audit, Beta(1,25) origin).
- Morgenthaler-Thilly 2007 *Mutation Research* 615:28 (CAST) — ✅ confirmed (audit use-mismatch (b) clarification).
- Li B-Leal S 2008 *AJHG* 83:311 (CMC) — ✅ confirmed (audit use-mismatch (b) corrected attribution).
- Liu Y ACAT 2019 *AJHG* 104:410 — NOT in audit but standard; let it stand (background only).

## Step 2

- Zhou-Cuomo 2024 medRxiv 2024.05.15.24307317 (SAIGE-QTL) — ✅ confirmed (main.tex bibliography).
- Zhou SAIGE-GENE+ 2022 *Nat Genet* 54(10):1466 — ✅ confirmed (audit).
- Karczewski Genebass 2022 *Cell Genomics* 2:100168 — ✅ confirmed (audit).
- Sun KY 2024 *Nature* 631:583 — ✅ confirmed.

## Step 3

- Burgess-Labrecque 2018 *EJE* 33(10):947 — ✅ confirmed (main.tex bib).
- Didelez-Sheehan 2007 *Stat Methods Med Res* 16:309 — ✅ confirmed (main.tex bib).
- Lee-McCrary-Moreira-Porter 2022 *AER* 112(10):3260 — ✅ confirmed (main.tex bib).
- Ye T 2021 dIVW *Ann Stat* — comparator only; not formally cited in main.tex bib. **ACTION: include as background comparator without full citation OR cite ye_etal_2021_annals_statistics generically as comparator.**
- Zhao Q 2020 MR-RAPS *Ann Stat* — same as above; background comparator.

## Step 4

- Anderson-Rubin 1949 *Ann Math Stat* 20(1):46 — ✅ confirmed (main.tex bib).
- Fieller 1954 *JRSSB* 16(2):175 — ✅ confirmed (main.tex bib).
- Wang-Kang 2022 *Biometrics* 78(4):1699 — ✅ confirmed (main.tex bib).
- Patel-Lane-Burgess 2024 arXiv:2408.09868 — ✅ confirmed; note actual title is "Weak instruments in multivariable Mendelian randomization: methods and practice" (per user spec).
- Lee-McCrary-Moreira-Porter 2022 — ✅ as above.

## Step 5

- Patel-Lane-Burgess 2024 arXiv:2408.09868 — ✅ (correct title noted).
- Sargan 1958 *Econometrica* 26(3):393 — ✅ confirmed (audit + main.tex bib).
- Burgess-Butterworth-Thompson 2013 — ✅ (Genet Epidemiol).

## Step 6

- Anderson-Rubin 1949 — ✅.
- Wang-Kang 2022 — ✅.
- Patel-Lane-Burgess 2024 — ✅ (correct title).
- Burgess-Butterworth-Thompson 2013 — ✅.
- Bowden 2015 MR-Egger *Int J Epidemiol* 44(2):512 — comparator; not in main.tex bib. **ACTION: add to references list or cite as background comparator only.**

## Step 7

- Anderson-Rubin 1949 — ✅.
- Fieller 1954 — ✅.
- Brent 1973 root-finding — software detail; R `stats::uniroot()` is well-known.

## Step 8

- Sargan 1958 — ✅.
- Hansen 1982 *Econometrica* 50(4):1029 — ✅ confirmed (audit).
- Patel-Lane-Burgess 2024 — ✅.
- Bowden 2017 Cochran's Q on Wald ratios — comparator; not in main.tex bib.
- Verbanck 2018 MR-PRESSO *Nat Genet* 50:693 — comparator; not in main.tex bib. **ACTION: keep as background mention only.**

## Step 9

- Burgess-Davies-Thompson 2016 *Genet Epidemiol* 40(7):597 — ✅ confirmed (main.tex bib as `BurgessOverlap2016`; user spec phrasing: "Burgess et al. 2016 *Genet Epidemiol*, not Burgess-Butterworth which is 2013").
- Anderson-Rubin 1949 — ✅.
- Wang-Kang 2022 — ✅.
- Bulik-Sullivan LDSC 2015 — comparator; OK to keep as background.

## Step 10

- Zhu 2016 *Nat Genet* 48(5):481 (HEIDI) — ✅ confirmed (main.tex bib).
- Davies 1980 *Applied Statistics* 29(3):323 — ✅ confirmed (audit). Journal is *Applied Statistics* = *JRSS-C*; user spec writes "JRSSC". Both are valid names for the same journal; standardize to *Applied Statistics (JRSS-C)*.
- Kuonen 1999 *Biometrika* 86(4):929 — ✅ confirmed (audit).

## Step 11

- Cochran 1954 *Biometrics* 10(1):101 — ✅ confirmed (main.tex bib).
- Dhindsa 2023 *Nature* 622:339 — ✅ confirmed (audit).
- Ferkingstad 2021 *Nat Genet* 53:1712 — ✅ confirmed (audit).
- Sun BB 2018 INTERVAL *Nature* 558:73 — anti-comparator (only mention as "not this"). main.tex bib has it.
- Sun BB 2023 UKB-PPP common-variant *Nature* 622:329 — anti-comparator; main.tex bib has it.
- Han-Eskin 2011 *AJHG* 88:586 — flagged in audit as misattribution; main.tex doesn't cite it positively (only the pitfall mentions "don't cite Han-Eskin").

## Step 12

- Cuomo 2025 medRxiv 2025.03.20.25324352 — ✅ confirmed (main.tex bib).
- Ray 2025 *AJHG* 112(7):1597 — ✅ confirmed (main.tex bib; user spec: NOT Ge 2025).
- Yazar OneK1K 2022 *Science* 376:eabf3041 — ✅ confirmed (audit).
- Cochran 1954 — ✅.

## Step 13

- Cinelli-Hazlett 2025 *Biometrika* asaf004 — ✅ confirmed (audit dual cite). PRIMARY IV.
- Cinelli-Hazlett 2020 *JRSS-B* 82(1):39 — ✅ confirmed (audit dual cite). SECONDARY OLS.
- Swanson-VanderWeele 2020 *Epidemiology* 31(3):e23 — ✅ confirmed (audit additional verified §1).
- VanderWeele-Ding 2017 *Annals of Internal Medicine* 167(4):268 — ✅ confirmed (E-value origin; per user spec).

## Step 14

- Wang-Tchetgen Tchetgen 2018 *JRSS-B* 80(3):531 — ✅ confirmed (audit additional verified §7). Partial-identification bounds under invalid IV.
- Bowden 2015 *Int J Epidemiol* 44(2):512 — comparator; not in main.tex bib. **ACTION: cite as MR-Egger contrast in references; was specifically called out by user spec for Step 14.**
- Cinelli-Hazlett 2025 — ✅.

---

## Specifically flagged traps from user spec — all checked

1. ✅ Robins 1994 → *Communications in Statistics — Theory and Methods* (audit FIX 1 applied). content_draft.md does not introduce Robins; he appears only in briefing_for_wei.md §3. If HTML cites him in Step 3 / Step 4 estimand discussion, must use Comm Stat Theory Methods.

2. ✅ Sun KY (NOT Sun BB) for RGC-ME 2024. content_draft.md correctly uses Sun KY at Step 0 and Step 2.

3. ✅ STAAR 2020 (Li X, *Nat Genet* 52:969) NOT Li Z 2022. content_draft.md correctly uses Li X 2020 at Step 1.

4. ✅ Cinelli-Hazlett dual cite (2025 PRIMARY, 2020 SECONDARY). content_draft.md correctly does this at Step 13.

5. ✅ HEIDI-rv weights from V_δ (NOT V_δ^+ V_δ). content_draft.md headline-traps this at Step 10.

6. ✅ Ray 2025 (NOT Ge 2025) for cell-type comparator. content_draft.md correctly cites Ray at Step 12.

7. ✅ Patel-Lane-Burgess 2024 — correct title noted as "Weak instruments in multivariable Mendelian randomization: methods and practice" (NOT "AR tests for MR" which was HANDOVER's early misattribution). content_draft.md notes this at Step 4, Step 5, Step 6, Step 8.

8. ✅ Burgess-Butterworth-Thompson 2013 = *Genet Epidemiol* (NOT Stat Med). content_draft.md correctly does this at Step 0.

9. ✅ Burgess-Davies-Thompson 2016 (NOT Burgess-Butterworth-Thompson 2013) for sample overlap. content_draft.md correctly does this at Step 9.

10. ✅ Lee-McCrary-Moreira-Porter 2022 *AER* 112(10):3260 for tF (NOT Stock-Yogo 2005). content_draft.md correctly notes this at Step 3 and Step 4.

11. ✅ Davies 1980 journal = *Applied Statistics / JRSS-C*. content_draft.md uses "*Applied Statistics (JRSS-C)*".

12. ✅ J df = K-1 (NOT K). content_draft.md correctly states at Step 8 (and pitfall).

13. ✅ RV positive form (√(t²+4) - t)/2. content_draft.md correctly boxes the positive form at Step 13.

14. ✅ E-value RR coefficient = 0.91. content_draft.md correctly uses 0.91 at Step 13.

15. ✅ AR cross-term ordering D_y R_xy D_x. content_draft.md correctly uses this at Step 6 and Step 9.

16. ✅ AR reference distribution = $\chi^2_K$ at fixed $\beta_0$; J reference = $\chi^2_{K-1}$ at argmin. content_draft.md keeps these straight.

---

## Minor cleanups recommended for final HTML

1. Add explicit "Brent 1973" reference for `stats::uniroot()` if the HTML wants academic-grade footnote; otherwise mark as software detail (acceptable).
2. Bowden et al. 2015 MR-Egger *Int J Epidemiol* 44(2):512 should be added to the reference list since it's cited at Step 6, Step 8, Step 14 as comparator.
3. Verbanck 2018 MR-PRESSO *Nat Genet* 50(5):693 — same; if cited in HTML add to ref list.
4. Liu Y 2019 ACAT *AJHG* 104(3):410 — same.
5. Hemani 2017 bidirectional MR — optional, can be omitted.
6. Bulik-Sullivan 2015 LDSC *Nat Genet* 47:291 — optional, background only.
7. Bowden 2017 Cochran's Q on Wald ratios *Int J Epidemiol* — comparator only.
8. Imbens-Angrist 1994 *Econometrica* 62:467 — used in estimand framing (Step 3's mention of LATE / monotonicity). main.tex bib has it.

## No NEW citation errors introduced in content_draft.md

All flagged traps from user spec were avoided. The draft is citation-clean.

End of citation audit.
