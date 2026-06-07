# package_diff/ — snapshot of changed rvMR files

The rvMR R package lives at `/home/francisfenglu4/rvSMR/May_30md/rvMR/`, which is not under git version control. To keep this round's package changes auditable in this repo's commit, the two changed files are copied here verbatim:

- `mrAR_multi.R` ← `/home/francisfenglu4/rvSMR/May_30md/rvMR/R/mrAR_multi.R`
- `test-mrAR_multi.R` ← `/home/francisfenglu4/rvSMR/May_30md/rvMR/tests/testthat/test-mrAR_multi.R`

Diff summary:

- `mrAR_multi.R`: expanded `R_xy` `@param` roxygen block with the explicit index convention `R_xy[i,j] = cor(b_y_i, b_x_j)` and a "Symmetric-part invariance" paragraph documenting that AR depends only on the symmetric part of V_xy. Added explanatory inline comment at the V_xy assembly. No behavior change.
- `test-mrAR_multi.R`: appended three new `test_that` blocks (R_xy transpose invariance, symmetric-part materiality, Dy/Dx swap detection) totaling 14 new assertions. All pass.

See `../IMPLEMENTATION_NOTES_v3b.md` §1 for the full reasoning.
