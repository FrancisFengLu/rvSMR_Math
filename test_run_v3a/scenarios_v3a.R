# =====================================================================
# scenarios_v3a.R  (Round 3 / Worker A)
#
# Round 3 scope (Worker A):
#   (a) F-sweep at F in {0.25, 0.5, 1} -- the three weak-IV cells
#       CRITIQUE_v2 §S1.v2.1 identified as missing from Round 2.
#       Homogeneous-sign alpha_k (Round 2 baseline) so the row is
#       directly comparable with the Round-2 v2 table.
#   (b) Sign-alternated F-sweep at F in {0.25, 0.5, 1, 2, 5, 20}.
#       New `dgp = "signalt"` from generate_test_data_v3a.R. This is
#       the canonical Wang-Kang 2022 §3 / Fig 6 weak-IV stress
#       scenario for non-robust summary-IVW. Side-by-side with (a)
#       gives the homogeneous-vs-alternating sign comparison the
#       CRITIQUE_v2 §S1.v2.2 recommendation asked for.
#   (c) Confounder strength sweep at F=20, cs in {0.1, 0.3, 0.5, 0.7,
#       1.0, 1.5}. Round 2 measured only cs=0.5 (AR coverage 0.987,
#       over-cover). This sweep is the genuine non-AR-DGP finding
#       worth characterising (CRITIQUE_v2 §S2.1 / §S2.v2.5 follow-up).
#       Note that cs=1.5 > 1 will trigger the abs(cs)<=1 stopifnot in
#       the generator; we therefore include cs=1.0 as the upper
#       admissible point and report cs=1.5 as a "skipped" cell.
#
# Notes:
#   - All scenarios use K=3, n_x=n_y=10000, beta_true=0.4.
#   - n_reps_default_v3a = 1000 per scenario per master seed (Round 3
#     spec).
#   - Three master seeds for cross-seed SE (replace Round-2 near-
#     adjacent seeds with a more spread-out triple, CRITIQUE_v2 §S3.v2.3
#     recommendation).
# =====================================================================

scenarios_v3a <- list()

# --- (a) Homogeneous-sign weak-IV F-sweep (the missing cells) --------
F_grid_weak <- c(0.25, 0.5, 1)
for (F_t in F_grid_weak) {
  nm <- sprintf("Fsweep_F%g", F_t)
  scenarios_v3a[[nm]] <- list(
    label        = sprintf("F-sweep (homog sign): F=%g (lambda=%.2f)",
                            F_t, 3 * F_t),
    family       = "F_sweep_homog",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = F_t,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    pleio_size_mult = 0,
    dgp          = "ar",
    conf_strength = 0,
    rho_xx       = 0,
    rho_xy_diag  = 0
  )
}

# --- (b) Sign-alternated F-sweep (Wang-Kang-style collapse demo) -----
# Reproduces the actual Wang-Kang 2022 §3 weak-IV stress: alpha_k
# alternates +/-/+ across K=3 masks, so the pooled b_x summary statistic
# is near zero in expectation, amplifying the (b_x' W b_x)^{-1} factor
# in scalar IVW. AR is sign-invariant in the moment
# m_k = b_y_k - beta * b_x_k and should hold at nominal.
F_grid_alt <- c(0.25, 0.5, 1, 2, 5, 20)
for (F_t in F_grid_alt) {
  nm <- sprintf("SignAlt_F%g", F_t)
  scenarios_v3a[[nm]] <- list(
    label        = sprintf("Sign-alt alpha (Wang-Kang style): F=%g (lambda=%.2f)",
                            F_t, 3 * F_t),
    family       = "F_sweep_signalt",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = F_t,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    pleio_size_mult = 0,
    dgp          = "signalt",
    conf_strength = 0,
    rho_xx       = 0,
    rho_xy_diag  = 0
  )
}

# --- (c) Confounder-strength sweep at F=20 ---------------------------
# cs > 1 violates the generator's `abs(cs) <= 1` constraint, so we
# include cs=1.5 in the schedule but mark it as expected-to-skip
# (handled in the driver with a tryCatch and a "skipped" marker).
cs_grid <- c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5)
for (cs in cs_grid) {
  nm <- sprintf("ConfSweep_cs%g", cs)
  scenarios_v3a[[nm]] <- list(
    label        = sprintf("Confounder sweep: cs=%g, F=20", cs),
    family       = "Conf_sweep",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    pleio_size_mult = 0,
    dgp          = "confounder",
    conf_strength = cs,
    rho_xx       = 0,
    rho_xy_diag  = 0
  )
}

n_reps_default_v3a  <- 1000L
# Spread-out master seeds (CRITIQUE_v2 §S3.v2.3): break any near-adjacent
# integer correlation in Mersenne-Twister streams.
master_seeds_v3a    <- c(20260603L, 99887766L, 1234567L)
