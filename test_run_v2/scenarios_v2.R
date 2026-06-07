# =====================================================================
# scenarios_v2.R  (Round 2)
#
# Round 2 grid built to address CRITIQUE:
#   (a) F sweep 7 cells {0.25, 0.5, 1, 2, 5, 10, 20} at beta=0.4, no plei,
#       no conf. This is the headline Wang-Kang / PLB Fig 6 / Fig 2
#       coverage-vs-lambda figure. Cells are labelled by lambda_joint =
#       K * F_target (CRITIQUE S1.4 fix: lambda axis, not F).
#   (b) Pleiotropy magnitude sweep 5 cells at F=20, pleio_size_mult in
#       {0, 0.5, 1, 2, 5} * SE_y on 1/3 invalid masks (CRITIQUE S2.3).
#   (c) Confounder scenario (Tier 4, CRITIQUE S2.1): non-AR DGP via
#       shared latent u, conf_strength = 0.5.
#   (d) LD-between-masks scenario (CRITIQUE S2.2): R_xx not I, rho_xx=0.3,
#       inference call given correct R_xx.
#   (e) Anchors: A null (beta=0, F=20), B strong-IV signal (beta=0.4,
#       F=20).
#   (f) Sample-overlap honest scenario (CRITIQUE S1.5 / S3): full
#       block R_xx, R_yy, R_xy, all passed to mrAR_multi.
# =====================================================================

scenarios_v2 <- list()

# --- (a) F sweep -----------------------------------------------------
F_grid <- c(0.25, 0.5, 1, 2, 5, 10, 20)
for (F_t in F_grid) {
  nm <- sprintf("Fsweep_F%g", F_t)
  scenarios_v2[[nm]] <- list(
    label        = sprintf("Fsweep F=%g (lambda=%.1f)", F_t, 3 * F_t),
    family       = "F_sweep",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = F_t,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    pleio_size_mult = 5,
    dgp          = "ar",
    conf_strength = 0,
    rho_xx       = 0,
    rho_xy_diag  = 0
  )
}

# --- (b) Pleiotropy magnitude sweep at F=20 ---------------------------
plei_mults <- c(0, 0.5, 1, 2, 5)
for (pm in plei_mults) {
  nm <- sprintf("Plei_mult%g", pm)
  scenarios_v2[[nm]] <- list(
    label        = sprintf("Plei sweep: pleio=%g*SE_y, 1/3 invalid, F=20", pm),
    family       = "Plei_sweep",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 1/3,
    pleio_size_mult = pm,
    dgp          = "ar",
    conf_strength = 0,
    rho_xx       = 0,
    rho_xy_diag  = 0
  )
}

# --- (c) Confounder (non-AR DGP, Tier 4) ------------------------------
scenarios_v2$Conf_strong <- list(
  label        = "Confounder: cs=0.5, beta=0.4, F=20",
  family       = "Confounder",
  K            = 3L,
  beta_true    = 0.4,
  F_target     = 20,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "confounder",
  conf_strength = 0.5,
  rho_xx       = 0,
  rho_xy_diag  = 0
)
scenarios_v2$Conf_weakF <- list(
  label        = "Confounder + weak IV: cs=0.5, F=1",
  family       = "Confounder",
  K            = 3L,
  beta_true    = 0.4,
  F_target     = 1,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "confounder",
  conf_strength = 0.5,
  rho_xx       = 0,
  rho_xy_diag  = 0
)

# --- (d) LD between masks --------------------------------------------
scenarios_v2$LD_xx <- list(
  label        = "LD between masks: R_xx rho=0.3, beta=0.4, F=20",
  family       = "LD_xx",
  K            = 3L,
  beta_true    = 0.4,
  F_target     = 20,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "ld_xx",
  conf_strength = 0,
  rho_xx       = 0.3,
  rho_xy_diag  = 0
)

# --- (e) Anchors -----------------------------------------------------
scenarios_v2$Anchor_null <- list(
  label        = "Anchor A: Null (beta=0, F=20)",
  family       = "Anchor",
  K            = 3L,
  beta_true    = 0.0,
  F_target     = 20,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "ar",
  conf_strength = 0,
  rho_xx       = 0,
  rho_xy_diag  = 0
)
scenarios_v2$Anchor_strong <- list(
  label        = "Anchor B: Strong IV (beta=0.4, F=20)",
  family       = "Anchor",
  K            = 3L,
  beta_true    = 0.4,
  F_target     = 20,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "ar",
  conf_strength = 0,
  rho_xx       = 0,
  rho_xy_diag  = 0
)

# --- (f) Honest sample overlap ---------------------------------------
scenarios_v2$Overlap_honest <- list(
  label        = "Honest overlap (full R_xx,R_yy,R_xy block, rho=0.3, F=20)",
  family       = "Overlap",
  K            = 3L,
  beta_true    = 0.4,
  F_target     = 20,
  n_x          = 10000L,
  n_y          = 10000L,
  pleiotropy_frac = 0,
  pleio_size_mult = 0,
  dgp          = "overlap",
  conf_strength = 0,
  rho_xx       = 0,
  rho_xy_diag  = 0.3
)

n_reps_default_v2  <- 500L
master_seeds_v2    <- c(20260603L, 20260604L, 20260605L)
