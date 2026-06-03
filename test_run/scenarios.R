# =====================================================================
# scenarios.R
#
# Defines the 6 validation scenarios (Track 1 of VALIDATION_PLAN.md).
# Each scenario is a list with parameters consumed by
# simulate_burden_mr() (from generate_test_data.R).  The driver
# (run_tests.R) loops over `scenarios` × `n_reps`.
# =====================================================================

scenarios <- list(
  A_null = list(
    label        = "A: Null (beta=0, F=20)",
    K            = 3L,
    beta_true    = 0.0,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    R_xy_block   = 0,
    expected_coverage_target = "alpha-level Type-I in non-coverage test of beta=0"
  ),
  B_strong = list(
    label        = "B: Strong IV (beta=0.4, F=20)",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    R_xy_block   = 0,
    expected_coverage_target = "95% AR CI coverage of beta=0.4"
  ),
  C_weak = list(
    label        = "C: Weak IV (beta=0.4, F=1)",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 1,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    R_xy_block   = 0,
    expected_coverage_target = "95% coverage maintained (AR's claim to fame)"
  ),
  D_very_weak = list(
    label        = "D: Very weak IV (beta=0.4, F=0.5)",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 0.5,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    R_xy_block   = 0,
    expected_coverage_target = "95% coverage; whole_line CIs prevalent"
  ),
  E_pleiotropy = list(
    label        = "E: Pleiotropy (beta=0.4, F=20, 1/3 IVs invalid)",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 1/3,
    R_xy_block   = 0,
    expected_coverage_target = "Sargan-J power; coverage no longer guaranteed (AR misspecified)"
  ),
  F_overlap = list(
    label        = "F: Sample overlap (beta=0.4, F=20, R_xy=0.3)",
    K            = 3L,
    beta_true    = 0.4,
    F_target     = 20,
    n_x          = 10000L,
    n_y          = 10000L,
    pleiotropy_frac = 0,
    R_xy_block   = 0.3,
    expected_coverage_target = "95% coverage when R_xy supplied to mrAR_multi"
  )
)

n_reps_default <- 1000L
