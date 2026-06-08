# rvMR Round-3 / Worker A simulation results

- Replicates per scenario per seed: **1000**
- Master seeds: 20260603, 99887766, 1234567
- Total wall time: 7431.3s
- R version: R version 4.1.2 (2021-11-01)
- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR (Worker B territory; not modified)

Comparators reported:
  - **IVW-of-ratios** (delta-method SE): per-mask Wald ratio + IVW mean. Partially weak-IV-buffered (Round 2 confirmed over-cover at low lambda).
  - **scalar IVW (`ivw_summary`)**: (b_x' W b_x)^{-1} b_x' W b_y with W = diag(1/SE_y^2). Equivalent to `MendelianRandomization::mr_ivw(method="default")`. NOT 2SLS (CRITIQUE_v2 §S2.v2.3 rename from v2 `tsls_*`).

## (a) Homogeneous-sign weak-IV F-sweep (the cells missing from Round 2)

alpha_k = +sqrt(F)*SE_x for all k. K=3, lambda_joint = K*F. This is the row that
directly extends the Round-2 v2 headline table at F in {0.25, 0.5, 1}.

| F | lambda | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | F_mean | n_bounded | bias|bounded |
|---:|---:|---|---|---|---:|---:|---|
| 0.25 | 0.75 | 0.956 +/- 0.004 | 0.976 +/- 0.003 | 0.901 +/- 0.006 | 1.26 | 99 | -0.3197 |
| 0.5 | 1.50 | 0.951 +/- 0.003 | 0.971 +/- 0.002 | 0.902 +/- 0.005 | 1.53 | 159 | -0.2405 |
| 1 | 3.00 | 0.954 +/- 0.006 | 0.972 +/- 0.001 | 0.924 +/- 0.004 | 2.03 | 279 | -0.1830 |

## (b) Sign-alternated alpha_k weak-IV sweep (Wang-Kang 2022 §3 style)

alpha_k alternates sign across K=3 masks (+ - +) so the pooled b_x summary statistic
is near zero in expectation. Scalar IVW (b_x' W b_x)^{-1} factor blows up at low F.
AR is sign-invariant in the moment m_k = b_y_k - beta * b_x_k and should hold at nominal.

Reference: Wang & Kang 2022 Biometrics 78(4):1699-1713, §3 (DOI 10.1111/biom.13524).

| F | lambda | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | F_mean | n_bounded | bias|bounded |
|---:|---:|---|---|---|---:|---:|---|
| 0.25 | 0.75 | 0.950 +/- 0.006 | 0.979 +/- 0.005 | 0.896 +/- 0.006 | 1.24 | 98 | -0.3410 |
| 0.5 | 1.50 | 0.943 +/- 0.003 | 0.970 +/- 0.003 | 0.904 +/- 0.002 | 1.50 | 155 | -0.2571 |
| 1 | 3.00 | 0.950 +/- 0.005 | 0.969 +/- 0.004 | 0.916 +/- 0.005 | 1.99 | 273 | -0.2068 |
| 2 | 6.00 | 0.953 +/- 0.007 | 0.957 +/- 0.004 | 0.925 +/- 0.001 | 3.02 | 520 | -0.0959 |
| 5 | 15.00 | 0.954 +/- 0.003 | 0.946 +/- 0.002 | 0.927 +/- 0.002 | 6.10 | 913 | -0.0051 |
| 20 | 60.00 | 0.957 +/- 0.000 | 0.958 +/- 0.002 | 0.941 +/- 0.002 | 20.93 | 982 | +0.0090 |

## (c) Confounder-strength sweep at F=20

Non-AR DGP via shared latent u: b_x_k = alpha_k + cs*SE_x*u + s_idio*SE_x*eps_x,
b_y_k = beta*alpha_k + cs*SE_y*u + s_idio*SE_y*eps_y; s_idio = sqrt(1-cs^2).
Inference call uses default R_xx = R_yy = I, R_xy = 0 (canonical "user does not know about u" stress test).
cs > 1 violates s_idio's domain and is reported as DGP-error (skipped).

| cs | AR cov (mean +/- SE) | IVW-ratios cov | scalar-IVW cov | rej_zero (AR) | F_mean | n_bounded | n_dgp_err |
|---:|---|---|---|---:|---:|---:|---:|
| 0.1 | 0.946 +/- 0.003 | 0.946 +/- 0.004 | 0.924 +/- 0.005 | 0.739 | 21.11 | 979 | 0 |
| 0.3 | 0.962 +/- 0.002 | 0.958 +/- 0.004 | 0.940 +/- 0.003 | 0.718 | 21.03 | 988 | 0 |
| 0.5 | 0.975 +/- 0.002 | 0.953 +/- 0.005 | 0.937 +/- 0.006 | 0.703 | 21.04 | 995 | 0 |
| 0.7 | 0.988 +/- 0.001 | 0.956 +/- 0.004 | 0.932 +/- 0.006 | 0.631 | 20.97 | 999 | 0 |
| 1 | 0.997 +/- 0.001 | 0.960 +/- 0.001 | 0.943 +/- 0.001 | 0.571 | 21.22 | 998 | 0 |
| 1.5 | NaN +/- NA | NaN +/- NA | NaN +/- NA | NaN | NaN | 0 | 3000 |

## Per-seed AR coverage (each row = scenario; columns = 3 seeds)

Seeds: 20260603, 99887766, 1234567 (spread-out, CRITIQUE_v2 §S3.v2.3).

| Scenario | seed1 (20260603) | seed2 (99887766) | seed3 (1234567) |
|---|---:|---:|---:|
| F-sweep (homog sign): F=0.25 (lambda=0.75) | 0.964 | 0.955 | 0.950 |
| F-sweep (homog sign): F=0.5 (lambda=1.50) | 0.951 | 0.955 | 0.946 |
| F-sweep (homog sign): F=1 (lambda=3.00) | 0.965 | 0.944 | 0.952 |
| Sign-alt alpha (Wang-Kang style): F=0.25 (lambda=0.75) | 0.961 | 0.944 | 0.944 |
| Sign-alt alpha (Wang-Kang style): F=0.5 (lambda=1.50) | 0.947 | 0.938 | 0.945 |
| Sign-alt alpha (Wang-Kang style): F=1 (lambda=3.00) | 0.956 | 0.941 | 0.953 |
| Sign-alt alpha (Wang-Kang style): F=2 (lambda=6.00) | 0.962 | 0.940 | 0.956 |
| Sign-alt alpha (Wang-Kang style): F=5 (lambda=15.00) | 0.959 | 0.948 | 0.956 |
| Sign-alt alpha (Wang-Kang style): F=20 (lambda=60.00) | 0.957 | 0.957 | 0.958 |
| Confounder sweep: cs=0.1, F=20 | 0.944 | 0.942 | 0.953 |
| Confounder sweep: cs=0.3, F=20 | 0.963 | 0.958 | 0.966 |
| Confounder sweep: cs=0.5, F=20 | 0.971 | 0.974 | 0.979 |
| Confounder sweep: cs=0.7, F=20 | 0.987 | 0.990 | 0.986 |
| Confounder sweep: cs=1, F=20 | 0.995 | 0.997 | 0.998 |
| Confounder sweep: cs=1.5, F=20 | NA | NA | NA |

