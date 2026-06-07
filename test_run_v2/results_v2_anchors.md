# rvMR Round-2 simulation results

- Replicates per scenario per seed: **300**
- Master seeds: 20260603, 20260604, 20260605
- Total wall time: 1163.8s
- R version: R version 4.1.2 (2021-11-01)
- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR

## Headline: coverage vs lambda (F-sweep)

`lambda_joint = K * F_target` at K=3, R_xx=I (CRITIQUE S1.4).
Two comparators are reported:
  - **IVW-of-ratios** (CRITIQUE S1.1 canonical form): per-mask Wald ratio
    + delta-method SE + IVW. Partially weak-IV-buffered via SE inflation.
  - **TSLS** (summary form): (b_x' W b_x)^{-1} b_x' W b_y, W=diag(1/SE_y^2).
    Equivalent to scalar IVW from `MendelianRandomization::mr_ivw`.
    This is the comparator that Wang-Kang Fig 6 reports collapsing.

| F | lambda | AR cov (mean +/- SE) | IVW cov (mean +/- SE) | TSLS cov (mean +/- SE) | F_mean | n_bounded | bias|bounded |
|---:|---:|---|---|---|---:|---:|---|
| 2 | 6.0 | 0.943 +/- 0.002 | 0.964 +/- 0.008 | 0.924 +/- 0.009 | 2.92 | 152 | -0.1114 |
| 5 | 15.0 | 0.944 +/- 0.003 | 0.948 +/- 0.003 | 0.931 +/- 0.003 | 5.88 | 267 | +0.0019 |
| 10 | 30.0 | 0.944 +/- 0.003 | 0.941 +/- 0.009 | 0.918 +/- 0.017 | 10.94 | 295 | +0.0150 |
| 20 | 60.0 | 0.949 +/- 0.005 | 0.942 +/- 0.006 | 0.921 +/- 0.019 | 21.10 | 295 | +0.0005 |

## Pleiotropy magnitude sweep (at F=20, 1/3 invalid)

| pleio mult | AR cov mean | AR cov SE | J<0.05 rate | F_mean |
|---:|---:|---:|---:|---:|
| 0 | 0.944 | 0.0062 | 0.039 | 21.21 |
| 0.5 | 0.949 | 0.0087 | 0.059 | 21.05 |
| 1 | 0.901 | 0.0040 | 0.092 | 20.87 |
| 2 | 0.698 | 0.0062 | 0.211 | 20.79 |
| 5 | 0.020 | 0.0069 | 0.822 | 20.87 |

## Anchor / confounder / LD / overlap cells

| Scenario | AR cov | IVW cov | TSLS cov | rej_zero (AR) | J<0.05 | F_mean | n_bd |
|---|---:|---:|---:|---:|---:|---:|---:|
| Confounder: cs=0.5, beta=0.4, F=20 | 0.987 +/- 0.003 | 0.971 | 0.947 | 0.679 | 0.012 | 20.90 | 299 |
| LD between masks: R_xx rho=0.3, beta=0.4, F=20 | 0.950 +/- 0.009 | 0.936 | 0.910 | 0.722 | 0.050 | 21.06 | 295 |
| Anchor A: Null (beta=0, F=20) | 0.952 +/- 0.004 | 0.960 | 0.937 | 0.048 | 0.030 | 20.95 | 296 |
| Anchor B: Strong IV (beta=0.4, F=20) | 0.962 +/- 0.007 | 0.958 | 0.932 | 0.750 | 0.030 | 20.97 | 297 |
| Honest overlap (full R_xx,R_yy,R_xy block, rho=0.3, F=20) | 0.957 +/- 0.011 | 0.944 | 0.913 | 0.531 | 0.044 | 21.10 | 295 |

## Per-seed AR coverage (each row = scenario; columns = 3 seeds)

| Scenario | seed1 | seed2 | seed3 |
|---|---:|---:|---:|
| Fsweep F=2 (lambda=6.0) | 0.947 | 0.940 | 0.943 |
| Fsweep F=5 (lambda=15.0) | 0.940 | 0.943 | 0.950 |
| Fsweep F=10 (lambda=30.0) | 0.943 | 0.950 | 0.940 |
| Fsweep F=20 (lambda=60.0) | 0.950 | 0.940 | 0.957 |
| Plei sweep: pleio=0*SE_y, 1/3 invalid, F=20 | 0.940 | 0.957 | 0.937 |
| Plei sweep: pleio=0.5*SE_y, 1/3 invalid, F=20 | 0.963 | 0.933 | 0.950 |
| Plei sweep: pleio=1*SE_y, 1/3 invalid, F=20 | 0.903 | 0.893 | 0.907 |
| Plei sweep: pleio=2*SE_y, 1/3 invalid, F=20 | 0.690 | 0.693 | 0.710 |
| Plei sweep: pleio=5*SE_y, 1/3 invalid, F=20 | 0.010 | 0.017 | 0.033 |
| Confounder: cs=0.5, beta=0.4, F=20 | 0.980 | 0.990 | 0.990 |
| LD between masks: R_xx rho=0.3, beta=0.4, F=20 | 0.967 | 0.937 | 0.947 |
| Anchor A: Null (beta=0, F=20) | 0.943 | 0.957 | 0.957 |
| Anchor B: Strong IV (beta=0.4, F=20) | 0.950 | 0.963 | 0.973 |
| Honest overlap (full R_xx,R_yy,R_xy block, rho=0.3, F=20) | 0.953 | 0.977 | 0.940 |

