# rvMR Phase-3 simulation results

- Replicates per scenario: **1000**
- Master seed: 20260603
- Total wall time: 874.0s
- R version: R version 4.1.2 (2021-11-01)
- rvMR loaded from /home/francisfenglu4/rvSMR/May_30md/rvMR

## Headline coverage table

| Scenario | n_ok | AR cov(β) | naive cov(β) | mean F | bias β̂ | rej β=0 (AR) | J<0.05 rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| A: Null (beta=0, F=20) | 1000 | 0.951 | 0.955 | 21.15 | +0.0005 | 0.049 | 0.038 |
| B: Strong IV (beta=0.4, F=20) | 1000 | 0.954 | 0.959 | 21.04 | +0.0016 | 0.724 | 0.047 |
| C: Weak IV (beta=0.4, F=1) | 1000 | 0.950 | 0.985 | 1.97 | +0.6057 | 0.089 | 0.015 |
| D: Very weak IV (beta=0.4, F=0.5) | 1000 | 0.962 | 0.988 | 1.48 | -0.1834 | 0.061 | 0.007 |
| E: Pleiotropy (beta=0.4, F=20, 1/3 IVs invalid) | 1000 | 0.021 | 0.318 | 21.00 | +0.5529 | 1.000 | 0.807 |
| F: Sample overlap (beta=0.4, F=20, R_xy=0.3) | 1000 | 0.943 | 0.972 | 21.06 | -0.0018 | 0.711 | 0.052 |

## CI shape distribution per scenario

### A: Null (beta=0, F=20)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.983 |
| empty | 0.017 |

### B: Strong IV (beta=0.4, F=20)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.979 |
| empty | 0.021 |

### C: Weak IV (beta=0.4, F=1)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.263 |
| disconnected_union | 0.168 |
| empty | 0.006 |
| whole_line | 0.563 |

### D: Very weak IV (beta=0.4, F=0.5)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.158 |
| disconnected_union | 0.144 |
| empty | 0.001 |
| whole_line | 0.697 |

### E: Pleiotropy (beta=0.4, F=20, 1/3 IVs invalid)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.310 |
| empty | 0.690 |

### F: Sample overlap (beta=0.4, F=20, R_xy=0.3)

| CI shape | proportion |
|---|---:|
| bounded_interval | 0.979 |
| empty | 0.021 |

## Pass / fail interpretation

- **A: Null (beta=0, F=20)**: coverage_AR = 0.951, bias = +0.0005, J-rej = 0.038. Coverage of β=0 should be ≈ 0.95 (PASS if in [0.93, 0.97]). Type-I of rejecting β=0 is 0.049. **PASS**.
- **B: Strong IV (beta=0.4, F=20)**: coverage_AR = 0.954, bias = +0.0016, J-rej = 0.047. Coverage should be ≈ 0.95; bias should be ~0. **PASS**.
- **C: Weak IV (beta=0.4, F=1)**: coverage_AR = 0.950, bias = +0.6057, J-rej = 0.015. AR's headline regime: coverage **must hold** at F≈1. **PASS** -- AR is doing its job.
- **D: Very weak IV (beta=0.4, F=0.5)**: coverage_AR = 0.962, bias = -0.1834, J-rej = 0.007. F<1: expect many whole_line CIs, but coverage should still ≥ 0.95 (conservative). **PASS** (no under-coverage).
- **E: Pleiotropy (beta=0.4, F=20, 1/3 IVs invalid)**: coverage_AR = 0.021, bias = +0.5529, J-rej = 0.807. AR no longer valid (model misspecified); but Sargan-J should reject (high J p<0.05 rate). **PASS** (J detects pleiotropy).
- **F: Sample overlap (beta=0.4, F=20, R_xy=0.3)**: coverage_AR = 0.943, bias = -0.0018, J-rej = 0.052. Coverage should remain ≈ 0.95 when R_xy is supplied. **PASS**.
