# Run a fast subset of scenarios to anchor the strong-IV end of the
# coverage-vs-lambda curve. Covers strong-IV F-sweep cells (F=2,5,10,20),
# pleiotropy sweep (5 cells at F=20), and the special cells (anchors,
# confounder, LD, overlap). Skips the slow F=0.25, F=0.5, F=1 cells whose
# data is already in the log.
source("run_tests_v2.R")

# Override scenarios to only fast / non-slow cells.
fast_keep <- c(
  "Fsweep_F2", "Fsweep_F5", "Fsweep_F10", "Fsweep_F20",
  "Plei_mult0", "Plei_mult0.5", "Plei_mult1", "Plei_mult2", "Plei_mult5",
  "Conf_strong", "LD_xx",
  "Anchor_null", "Anchor_strong", "Overlap_honest"
)
scenarios_v2_fast <- scenarios_v2[fast_keep]

N <- 300L
SEEDS <- master_seeds_v2

cat(sprintf("Running anchors+pleio+special cells with n_reps=%d, %d cells, %d seeds\n",
            N, length(scenarios_v2_fast), length(SEEDS)))

results <- list(); summaries <- list(); per_seed <- list()
t0 <- Sys.time()
for (i in seq_along(scenarios_v2_fast)) {
  nm <- names(scenarios_v2_fast)[i]
  scn <- scenarios_v2_fast[[nm]]
  cat(sprintf("\n=== %s ===\n", scn$label))
  # use scenario index from the FULL scenarios list to keep seeds consistent
  i_full <- match(nm, names(scenarios_v2))
  results[[nm]] <- list(); per_seed[[nm]] <- list()
  for (s_idx in seq_along(SEEDS)) {
    ms <- SEEDS[s_idx]
    cat(sprintf("  seed %d (idx %d) ... ", ms, s_idx))
    tic <- Sys.time()
    df <- run_one_scenario_seed(scn, n_reps = N, master_seed = ms,
                                scenario_idx_offset = i_full)
    sm <- summarize_one_seed(scn, df)
    results[[nm]][[s_idx]] <- df
    per_seed[[nm]][[s_idx]] <- sm
    toc <- Sys.time()
    cat(sprintf("done (%.0fs): cov_AR=%.3f cov_IVW=%.3f cov_TSLS=%.3f F_mean=%.2f n_bd=%d\n",
                as.numeric(toc - tic, units = "secs"),
                sm$coverage_AR, sm$coverage_naive, sm$coverage_tsls,
                sm$F_mean, sm$n_bounded))
  }
  summaries[[nm]] <- aggregate_seeds(per_seed[[nm]])
}
t1 <- Sys.time()
dur <- as.numeric(t1 - t0, units = "secs")
cat(sprintf("\nTotal duration: %.1fs\n", dur))

R <- list(results = results, summaries = summaries, per_seed = per_seed,
          n_reps = N, master_seeds = SEEDS, duration_sec = dur)
saveRDS(R, "results_v2_anchors.rds")
write_results_v2_md(R, "results_v2_anchors.md")
cat("Wrote results_v2_anchors.rds and results_v2_anchors.md\n")
