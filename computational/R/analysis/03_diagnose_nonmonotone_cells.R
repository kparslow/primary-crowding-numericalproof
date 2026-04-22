# computational/R/analysis/04_diagnose_nonmonotone_cells.R
# Diagnose (logS, delta) cells where the best response a_star(v) is non-monotone
# (i.e., switches 2+ times as v increases).
#
# What this script does:
#   1) Reads the derived summary cutoffs_and_action_rates.csv and finds multi-switch cells
#   2) Reads the raw grid single_agent_grid.csv and (re)computes switch counts under eps tolerances
#   3) Recomputes Delta(v) on a dense v-grid for a few representative cells (to check that it's not numerical noise)
#   4) Writes a compact diagnostics CSV and saves simple plots for inspected cells
#
# Assumptions:
#   - 01_single_agent_grid.R produced columns: v, logS, delta, sigma, Delta, a_star
#   - model functions exist at computational/R/functions/model.R, including Delta_incentive()

suppressWarnings(suppressMessages({
  have_dt <- requireNamespace("data.table", quietly = TRUE)
}))

read_any <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (have_dt) return(data.table::fread(path))
  utils::read.csv(path, check.names = TRUE)
}

# ---- Paths ----
summary_file <- file.path("computational", "output", "derived", "cutoffs_and_action_rates.csv")
grid_file    <- file.path("computational", "output", "datasets", "single_agent_grid.csv")

base_out <- file.path("computational", "output", "diagnostics")
out_dir  <- base_out
fig_dir  <- file.path(base_out, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
summ <- read_any(summary_file)
grid <- read_any(grid_file)

# Ensure data.table behavior if available
if (have_dt) {
  data.table::setDT(summ)
  data.table::setDT(grid)
}

# ---- Basic checks ----
req_summ <- c("logS", "delta", "n_switches_a_star")
miss_summ <- setdiff(req_summ, names(summ))
if (length(miss_summ) > 0) stop("Summary missing: ", paste(miss_summ, collapse = ", "))

req_grid <- c("v", "sigma", "logS", "delta", "Delta", "a_star")
miss_grid <- setdiff(req_grid, names(grid))
if (length(miss_grid) > 0) stop("Grid missing: ", paste(miss_grid, collapse = ", "))

# ---- Identify multi-switch cells ----
multi <- summ[summ$n_switches_a_star >= 2, ]
multi <- multi[order(multi$n_switches_a_star, multi$logS, multi$delta), ]

cat("Total cells in summary:", nrow(summ), "\n")
cat("Cells with 2+ switches:", nrow(multi), "\n")
if (nrow(multi) == 0) {
  cat("No nonmonotone cells found. Exiting.\n")
  quit(save = "no", status = 0)
}

# Print the full list
print(multi)

# ---- Helper functions ----
count_switches_from_Delta <- function(Delta, eps = 0) {
  a <- as.integer(Delta >= eps)
  if (length(a) < 2) return(0L)
  sum(a[-1L] != a[-length(a)])
}

count_switches_from_a <- function(a) {
  a <- as.integer(a)
  if (length(a) < 2) return(0L)
  sum(a[-1L] != a[-length(a)])
}

intervals_from_indicator <- function(v, a) {
  # Returns a data.frame of contiguous intervals where a == 1:
  # columns: v_enter, v_exit
  a <- as.integer(a)
  if (length(a) == 0) return(data.frame(v_enter = numeric(0), v_exit = numeric(0)))
  
  r <- rle(a)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  
  keep <- which(r$values == 1L)
  if (length(keep) == 0) return(data.frame(v_enter = numeric(0), v_exit = numeric(0)))
  
  data.frame(
    v_enter = v[starts[keep]],
    v_exit  = v[ends[keep]]
  )
}

plot_cell <- function(sub, title, file) {
  # sub must contain v, Delta, a_star
  grDevices::png(file, width = 900, height = 700)
  op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  plot(sub$v, sub$Delta, type = "l",
       xlab = "v", ylab = "Delta = u1 - u0",
       main = paste0(title, " | Delta(v)"))
  abline(h = 0, col = "red", lty = 2)
  
  plot(sub$v, sub$a_star, type = "s",
       xlab = "v", ylab = "a_star",
       ylim = c(-0.05, 1.05),
       main = paste0(title, " | a_star(v)"))
  
  par(op)
  grDevices::dev.off()
}

# ---- Epsilon robustness on the stored grid ----
eps_vec <- c(0, 1e-12, 1e-10, 1e-8, 1e-6)

# Join raw grid to multi-switch cells
if (have_dt) {
  cells <- unique(multi[, .(logS, delta)])
  data.table::setkey(cells, logS, delta)
  data.table::setkey(grid, logS, delta)
  
  grid_multi <- grid[cells, nomatch = 0L]
  data.table::setorderv(grid_multi, c("logS", "delta", "v"))
  
  eps_dt <- data.table::rbindlist(lapply(eps_vec, function(eps) {
    tmp <- grid_multi[, .(
      n_switches_eps = count_switches_from_Delta(Delta, eps = eps),
      n_switches_astar = count_switches_from_a(a_star)
    ), by = .(logS, delta)]
    tmp[, eps := eps]
    tmp
  }))
  
  # Summarize how many cells remain multi-switch as eps increases
  cat("\nEpsilon robustness summary (stored grid):\n")
  print(eps_dt[, .N, by = .(eps, n_switches_eps)][order(eps, n_switches_eps)])
  print(eps_dt[n_switches_eps >= 2, .N, by = eps][order(eps)])
  
} else {
  # Base R version
  cells <- unique(multi[, c("logS", "delta")])
  eps_dt_list <- list()
  
  for (eps in eps_vec) {
    out <- data.frame(logS = numeric(0), delta = numeric(0),
                      n_switches_eps = integer(0),
                      n_switches_astar = integer(0),
                      eps = eps)
    
    for (i in seq_len(nrow(cells))) {
      logS0 <- cells$logS[i]
      delta0 <- cells$delta[i]
      sub <- grid[grid$logS == logS0 & grid$delta == delta0, ]
      sub <- sub[order(sub$v), ]
      
      out <- rbind(out, data.frame(
        logS = logS0,
        delta = delta0,
        n_switches_eps = count_switches_from_Delta(sub$Delta, eps = eps),
        n_switches_astar = count_switches_from_a(sub$a_star),
        eps = eps
      ))
    }
    eps_dt_list[[as.character(eps)]] <- out
  }
  eps_dt <- do.call(rbind, eps_dt_list)
}

# ---- Dense recomputation for a few representative cells ----
n_dense <- 6L
dense_n_v <- 50001L
dense_vmin <- -4
dense_vmax <- 4

source("computational/R/functions/model.R")

# Choose which cells to densify:
# - all 3-switch cells first
# - then enough 2-switch cells to reach n_dense
if (have_dt) {
  pick <- rbind(
    multi[n_switches_a_star >= 3, .(logS, delta, n_switches_a_star)],
    multi[n_switches_a_star == 2, .(logS, delta, n_switches_a_star)]
  )
  pick <- unique(pick)[order(-n_switches_a_star)][1:min(n_dense, .N)]
} else {
  pick <- multi[order(-multi$n_switches_a_star), c("logS", "delta", "n_switches_a_star")]
  pick <- pick[!duplicated(paste(pick$logS, pick$delta)), ]
  pick <- pick[seq_len(min(n_dense, nrow(pick))), ]
}

dense_results <- vector("list", nrow(pick))

refine_cell <- function(logS0, delta0, sigma0, n = 50001L, vmin = -4, vmax = 4) {
  S0 <- exp(logS0)
  v <- seq(vmin * sigma0, vmax * sigma0, length.out = n)
  Delta <- Delta_incentive(v, S = S0, delta = delta0)
  a <- as.integer(Delta >= 0)
  list(v = v, Delta = Delta, a = a,
       n_switches = count_switches_from_a(a),
       intervals = intervals_from_indicator(v, a))
}

cat("\nDense recomputation for selected cells:\n")
for (i in seq_len(nrow(pick))) {
  logS0 <- pick$logS[i]
  delta0 <- pick$delta[i]
  
  # sigma should be constant within a cell; grab from stored grid
  if (have_dt) {
    sigma0 <- grid[logS == logS0 & delta == delta0, sigma][1L]
  } else {
    sigma0 <- grid$sigma[grid$logS == logS0 & grid$delta == delta0][1]
  }
  if (is.na(sigma0)) sigma0 <- 1
  
  dense <- refine_cell(logS0, delta0, sigma0, n = dense_n_v, vmin = dense_vmin, vmax = dense_vmax)
  
  dense_results[[i]] <- list(
    logS = logS0,
    delta = delta0,
    sigma = sigma0,
    n_switches_dense = dense$n_switches,
    n_intervals_dense = nrow(dense$intervals),
    intervals_dense = dense$intervals
  )
  
  cat(sprintf("  logS=%s, delta=%s | dense switches=%d, intervals=%d\n",
              format(logS0), format(delta0), dense$n_switches, nrow(dense$intervals)))
  
  # Save plots using the stored grid (fast) for visual inspection
  if (have_dt) {
    sub <- grid[logS == logS0 & delta == delta0][order(v)]
    sub_df <- as.data.frame(sub[, .(v, Delta, a_star)])
  } else {
    sub <- grid[grid$logS == logS0 & grid$delta == delta0, ]
    sub <- sub[order(sub$v), ]
    sub_df <- sub[, c("v", "Delta", "a_star")]
  }
  
  title <- paste0("logS=", signif(logS0, 4), ", delta=", signif(delta0, 4))
  plot_file <- file.path(fig_dir, paste0("cell_logS_", signif(logS0, 4),
                                         "__delta_", signif(delta0, 4), ".png"))
  plot_cell(sub_df, title, plot_file)
}

# Flatten dense results into a compact table
dense_tbl <- do.call(rbind, lapply(dense_results, function(x) {
  # Save up to 3 intervals as columns (enough for up to 3 switches => up to 2 intervals),
  # but we allow up to 3 to be safe.
  ints <- x$intervals_dense
  get_int <- function(k, field) {
    if (nrow(ints) >= k) ints[[field]][k] else NA_real_
  }
  
  data.frame(
    logS = x$logS,
    delta = x$delta,
    sigma = x$sigma,
    n_switches_dense = x$n_switches_dense,
    n_intervals_dense = x$n_intervals_dense,
    v_enter_1 = get_int(1, "v_enter"),
    v_exit_1  = get_int(1, "v_exit"),
    v_enter_2 = get_int(2, "v_enter"),
    v_exit_2  = get_int(2, "v_exit"),
    v_enter_3 = get_int(3, "v_enter"),
    v_exit_3  = get_int(3, "v_exit")
  )
}))

# Merge: coarse summary + eps robustness at eps=0 + dense results
if (have_dt) {
  dense_dt <- data.table::as.data.table(dense_tbl)
  
  eps0 <- eps_dt[eps == 0, .(logS, delta, n_switches_eps_0 = n_switches_eps,
                             n_switches_astar_grid = n_switches_astar)]
  diag_dt <- merge(multi[, .(logS, delta, n_switches_a_star)], eps0, by = c("logS", "delta"), all.x = TRUE)
  diag_dt <- merge(diag_dt, dense_dt, by = c("logS", "delta"), all.x = TRUE)
  
  diag_outfile <- file.path(out_dir, "nonmonotone_cells_diagnostics.csv")
  data.table::fwrite(diag_dt, diag_outfile)
} else {
  # base R merge
  eps0 <- eps_dt[eps_dt$eps == 0, c("logS", "delta", "n_switches_eps", "n_switches_astar")]
  names(eps0) <- c("logS", "delta", "n_switches_eps_0", "n_switches_astar_grid")
  
  diag_df <- merge(multi[, c("logS", "delta", "n_switches_a_star")], eps0,
                   by = c("logS", "delta"), all.x = TRUE)
  diag_df <- merge(diag_df, dense_tbl, by = c("logS", "delta"), all.x = TRUE)
  
  diag_outfile <- file.path(out_dir, "nonmonotone_cells_diagnostics.csv")
  utils::write.csv(diag_df, diag_outfile, row.names = FALSE)
}

cat("\nWrote diagnostics CSV:\n  ", diag_outfile, "\n", sep = "")
cat("Wrote plots to:\n  ", fig_dir, "\n", sep = "")
cat("\nDone.\n")
