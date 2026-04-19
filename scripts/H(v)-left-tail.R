# ============================================================================
# REFINED GRID SEARCH FOR nu_min
# Focus on delta dependency (nu_min varies with delta, not S)
# ============================================================================

library(tidyverse)
library(ggplot2)
library(gridExtra)

# ============================================================================
# LOG-SPACE FUNCTIONS
# ============================================================================

log_G <- function(v, S) {
  log_p1 <- (v + 1) - log(exp(v + 1) + S)
  log_p0 <- v - log(exp(v) + S)
  return(log_p1 - log_p0)
}

log_C <- function(v, delta) {
  log_Phi_v <- pnorm(v, log.p = TRUE)
  log_Phi_v_minus_delta <- pnorm(v - delta, log.p = TRUE)
  return(log_Phi_v - log_Phi_v_minus_delta)
}

H_log <- function(v, S, delta) {
  return(log_G(v, S) - log_C(v, delta))
}

# ============================================================================
# SEARCH FUNCTION
# ============================================================================

find_nu_min_precise <- function(delta, S, v_range = c(-200, -5),
                                coarse_step = 10, fine_step = 0.1) {
  # Search for nu_min using three-stage refinement
  
  v_left <- v_range[1]
  v_right <- v_range[2]
  
  # Stage 1: Coarse search
  v_coarse <- seq(v_left, v_right, by = coarse_step)
  H_log_coarse <- rep(NA, length(v_coarse))
  
  for (i in seq_along(v_coarse)) {
    tryCatch({
      H_log_coarse[i] <- H_log(v_coarse[i], S, delta)
    }, error = function(e) {
      # Keep as NA
    })
  }
  
  # Find sign changes
  sign_coarse <- sign(H_log_coarse)
  crossings <- which(diff(sign_coarse) != 0)
  
  if (length(crossings) == 0) {
    return(list(nu_min = NA, found = FALSE))
  }
  
  # Use rightmost crossing
  crossing_idx <- crossings[length(crossings)]
  v_left_coarse <- v_coarse[crossing_idx]
  v_right_coarse <- v_coarse[crossing_idx + 1]
  
  # Stage 2: Medium search
  v_medium <- seq(v_left_coarse, v_right_coarse, by = coarse_step/5)
  H_log_medium <- rep(NA, length(v_medium))
  
  for (i in seq_along(v_medium)) {
    tryCatch({
      H_log_medium[i] <- H_log(v_medium[i], S, delta)
    }, error = function(e) {
      # Keep as NA
    })
  }
  
  sign_medium <- sign(H_log_medium)
  crossings_medium <- which(diff(sign_medium) != 0)
  
  if (length(crossings_medium) == 0) {
    return(list(nu_min = NA, found = FALSE))
  }
  
  crossing_idx_m <- crossings_medium[length(crossings_medium)]
  v_left_medium <- v_medium[crossing_idx_m]
  v_right_medium <- v_medium[crossing_idx_m + 1]
  
  # Stage 3: Fine search
  v_fine <- seq(v_left_medium, v_right_medium, by = fine_step)
  H_log_fine <- rep(NA, length(v_fine))
  
  for (i in seq_along(v_fine)) {
    tryCatch({
      H_log_fine[i] <- H_log(v_fine[i], S, delta)
    }, error = function(e) {
      # Keep as NA
    })
  }
  
  sign_fine <- sign(H_log_fine)
  crossings_fine <- which(diff(sign_fine) != 0)
  
  if (length(crossings_fine) == 0) {
    nu_min <- v_fine[which.min(abs(H_log_fine))]
    return(list(nu_min = nu_min, found = TRUE))
  }
  
  # Linear interpolation
  idx <- crossings_fine[length(crossings_fine)]
  v1 <- v_fine[idx]
  v2 <- v_fine[idx + 1]
  H1 <- H_log_fine[idx]
  H2 <- H_log_fine[idx + 1]
  
  nu_min <- v1 - H1 * (v2 - v1) / (H2 - H1)
  
  return(list(nu_min = nu_min, found = TRUE))
}

# ============================================================================
# MAIN: REFINED GRID SEARCH (FINER delta, COARSER S)
# ============================================================================

cat(strrep("=", 70), "\n", sep = "")
cat("REFINED GRID SEARCH FOR nu_min\n")
cat("Focus: Fine delta grid [0.01, 0.25], Coarse S grid\n")
cat(strrep("=", 70), "\n\n")

# Refined parameter grid
delta_grid <- seq(0.01, 0.25, by = 0.01)  # 25 fine values
S_grid <- seq(2, 200, by = 10)             # 20 coarse values
# Total: 500 combinations (much faster)

cat(sprintf("Computing nu_min for %d x %d = %d combinations...\n",
            length(delta_grid), length(S_grid), 
            length(delta_grid) * length(S_grid)))
cat("This will take 2-3 minutes. Progress:\n\n")

grid_results <- tibble()

for (i in seq_along(delta_grid)) {
  delta <- delta_grid[i]
  
  # Progress indicator
  if (i %% 5 == 0) {
    cat(sprintf("  delta = %.2f (%d/%d completed)\n", delta, i, length(delta_grid)))
  }
  
  for (S in S_grid) {
    result <- find_nu_min_precise(
      delta = delta,
      S = S,
      v_range = c(-200, -5),
      coarse_step = 10,
      fine_step = 0.1
    )
    
    if (result$found) {
      nu_min <- result$nu_min
    } else {
      nu_min <- NA
    }
    
    grid_results <- bind_rows(grid_results,
                              tibble(delta = delta, S = S, nu_min = nu_min))
  }
}

cat("\n✓ Grid search complete!\n\n")

# ============================================================================
# ANALYSIS OF RESULTS
# ============================================================================

cat(strrep("=", 70), "\n", sep = "")
cat("ANALYSIS: HOW nu_min VARIES WITH delta\n")
cat(strrep("=", 70), "\n\n")

nu_min_valid <- grid_results %>% filter(!is.na(nu_min))

cat(sprintf("Total parameter combinations: %d\n", nrow(grid_results)))
cat(sprintf("Combinations with nu_min found: %d\n", nrow(nu_min_valid)))
cat(sprintf("Combinations with NA: %d\n\n", nrow(grid_results) - nrow(nu_min_valid)))

cat("OVERALL nu_min STATISTICS:\n")
cat(sprintf("  Minimum: %.3f\n", min(nu_min_valid$nu_min, na.rm = TRUE)))
cat(sprintf("  25th percentile: %.3f\n", quantile(nu_min_valid$nu_min, 0.25)))
cat(sprintf("  Median: %.3f\n", median(nu_min_valid$nu_min, na.rm = TRUE)))
cat(sprintf("  75th percentile: %.3f\n", quantile(nu_min_valid$nu_min, 0.75)))
cat(sprintf("  Maximum: %.3f\n", max(nu_min_valid$nu_min, na.rm = TRUE)))
cat(sprintf("  Mean: %.3f\n", mean(nu_min_valid$nu_min, na.rm = TRUE)))
cat(sprintf("  Std Dev: %.3f\n\n", sd(nu_min_valid$nu_min, na.rm = TRUE)))

# KEY FINDING: How many are > -10?
n_greater_than_minus_10 <- sum(nu_min_valid$nu_min > -10)
pct_greater_than_minus_10 <- 100 * n_greater_than_minus_10 / nrow(nu_min_valid)

cat("KEY FINDING:\n")
cat(sprintf("  nu_min > -10 (within economically relevant range): %d/%d (%.1f%%)\n",
            n_greater_than_minus_10, nrow(nu_min_valid), pct_greater_than_minus_10))
cat(sprintf("  nu_min <= -10 (outside economically relevant range): %d/%d (%.1f%%)\n\n",
            nrow(nu_min_valid) - n_greater_than_minus_10, nrow(nu_min_valid),
            100 - pct_greater_than_minus_10))

# ============================================================================
# HOW nu_min VARIES WITH delta (PRIMARY FOCUS)
# ============================================================================

cat(strrep("=", 70), "\n", sep = "")
cat("HOW nu_min VARIES WITH delta (PRIMARY FINDING)\n")
cat(strrep("=", 70), "\n\n")

by_delta <- nu_min_valid %>%
  group_by(delta) %>%
  summarize(
    n_cases = n(),
    nu_min_min = min(nu_min),
    nu_min_max = max(nu_min),
    nu_min_mean = mean(nu_min),
    nu_min_sd = sd(nu_min),
    n_greater_than_10 = sum(nu_min > -10),
    .groups = 'drop'
  )

cat(sprintf("%-8s %-10s %-12s %-12s %-12s %-12s %-10s\n",
            "delta", "N_S", "Min", "Max", "Mean", "Std Dev", "N>-10"))
cat(strrep("-", 80), "\n")

for (i in 1:nrow(by_delta)) {
  row <- by_delta[i, ]
  cat(sprintf("%.2f     %-10d %-12.3f %-12.3f %-12.3f %-12.3f %-10d\n",
              row$delta, row$n_cases, row$nu_min_min, row$nu_min_max,
              row$nu_min_mean, row$nu_min_sd, row$n_greater_than_10))
}

# ============================================================================
# HOW nu_min VARIES WITH S (SECONDARY CHECK)
# ============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("HOW nu_min VARIES WITH S (SECONDARY CHECK)\n")
cat(strrep("=", 70), "\n\n")

by_S <- nu_min_valid %>%
  group_by(S) %>%
  summarize(
    n_cases = n(),
    nu_min_min = min(nu_min),
    nu_min_max = max(nu_min),
    nu_min_mean = mean(nu_min),
    nu_min_sd = sd(nu_min),
    .groups = 'drop'
  )

cat(sprintf("%-8s %-10s %-12s %-12s %-12s %-12s\n",
            "S", "N_delta", "Min", "Max", "Mean", "Std Dev"))
cat(strrep("-", 70), "\n")

for (i in 1:nrow(by_S)) {
  row <- by_S[i, ]
  cat(sprintf("%-8d %-10d %-12.3f %-12.3f %-12.3f %-12.3f\n",
              as.integer(row$S), row$n_cases, row$nu_min_min, row$nu_min_max,
              row$nu_min_mean, row$nu_min_sd))
}

# ============================================================================
# IDENTIFY PROBLEMATIC REGIONS
# ============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("PARAMETER COMBINATIONS WHERE nu_min > -10\n")
cat(strrep("=", 70), "\n\n")

problematic <- nu_min_valid %>%
  filter(nu_min > -10) %>%
  arrange(desc(nu_min))

if (nrow(problematic) > 0) {
  cat(sprintf("Found %d problematic cases (nu_min > -10):\n\n", nrow(problematic)))
  
  cat(sprintf("%-8s %-8s %-10s\n", "delta", "S", "nu_min"))
  cat(strrep("-", 30), "\n")
  
  for (i in 1:min(30, nrow(problematic))) {
    row <- problematic[i, ]
    cat(sprintf("%.2f     %-8d %-10.3f\n", row$delta, row$S, row$nu_min))
  }
  
  if (nrow(problematic) > 30) {
    cat(sprintf("... and %d more\n", nrow(problematic) - 30))
  }
} else {
  cat("✓ No problematic cases found - all nu_min <= -10\n")
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

write_csv(grid_results, "nu_min_refined_grid.csv")
cat("\n\n✓ Full results saved to nu_min_refined_grid.csv\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

cat("Generating visualizations...\n\n")

# Plot 1: nu_min vs delta (primary relationship)
p1 <- by_delta %>%
  ggplot(aes(x = delta)) +
  geom_ribbon(aes(ymin = nu_min_min, ymax = nu_min_max), 
              fill = "lightblue", alpha = 0.5, label = "Range across S") +
  geom_line(aes(y = nu_min_mean), color = "darkblue", size = 1.2, label = "Mean") +
  geom_hline(yintercept = -10, linetype = "dashed", color = "red", size = 1.2,
             label = "Economics boundary (v=-10)") +
  labs(
    title = "PRIMARY: nu_min varies with delta (penalty)",
    x = "General-election penalty (delta)",
    y = "nu_min",
    subtitle = "Range shows min/max across S values; mean shown as line"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Plot 2: nu_min vs S (secondary check)
p2 <- by_S %>%
  ggplot(aes(x = S)) +
  geom_ribbon(aes(ymin = nu_min_min, ymax = nu_min_max), 
              fill = "lightgreen", alpha = 0.5, label = "Range across delta") +
  geom_line(aes(y = nu_min_mean), color = "darkgreen", size = 1.2, label = "Mean") +
  geom_hline(yintercept = -10, linetype = "dashed", color = "red", size = 1.2) +
  labs(
    title = "SECONDARY: nu_min varies minimally with S (competition)",
    x = "Competition intensity (S)",
    y = "nu_min",
    subtitle = "Range shows min/max across delta values; mean shown as line"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Plot 3: Heatmap
p3 <- nu_min_valid %>%
  ggplot(aes(x = delta, y = S, fill = nu_min)) +
  geom_raster() +
  geom_contour(aes(z = nu_min), color = "black", alpha = 0.2, size = 0.3) +
  scale_fill_gradient(low = "red", high = "blue", name = "nu_min") +
  labs(
    title = "nu_min landscape (delta vs S)",
    x = "General-election penalty (delta)",
    y = "Competition intensity (S)",
    fill = "nu_min"
  ) +
  theme_minimal()

# Plot 4: Distribution
p4 <- nu_min_valid %>%
  ggplot(aes(x = nu_min)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = -10, linetype = "dashed", color = "red", size = 1.2,
             label = "Boundary at -10") +
  geom_vline(xintercept = mean(nu_min_valid$nu_min, na.rm = TRUE),
             linetype = "solid", color = "green", size = 1.2,
             label = sprintf("Mean: %.2f", mean(nu_min_valid$nu_min, na.rm = TRUE))) +
  labs(
    title = "Distribution of nu_min",
    x = "nu_min",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Combine
combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
ggsave("nu_min_refined_grid_analysis.png", combined, width = 14, height = 10, dpi = 300)
cat("✓ Figure saved to nu_min_refined_grid_analysis.png\n\n")

# ============================================================================
# CONCLUSION
# ============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("CONCLUSION\n")
cat(strrep("=", 70), "\n\n")

# Calculate variation in nu_min due to S vs delta
sd_within_delta <- mean(by_delta$nu_min_sd, na.rm = TRUE)
sd_across_delta <- sd(by_delta$nu_min_mean, na.rm = TRUE)

cat("KEY INSIGHT:\n")
cat(sprintf("  Standard deviation WITHIN delta (across S): %.4f\n", sd_within_delta))
cat(sprintf("  Standard deviation ACROSS delta (mean values): %.4f\n", sd_across_delta))
cat(sprintf("  Ratio (across/within): %.1f\n\n", sd_across_delta / sd_within_delta))

cat("INTERPRETATION:\n")
if (sd_across_delta > 5 * sd_within_delta) {
  cat("  ✓ nu_min varies PRIMARILY with delta, minimally with S\n")
  cat("  ✓ Paper's claim supported: nu_min depends on penalty, not crowding\n")
} else if (sd_across_delta > 2 * sd_within_delta) {
  cat("  ≈ nu_min varies MORE with delta than S, but S has non-negligible effect\n")
} else {
  cat("  ⚠ nu_min varies similarly with both delta and S\n")
}

cat("\n")

if (pct_greater_than_minus_10 == 0) {
  cat("✓ ALL nu_min < -10: Paper's bound is VALID\n")
} else if (pct_greater_than_minus_10 < 5) {
  cat(sprintf("⚠ Few exceptions (%.1f%% have nu_min > -10): Minor caveats needed\n", 
              pct_greater_than_minus_10))
} else {
  cat(sprintf("✗ SIGNIFICANT violations (%.1f%% have nu_min > -10): Contradicts paper\n", 
              pct_greater_than_minus_10))
}

cat("\n")