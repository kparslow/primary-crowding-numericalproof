################################################################################
################## Comparative Statics: Action Cutoffs ########################
################################################################################
# Date created: 4/13/2026, 8:42 AM
# Author: Katherine Parslow
#
# PURPOSE: Analyze how the action cutoff v*(delta,S) varies with parameters
#          within the single-tail regime (delta > delta* approx 0.10)
#
# From paper: "As S rises, the nomination advantage becomes more valuable, 
# expanding the set of types for whom a = 1 is optimal. As delta rises, the 
# electability cost becomes more severe, shrinking that set"
# [Page 7, Primary_Competition_(9).pdf]
#
# REQUIRES: Must run numerical_proof.R first to generate results/
################################################################################

cat("\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("COMPARATIVE STATICS: ACTION CUTOFFS v*(delta,S)\n")
cat("Analysis conditional on single-tail regime (delta > 0.10)\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

# ==============================================================================
# SETUP
# ==============================================================================

if (!file.exists("results/numerical_proof_results.RData")) {
  stop("Error: Run numerical_proof.R first to generate results/")
}

load("results/numerical_proof_results.RData")

# Load functions
G <- function(v, S) {
  numerator <- (exp(1) - 1) * S
  denominator <- exp(v + 1) + S
  return(1 + numerator / denominator)
}

C <- function(v, delta) {
  numerator <- pnorm(v)
  denominator <- pnorm(v - delta)
  return(numerator / denominator)
}

H <- function(v, S, delta) {
  return(G(v, S) - C(v, delta))
}

cat("Loaded results from numerical_proof.R\n")
cat(sprintf("  %d parameter combinations analyzed\n", nrow(results)))
cat(sprintf("  Testing range: v in [%d, %d]\n\n", v_min, v_max))

# ==============================================================================
# PART 1: REGIME IDENTIFICATION
# ==============================================================================

cat("Part 1: Identifying regime boundary...\n\n")

# From paper: "Moreover, in parameter ranges relevant for the analysis below, 
# a*(.; S) takes one of the following two forms: (i) Single cutoff... 
# (ii) Union..." [Page 9, Primary_Competition_(9).pdf, Proposition 5]

all_action_cases <- results[results$structure == "All action", ]
single_tail_cases <- results[results$structure == "Single upper tail", ]

cat(sprintf("All action regime: %d cases (%.1f%%)\n",
            nrow(all_action_cases), 100*nrow(all_action_cases)/nrow(results)))
cat(sprintf("Single upper tail regime: %d cases (%.1f%%)\n\n",
            nrow(single_tail_cases), 100*nrow(single_tail_cases)/nrow(results)))

# Find transition boundary
S_unique <- sort(unique(results$S))
transition_deltas <- data.frame(S = numeric(), delta_star = numeric())

for (S_val in S_unique) {
  S_data <- results[results$S == S_val, ]
  S_data <- S_data[order(S_data$delta), ]
  
  first_tail <- which(S_data$structure == "Single upper tail")[1]
  
  if (!is.na(first_tail) && first_tail > 1) {
    delta_star <- S_data$delta[first_tail]
    transition_deltas <- rbind(transition_deltas, 
                               data.frame(S = S_val, delta_star = delta_star))
  }
}

if (nrow(transition_deltas) > 0) {
  delta_star_mean <- mean(transition_deltas$delta_star)
  delta_star_sd <- sd(transition_deltas$delta_star)
  delta_star_min <- min(transition_deltas$delta_star)
  delta_star_max <- max(transition_deltas$delta_star)
  
  cat("Regime transition analysis:\n")
  cat(sprintf("  Transition delta*: mean = %.4f, SD = %.4f\n",
              delta_star_mean, delta_star_sd))
  cat(sprintf("  Range: [%.4f, %.4f]\n", delta_star_min, delta_star_max))
  
  # Check for variance before computing correlation
  if (delta_star_sd < 1e-10) {
    # delta* is essentially constant
    cat(sprintf("  delta* is constant across all S (SD < 0.0001)\n\n"))
    
    cat("KEY FINDING: delta* is S-invariant\n")
    cat(sprintf("  The regime boundary occurs at delta* approx %.2f for ALL S in [%d, %d]\n",
                delta_star_mean, min(transition_deltas$S), max(transition_deltas$S)))
    cat("  From paper: 'Numerical exploration and the monotone single-index\n")
    cat("  structure of G suggest that economically relevant parameter combinations\n")
    cat("  generate either a single upper tail or the union of an interior interval\n")
    cat("  and an upper tail' [Page 14, Primary_Competition_(9).pdf, Appendix]\n\n")
    cat("  Our finding: The regime transition is determined purely by delta,\n")
    cat("  not by competition intensity S. This reveals that:\n")
    cat("    - delta controls the EXTENSIVE margin (whether any candidates abstain)\n")
    cat("    - S controls the INTENSIVE margin (which candidates abstain)\n\n")
    
  } else {
    # delta* varies with S - compute correlation
    cor_test <- cor.test(transition_deltas$S, transition_deltas$delta_star)
    cat(sprintf("  Correlation with S: %.4f (p = %.4f)\n\n",
                cor_test$estimate, cor_test$p.value))
    
    if (cor_test$p.value > 0.05) {
      cat("delta* is S-invariant (p > 0.05)\n")
      cat("  The regime boundary is determined purely by penalty delta,\n")
      cat("  not by competition intensity S\n\n")
    } else {
      cat("WARNING: delta* varies with S (p < 0.05)\n")
      cat(sprintf("  Effect: correlation = %.3f\n\n", cor_test$estimate))
    }
  }
  
} else {
  stop("No regime transition found. Check parameter ranges.")
}

# ==============================================================================
# PART 2: EXTRACT CUTOFFS (SINGLE-TAIL REGIME ONLY)
# ==============================================================================

cat("Part 2: Extracting action cutoffs from single-tail regime...\n\n")

# From paper: "Lemma 3 shows that the electability penalty is most severe, 
# in proportional terms, for weaker candidates and becomes negligible for 
# very strong candidates" [Page 7, Primary_Competition_(9).pdf]

# Filter: only analyze cases above transition threshold
single_tail_filtered <- single_tail_cases[single_tail_cases$delta > delta_star_mean, ]

cat(sprintf("Analyzing %d single-tail cases (delta > %.2f)\n",
            nrow(single_tail_filtered), delta_star_mean))
cat("Computing precise cutoffs (this may take 1-2 minutes)...\n")

cutoffs <- data.frame(
  delta = numeric(),
  S = numeric(),
  v_star = numeric(),
  stringsAsFactors = FALSE
)

pb_start <- Sys.time()

for (i in 1:nrow(single_tail_filtered)) {
  delta_i <- single_tail_filtered$delta[i]
  S_i <- single_tail_filtered$S[i]
  
  v_seq <- seq(v_min, v_max, by = 0.01)
  H_vals <- sapply(v_seq, function(v) H(v, S_i, delta_i))
  
  sign_changes <- which(diff(sign(H_vals)) != 0)
  
  if (length(sign_changes) > 0) {
    idx <- sign_changes[1]
    v_before <- v_seq[idx]
    v_after <- v_seq[idx + 1]
    H_before <- H_vals[idx]
    H_after <- H_vals[idx + 1]
    
    v_star <- v_before - H_before * (v_after - v_before) / (H_after - H_before)
    
    cutoffs <- rbind(cutoffs, data.frame(
      delta = delta_i,
      S = S_i,
      v_star = v_star
    ))
  }
  
  if (i %% 100 == 0) {
    cat(sprintf("  Progress: %d/%d (%.0f%%)\r", 
                i, nrow(single_tail_filtered), 100*i/nrow(single_tail_filtered)))
  }
}

runtime <- as.numeric(difftime(Sys.time(), pb_start, units = "secs"))
cat(sprintf("\n\nExtracted %d cutoffs in %.1f seconds\n\n", 
            nrow(cutoffs), runtime))

# ==============================================================================
# PART 3: REGRESSION ANALYSIS
# ==============================================================================

cat("Part 3: Regression analysis of cutoff determinants...\n\n")

# From paper: "As S rises, the nomination advantage becomes more valuable, 
# expanding the set of types for whom a = 1 is optimal. As delta rises, the 
# electability cost becomes more severe, shrinking that set"
# [Page 7, Primary_Competition_(9).pdf]

cat(paste(rep("-", 78), collapse = ""), "\n")
cat("LINEAR MODEL: v* ~ delta + S\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")

lm_linear <- lm(v_star ~ delta + S, data = cutoffs)
summary_linear <- summary(lm_linear)

print(summary_linear)

beta_delta <- coef(lm_linear)["delta"]
beta_S <- coef(lm_linear)["S"]

cat("\n")
cat(paste(rep("-", 78), collapse = ""), "\n")
cat("INTERPRETATION OF MARGINAL EFFECTS\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")

cat("dv*/ddelta (holding S fixed):\n")
cat(sprintf("  Coefficient: %.4f\n", beta_delta))
cat(sprintf("  Interpretation: A 0.01 increase in delta raises cutoff by %.4f\n",
            beta_delta * 0.01))
cat("  From paper: 'As delta rises, the electability cost becomes more severe,\n")
cat("  shrinking that set' [Page 7, Primary_Competition_(9).pdf]\n")
cat(sprintf("  -> Higher penalties screen out weaker candidates\n\n"))

cat("dv*/dS (holding delta fixed):\n")
cat(sprintf("  Coefficient: %.6f\n", beta_S))
cat(sprintf("  Interpretation: A 10-unit increase in S changes cutoff by %.4f\n",
            beta_S * 10))
cat("  From paper: 'As S rises, the nomination advantage becomes more valuable,\n")
cat("  expanding the set of types for whom a = 1 is optimal'\n")
cat("  [Page 7, Primary_Competition_(9).pdf]\n")
cat(sprintf("  -> More crowding brings weaker candidates into action\n\n"))

# Elasticities
mean_delta <- mean(cutoffs$delta)
mean_S <- mean(cutoffs$S)
mean_vstar <- mean(cutoffs$v_star)

elasticity_delta <- beta_delta * (mean_delta / mean_vstar)
elasticity_S <- beta_S * (mean_S / mean_vstar)

cat(paste(rep("-", 78), collapse = ""), "\n")
cat("ELASTICITIES (at sample means)\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")

cat(sprintf("Elasticity(v*, delta) = %.4f\n", elasticity_delta))
cat(sprintf("  -> 1%% increase in delta -> %.3f%% increase in v*\n\n",
            elasticity_delta * 100))

cat(sprintf("Elasticity(v*, S) = %.4f\n", elasticity_S))
cat(sprintf("  -> 1%% increase in S -> %.3f%% change in v*\n\n",
            elasticity_S * 100))

# Standardized effects
sd_delta <- sd(cutoffs$delta)
sd_S <- sd(cutoffs$S)
sd_vstar <- sd(cutoffs$v_star)

std_effect_delta <- beta_delta * sd_delta / sd_vstar
std_effect_S <- beta_S * sd_S / sd_vstar

cat(paste(rep("-", 78), collapse = ""), "\n")
cat("STANDARDIZED EFFECTS (effect of 1 SD change)\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")

cat(sprintf("1 SD increase in delta (%.4f):\n", sd_delta))
cat(sprintf("  -> %.3f SD change in v*\n\n", std_effect_delta))

cat(sprintf("1 SD increase in S (%.2f):\n", sd_S))
cat(sprintf("  -> %.3f SD change in v*\n\n", std_effect_S))

# Compare magnitudes
if (abs(std_effect_delta) > abs(std_effect_S)) {
  cat("Penalty effect dominates: delta has larger standardized impact on v*\n")
  cat("  From paper: 'larger general-election penalties work in the opposite direction'\n")
  cat("  [Page 7, Primary_Competition_(9).pdf]\n\n")
} else {
  cat("Crowding effect dominates: S has larger standardized impact on v*\n")
  cat("  From paper: 'Crowding, captured by a larger S, increases the proportional\n")
  cat("  value of the primary boost' [Page 7, Primary_Competition_(9).pdf]\n\n")
}

# ==============================================================================
# PART 4: INTERACTION EFFECTS
# ==============================================================================

cat(paste(rep("-", 78), collapse = ""), "\n")
cat("TESTING FOR INTERACTION: v* ~ delta x S\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")

# From paper: "Because both G(v; S) and C(v; delta) decreasing in v, the inequality 
# G(v; S) >= C(v; delta) need not exhibit single-crossing in v"
# [Page 8, Primary_Competition_(9).pdf]

lm_interaction <- lm(v_star ~ delta * S, data = cutoffs)

cat("Model comparison (ANOVA):\n")
anova_result <- anova(lm_linear, lm_interaction)
print(anova_result)

cat("\n")

# Extract p-value safely using column name without backticks
interaction_pval <- anova_result[2, "Pr(>F)"]

if (!is.na(interaction_pval) && interaction_pval < 0.05) {
  cat("WARNING: SIGNIFICANT INTERACTION DETECTED (p < 0.05)\n\n")
  cat("The effect of delta on v* depends on S (or vice versa)\n")
  cat("From paper: 'Primary crowding affects incentives only through G(v; S),\n")
  cat("while the general-election penalty affects incentives only through C(v; delta)'\n")
  cat("[Page 7, Primary_Competition_(9).pdf]\n\n")
  
  interaction_coef <- coef(lm_interaction)["delta:S"]
  cat(sprintf("Interaction coefficient: %.6f\n", interaction_coef))
  
  if (interaction_coef > 0) {
    cat("Interpretation: The penalty effect (dv*/ddelta) is STRONGER at higher S\n")
    cat("-> In crowded primaries, penalties have larger screening effects\n\n")
  } else {
    cat("Interpretation: The penalty effect (dv*/ddelta) is WEAKER at higher S\n")
    cat("-> Crowding partially offsets the screening effect of penalties\n\n")
  }
  
} else {
  cat("NO SIGNIFICANT INTERACTION (p > 0.05)\n\n")
  cat("The effects of delta and S on v* are roughly additive\n")
  cat("This supports the decomposition in the paper:\n")
  cat("'Primary crowding affects incentives only through G(v; S), while the\n")
  cat("general-election penalty affects incentives only through C(v; delta)'\n")
  cat("[Page 7, Primary_Competition_(9).pdf]\n\n")
}

# ==============================================================================
# PART 5: VISUALIZATIONS
# ==============================================================================

cat(paste(rep("=", 78), collapse = ""), "\n")
cat("PART 5: Creating visualizations...\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

# Figure 1: v*(delta,S) vs delta for fixed S values
# From paper: "As delta rises, the electability cost becomes more severe, 
# shrinking that set" [Page 7, Primary_Competition_(9).pdf]

png("figures/cutoff_vs_delta.png", width = 1200, height = 800, res = 120)

par(mar = c(5, 5, 4, 8))

S_plot <- c(10, 30, 50, 100, 150, 200)
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")

plot(0, 0, type = "n", 
     xlim = range(cutoffs$delta), 
     ylim = range(cutoffs$v_star),
     xlab = "General election penalty (delta)",
     ylab = "Action cutoff v*(delta,S)",
     main = "Action Cutoff Increases with Penalty\n(Higher delta screens out weaker candidates)",
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

grid(col = "gray80")

# Add regression lines and data
for (i in 1:length(S_plot)) {
  S_val <- S_plot[i]
  data_S <- cutoffs[cutoffs$S == S_val, ]
  if (nrow(data_S) > 0) {
    data_S <- data_S[order(data_S$delta), ]
    lines(data_S$delta, data_S$v_star, col = colors[i], lwd = 2.5)
    points(data_S$delta, data_S$v_star, col = colors[i], pch = 19, cex = 0.8)
    
    # Add trend line
    if (nrow(data_S) > 2) {
      lm_S <- lm(v_star ~ delta, data = data_S)
      abline(lm_S, col = colors[i], lty = 2, lwd = 1.5)
    }
  }
}

# Add annotation
text(max(cutoffs$delta) * 0.3, max(cutoffs$v_star) * 0.85,
     sprintf("Avg effect: dv*/ddelta approx %.2f", beta_delta),
     cex = 1.3, font = 2, col = "#2c3e50")

legend("right", inset = c(-0.25, 0), xpd = TRUE,
       legend = paste0("S = ", S_plot),
       col = colors, lwd = 2.5, pch = 19,
       title = "Competition", bty = "n", cex = 1.1)

dev.off()

cat("  Created: figures/cutoff_vs_delta.png\n")

# Figure 2: v*(delta,S) vs S for fixed delta values
# From paper: "As S rises, the nomination advantage becomes more valuable, 
# expanding the set of types for whom a = 1 is optimal"
# [Page 7, Primary_Competition_(9).pdf]

png("figures/cutoff_vs_S.png", width = 1200, height = 800, res = 120)

par(mar = c(5, 5, 4, 8))

delta_plot <- c(0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
colors2 <- c("#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ff9896", "#9edae5")

plot(0, 0, type = "n",
     xlim = range(cutoffs$S),
     ylim = range(cutoffs$v_star),
     xlab = "Competition intensity S",
     ylab = "Action cutoff v*(delta,S)",
     main = "Action Cutoff Decreases with Competition\n(Crowding expands participation to weaker candidates)",
     cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

grid(col = "gray80")

for (i in 1:length(delta_plot)) {
  delta_val <- delta_plot[i]
  data_delta <- cutoffs[abs(cutoffs$delta - delta_val) < 0.001, ]
  if (nrow(data_delta) > 0) {
    data_delta <- data_delta[order(data_delta$S), ]
    lines(data_delta$S, data_delta$v_star, col = colors2[i], lwd = 2.5)
    points(data_delta$S, data_delta$v_star, col = colors2[i], pch = 19, cex = 0.8)
    
    if (nrow(data_delta) > 2) {
      lm_delta <- lm(v_star ~ S, data = data_delta)
      abline(lm_delta, col = colors2[i], lty = 2, lwd = 1.5)
    }
  }
}

# Add annotation
text(max(cutoffs$S) * 0.3, max(cutoffs$v_star) * 0.85,
     sprintf("Avg effect: dv*/dS approx %.4f", beta_S),
     cex = 1.3, font = 2, col = "#2c3e50")

legend("right", inset = c(-0.25, 0), xpd = TRUE,
       legend = paste0("delta = ", delta_plot),
       col = colors2, lwd = 2.5, pch = 19,
       title = "Penalty", bty = "n", cex = 1.1)

dev.off()

cat("  Created: figures/cutoff_vs_S.png\n")

# Figure 3: Heatmap
png("figures/cutoff_heatmap.png", width = 1000, height = 800, res = 120)

delta_unique <- sort(unique(cutoffs$delta))
S_unique <- sort(unique(cutoffs$S))

cutoff_matrix <- matrix(NA, nrow = length(delta_unique), ncol = length(S_unique))

for (i in 1:nrow(cutoffs)) {
  row_idx <- which(delta_unique == cutoffs$delta[i])
  col_idx <- which(S_unique == cutoffs$S[i])
  cutoff_matrix[row_idx, col_idx] <- cutoffs$v_star[i]
}

par(mar = c(5, 5, 4, 6))

image(delta_unique, S_unique, cutoff_matrix,
      col = heat.colors(50, rev = TRUE),
      xlab = "General election penalty (delta)",
      ylab = "Competition intensity S",
      main = "Action Cutoff v*(delta,S) Across Parameter Space",
      cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

contour(delta_unique, S_unique, cutoff_matrix, add = TRUE, 
        col = "black", lwd = 1.5)

dev.off()

cat("  Created: figures/cutoff_heatmap.png\n")

# ==============================================================================
# PART 6: SAVE RESULTS
# ==============================================================================

cat("\nPart 6: Saving comparative statics results...\n\n")

save(cutoffs, delta_star_mean, transition_deltas, 
     lm_linear, lm_interaction,
     file = "results/comparative_statics_results.RData")

write.csv(cutoffs, "results/action_cutoffs.csv", row.names = FALSE)

# Summary table for Appendix B
summary_table <- data.frame(
  Statistic = c(
    "Sample size",
    "delta range",
    "S range",
    "v* range",
    "Regime boundary (delta*)",
    "dv*/ddelta",
    "dv*/dS",
    "Elasticity(v*,delta)",
    "Elasticity(v*,S)",
    "R-squared",
    "Interaction p-value"
  ),
  Value = c(
    sprintf("%d", nrow(cutoffs)),
    sprintf("[%.2f, %.2f]", min(cutoffs$delta), max(cutoffs$delta)),
    sprintf("[%.0f, %.0f]", min(cutoffs$S), max(cutoffs$S)),
    sprintf("[%.3f, %.3f]", min(cutoffs$v_star), max(cutoffs$v_star)),
    sprintf("%.4f (SD: %.4f)", delta_star_mean, delta_star_sd),
    sprintf("%.4f", beta_delta),
    sprintf("%.6f", beta_S),
    sprintf("%.4f", elasticity_delta),
    sprintf("%.4f", elasticity_S),
    sprintf("%.4f", summary_linear$r.squared),
    sprintf("%.4f", interaction_pval)
  )
)

write.csv(summary_table, "results/comparative_statics_summary.csv", row.names = FALSE)

cat("Results saved:\n")
cat("  - results/comparative_statics_results.RData\n")
cat("  - results/action_cutoffs.csv\n")
cat("  - results/comparative_statics_summary.csv\n\n")

cat("Comparative statics analysis complete\n\n")

cat("SUMMARY FOR APPENDIX B:\n")
cat(paste(rep("-", 78), collapse = ""), "\n\n")
print(summary_table, row.names = FALSE)

cat("\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")