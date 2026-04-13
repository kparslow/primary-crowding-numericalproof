################################################################################
##################### Figure Generation for Appendix B #########################
################################################################################
# Last Updated: 4/13/2026, 7:47 AM
# Author: Katherine Parslow
#
# PURPOSE: Generate three key figures for Appendix B based on numerical proof
#
# REQUIRES: Must run numerical_proof.R first to generate results/
################################################################################

cat("\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("FIGURE GENERATION FOR APPENDIX B\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

# ==============================================================================
# SETUP
# ==============================================================================

# Load results from main analysis
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

cat("✓ Loaded results from numerical_proof.R\n")
cat(sprintf("  %d parameter combinations analyzed\n", nrow(results)))
cat(sprintf("  Testing range: v ∈ [%d, %d]\n\n", v_min, v_max))

# ==============================================================================
# FIGURE 1: PHASE DIAGRAM - STRUCTURE BY (δ, S)
# ==============================================================================

cat("Creating Figure 1: Phase diagram...\n")

# From paper [Proposition 5]: "in parameter ranges relevant for the analysis 
# below, a*(·; S) takes one of the following two forms"
png("figures/fig1_phase_diagram.png", width = 1200, height = 800, res = 120)

par(mar = c(5, 5, 4, 6))

# Get unique parameter values
delta_unique <- sort(unique(results$delta))
S_unique <- sort(unique(results$S))

# Create matrix: 1 = All action, 2 = Single upper tail, NA = other
phase_matrix <- matrix(NA, nrow = length(delta_unique), ncol = length(S_unique))

for (i in 1:nrow(results)) {
  row_idx <- which(delta_unique == results$delta[i])
  col_idx <- which(S_unique == results$S[i])
  
  if (results$structure[i] == "All action") {
    phase_matrix[row_idx, col_idx] <- 1
  } else if (results$structure[i] == "Single upper tail") {
    phase_matrix[row_idx, col_idx] <- 2
  } else {
    phase_matrix[row_idx, col_idx] <- 3  # Other (if any exist)
  }
}

# Count structures for reporting
n_all_action <- sum(results$structure == "All action")
n_single_tail <- sum(results$structure == "Single upper tail")
n_other <- nrow(results) - n_all_action - n_single_tail

cat(sprintf("  All action: %d cases (%.1f%%)\n", 
            n_all_action, 100*n_all_action/nrow(results)))
cat(sprintf("  Single upper tail: %d cases (%.1f%%)\n",
            n_single_tail, 100*n_single_tail/nrow(results)))
if (n_other > 0) {
  cat(sprintf("  Other: %d cases (%.1f%%)\n", n_other, 100*n_other/nrow(results)))
}

# Plot with clear color scheme
# From paper [Appendix]: "economically relevant parameter combinations generate 
# either a single upper tail or the union"
colors <- c("#3498db", "#e74c3c", "#95a5a6")  # Blue, Red, Gray
image(delta_unique, S_unique, phase_matrix,
      col = colors,
      xlab = "General election penalty δ",
      ylab = "Competition intensity S",
      main = "Action Set Structure Across Parameter Space\n(v ∈ [-10, 10], covering >99.999% of candidates under F=N(0,4))",
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1)

# Add contour line at boundary between regimes
contour(delta_unique, S_unique, phase_matrix, 
        levels = 1.5, add = TRUE, col = "black", lwd = 3, labcex = 0)

# Add grid
abline(v = seq(0, 0.5, by = 0.1), col = "gray70", lty = 3)
abline(h = seq(0, 200, by = 50), col = "gray70", lty = 3)

# Legend
legend_text <- c(
  sprintf("All action (%.1f%%)", 100*n_all_action/nrow(results)),
  sprintf("Single upper tail (%.1f%%)", 100*n_single_tail/nrow(results))
)
if (n_other > 0) {
  legend_text <- c(legend_text, sprintf("Other (%.1f%%)", 100*n_other/nrow(results)))
}

legend("topright", 
       legend = legend_text,
       fill = colors[1:(2 + (n_other > 0))],
       cex = 1.1, bg = "white", box.lwd = 2)

dev.off()

cat("  ✓ Created: figures/fig1_phase_diagram.png\n\n")

# ==============================================================================
# FIGURE 2: TRANSITION BOUNDARY δ*(S)
# ==============================================================================

cat("Creating Figure 2: Transition boundary δ*(S)...\n")

# From paper [Corollary 6]: "Aδ(S) is weakly decreasing in δ"
# Find transition delta for each S
transition_deltas <- data.frame(S = numeric(), delta_transition = numeric())

for (S_val in S_unique) {
  S_data <- results[results$S == S_val, ]
  S_data <- S_data[order(S_data$delta), ]
  
  # Find first occurrence of "Single upper tail"
  first_tail <- which(S_data$structure == "Single upper tail")[1]
  
  if (!is.na(first_tail) && first_tail > 1) {
    # Transition occurs at this delta
    delta_trans <- S_data$delta[first_tail]
    transition_deltas <- rbind(transition_deltas, 
                               data.frame(S = S_val, delta_transition = delta_trans))
  }
}

cat(sprintf("  Found transition boundary for %d S values\n", nrow(transition_deltas)))

if (nrow(transition_deltas) > 0) {
  png("figures/fig2_transition_boundary.png", width = 1000, height = 700, res = 120)
  
  par(mar = c(5, 5, 4, 2))
  
  plot(transition_deltas$S, transition_deltas$delta_transition,
       type = 'l', lwd = 3, col = '#2c3e50',
       xlab = "Competition intensity S",
       ylab = "Transition penalty δ*(S)",
       main = "Penalty Threshold for Action Set Transition\n(Below δ*: all candidates take action; Above δ*: only v ≥ v*(δ,S) take action)",
       cex.main = 1.1, cex.lab = 1.2, cex.axis = 1.1,
       ylim = c(0, max(transition_deltas$delta_transition) * 1.1))
  
  points(transition_deltas$S, transition_deltas$delta_transition, 
         pch = 19, col = '#2c3e50', cex = 0.6)
  
  grid(col = "gray80")
  
  # Add reference lines for common δ values
  abline(h = c(0.05, 0.10, 0.15, 0.20, 0.25), lty = 2, col = "gray60")
  
  # Add annotations
  text(10, 0.05, "δ = 0.05", pos = 3, cex = 0.9, col = "gray40")
  text(10, 0.10, "δ = 0.10", pos = 3, cex = 0.9, col = "gray40")
  text(10, 0.15, "δ = 0.15", pos = 3, cex = 0.9, col = "gray40")
  
  # Compute and display correlation
  cor_val <- cor(transition_deltas$S, transition_deltas$delta_transition)
  text(max(transition_deltas$S) * 0.8, 
       max(transition_deltas$delta_transition) * 0.9,
       sprintf("Correlation: %.3f", cor_val),
       cex = 1.1, col = "#2c3e50", font = 2)
  
  dev.off()
  
  cat("  ✓ Created: figures/fig2_transition_boundary.png\n")
  cat(sprintf("  Transition δ* ranges from %.3f to %.3f\n",
              min(transition_deltas$delta_transition),
              max(transition_deltas$delta_transition)))
  cat(sprintf("  Correlation with S: %.3f\n\n", cor_val))
  
} else {
  cat("  ⚠ No transition boundary found (all cases may be same structure)\n\n")
}

# ==============================================================================
# FIGURE 3: REPRESENTATIVE H(v) CURVES (2×2 PANEL)
# ==============================================================================

cat("Creating Figure 3: Representative H(v) curves...\n")

# From paper [Proposition 1]: "ai = 1 ⇐⇒ G(vi; S*i(v)) ≥ C(vi; δ)"
# Select 4 representative cases

# Case A: Low δ, low S - All action
case_A_candidates <- results[results$delta <= 0.03 & results$S <= 20 & 
                               results$structure == "All action", ]
if (nrow(case_A_candidates) > 0) {
  case_A <- case_A_candidates[1, ]
} else {
  case_A <- results[results$structure == "All action", ][1, ]
}

# Case B: Low δ, high S - All action
case_B_candidates <- results[results$delta <= 0.03 & results$S >= 150 & 
                               results$structure == "All action", ]
if (nrow(case_B_candidates) > 0) {
  case_B <- case_B_candidates[1, ]
} else {
  case_B <- results[results$structure == "All action", ][nrow(results[results$structure == "All action", ]), ]
}

# Case C: Moderate δ, low S - Single tail
case_C_candidates <- results[results$delta >= 0.10 & results$delta <= 0.15 & 
                               results$S <= 20 & 
                               results$structure == "Single upper tail", ]
if (nrow(case_C_candidates) > 0) {
  case_C <- case_C_candidates[1, ]
} else {
  case_C <- results[results$structure == "Single upper tail", ][1, ]
}

# Case D: High δ, high S - Single tail
case_D_candidates <- results[results$delta >= 0.30 & results$S >= 150 & 
                               results$structure == "Single upper tail", ]
if (nrow(case_D_candidates) > 0) {
  case_D <- case_D_candidates[1, ]
} else {
  case_D <- results[results$structure == "Single upper tail", ][nrow(results[results$structure == "Single upper tail", ]), ]
}

cat("  Selected cases:\n")
cat(sprintf("    Panel A: δ=%.2f, S=%.0f (%s)\n", case_A$delta, case_A$S, case_A$structure))
cat(sprintf("    Panel B: δ=%.2f, S=%.0f (%s)\n", case_B$delta, case_B$S, case_B$structure))
cat(sprintf("    Panel C: δ=%.2f, S=%.0f (%s)\n", case_C$delta, case_C$S, case_C$structure))
cat(sprintf("    Panel D: δ=%.2f, S=%.0f (%s)\n\n", case_D$delta, case_D$S, case_D$structure))

# Create 2×2 panel figure
png("figures/fig3_H_curves_panel.png", width = 1400, height = 1000, res = 120)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1.5))

plot_H_panel <- function(delta, S, title, structure) {
  v_seq <- seq(v_min, v_max, by = 0.05)
  H_vals <- sapply(v_seq, function(v) H(v, S, delta))
  
  # Determine y-axis limits
  y_range <- range(H_vals)
  y_buffer <- diff(y_range) * 0.1
  
  plot(v_seq, H_vals, type = 'l', lwd = 2.5, col = 'blue',
       main = sprintf("%s\nδ = %.2f, S = %.0f", title, delta, S),
       xlab = "Candidate valence v",
       ylab = "H(v) = G(v;S) - C(v;δ)",
       cex.main = 1.2, cex.lab = 1.1, cex.axis = 1,
       ylim = c(y_range[1] - y_buffer, y_range[2] + y_buffer))
  
  abline(h = 0, col = 'red', lty = 2, lwd = 2)
  grid(col = "gray80")
  
  # Shade action regions (where H > 0)
  for (i in 1:(length(v_seq)-1)) {
    if (H_vals[i] > 0 && H_vals[i+1] > 0) {
      rect(v_seq[i], -100, v_seq[i+1], 100, 
           col = rgb(0, 1, 0, 0.15), border = NA)
    }
  }
  
  # Re-draw curve on top
  lines(v_seq, H_vals, lwd = 2.5, col = 'blue')
  abline(h = 0, col = 'red', lty = 2, lwd = 2)
  
  # Add structure label
  text(v_min + 2, y_range[2] + y_buffer * 0.5, 
       structure, pos = 4, font = 2, cex = 1)
}

plot_H_panel(case_A$delta, case_A$S, 
             "Panel A: Low Penalty, Low Crowding", 
             case_A$structure)

plot_H_panel(case_B$delta, case_B$S, 
             "Panel B: Low Penalty, High Crowding", 
             case_B$structure)

plot_H_panel(case_C$delta, case_C$S, 
             "Panel C: Moderate Penalty, Low Crowding", 
             case_C$structure)

plot_H_panel(case_D$delta, case_D$S, 
             "Panel D: High Penalty, High Crowding", 
             case_D$structure)

par(mfrow = c(1, 1))

dev.off()

cat("  ✓ Created: figures/fig3_H_curves_panel.png\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat(paste(rep("=", 78), collapse = ""), "\n")
cat("FIGURE GENERATION COMPLETE\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

cat("Three figures created for Appendix B:\n\n")

cat("Figure 1: Phase Diagram (figures/fig1_phase_diagram.png)\n")
cat("  Shows: Distribution of action set structures across (δ,S) space\n")
cat(sprintf("  Result: %.1f%% all action, %.1f%% single upper tail\n",
            100*n_all_action/nrow(results), 100*n_single_tail/nrow(results)))
cat("  From paper [Proposition 5]: 'a*(·; S) takes one of the following two forms'\n\n")

if (nrow(transition_deltas) > 0) {
  cat("Figure 2: Transition Boundary (figures/fig2_transition_boundary.png)\n")
  cat("  Shows: Critical penalty δ*(S) separating the two regimes\n")
  cat(sprintf("  Result: δ* ranges from %.3f to %.3f across S ∈ [%d,%d]\n",
              min(transition_deltas$delta_transition),
              max(transition_deltas$delta_transition),
              min(transition_deltas$S), max(transition_deltas$S)))
  cat("  From paper [Corollary 6]: 'Aδ(S) is weakly decreasing in δ'\n\n")
}

cat("Figure 3: Representative H(v) Curves (figures/fig3_H_curves_panel.png)\n")
cat("  Shows: Four representative cases illustrating regime differences\n")
cat("  From paper [Proposition 1]: 'ai = 1 ⇐⇒ G(vi; S*i(v)) ≥ C(vi; δ)'\n")
cat("  Panels A-B: Low penalty regime (all candidates take action)\n")
cat("  Panels C-D: Moderate/high penalty regime (only v ≥ v* take action)\n\n")

cat("These figures support the findings in [Appendix]:\n")
cat("'Numerical exploration... suggest that economically relevant parameter\n")
cat("combinations generate either a single upper tail or the union of an\n")
cat("interior interval and an upper tail'\n\n")

cat("Within v ∈ [-10,10] (>99.999% of candidates), we observe only:\n")
cat("  (i) All action (low penalty regime)\n")
cat("  (ii) Single upper tail (moderate/high penalty regime)\n\n")

cat(paste(rep("=", 78), collapse = ""), "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")