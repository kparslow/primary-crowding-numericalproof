# primary-crowding
Author: Katherine Parslow

Affiliation: Vanderbilt University, Department of Economics

## Overview ##
This repository houses on-going work on behavior distortion in two-stage electoral competition. The main artifact is the working paper titled "Crowded Primaries and Weakend Nominees: How Competition Distorts Candidate Behavior." The paper features a game-theoretic model of the trade-offs potential nominees face when campaign tactics boost their performance with their party base, but threaten their general-election viability. 

In addition to the formal analysis provided in the paper, this repo includes numerical analysis and simulations to demonstrate the model's key results.

## Model ##
- $N$ candidates compete in a primary election. Candidate $i \in \{1,...,N\}$ has valence $v_i \sim F$

## Repository Structure ##
```text
.
├── Code
│   ├── Figures.R
│   ├── comparative-statics.R
│   └── numerical_proof.R
├── README.md
├── docs
├── figures
│   ├── cutoff_heatmap.png
│   ├── cutoff_vs_S.png
│   ├── cutoff_vs_delta.png
│   ├── fig1_phase_diagram.png
│   ├── fig2_transition_boundary.png
│   └── fig3_H_curves_panel.png
├── primary-crowding-numerical-proof.Rproj
└── results
    ├── action_cutoffs.csv
    ├── comparative_statics_results.RData
    ├── comparative_statics_summary.csv
    ├── numerical_proof_results.RData
    ├── numerical_proof_table.csv
    └── summary_statistics.csv
```
