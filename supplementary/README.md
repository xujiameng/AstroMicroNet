# Supplementary numerical materials

Run `python supplementary/plot_all.py` from the repository root after installing this folder's requirements. Outputs are saved in `figures/`. S1–S8 retain the V13 source tables and drawings. The repository-specific S9 output contains quantitative B–F only; microscopy panel A is explicitly omitted.

| Figures | Source |
|---|---|
| S1, S8 | `data/method_*.csv`; 20 networks per condition; paired methods |
| S2 | `data/rho_membership.csv`; 13 cells, mean ± SEM |
| S3 | `data/lag_membership.csv`, `data/waveform_*.csv` |
| S4 | `data/heldout_*.csv`; cell values averaged across three held-out folds |
| S5 | `data/average40_network.npz`, `data/generated5_network.npz`; 241-ROI real-data-derived matrices |
| S6 | `data/threshold_membership.csv` |
| S7 | `data/cell_metrics.csv`, `data/cell_identifiers.csv` |
| S9 B–F | `data/H1R_*.csv`; paired slice summaries and statistical results |

`supporting_data/` retains alternative BE diagnostics, all tested time-segment schemes, waveform calibration, null-network inputs and outputs, source-mouse comparisons and sampling metadata. The `generated5` filename refers to sampling of a real recording; it is not synthetic imaging.

Plotting reads relative paths and saved values. It does not rerun detection. No parameter selection was changed for packaging. The source manuscript reports historical core counts of 18–49 while the supplementary archive reports 18–59; these version-specific results were not combined or silently reconciled. ROI pairs or cells in the same slice must not be counted as independent animals.
