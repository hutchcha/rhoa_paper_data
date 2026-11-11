# Cell-assay quantification (Figure 5 / SI Figure S3)

Quantification and statistics for the NIH3T3 RhoA immunofluorescence experiments.

```bash
python Analyze_experiments.py      # writes figs_out/{serum,starve}_stacked.pdf + per-metric ANOVA/Tukey stats
```

Source data (per-cell, expression-normalized intensities; wide CSV, one column per
condition):
- `phalloidin_serum.csv`, `phalloidin_starve.csv` — F-actin (phalloidin).
- `pmlc2_serum.csv`, `pmlc2_starve.csv` — phospho-myosin light chain 2.

For each feeding condition the script builds a stacked panel per metric, runs one-way
ANOVA with Tukey HSD post-hoc, and writes the posthoc tables and ANOVA summaries.

> Note: a stray empty argument in the `plt.subplots(...)` call (a dangling comma that
> raised `SyntaxError`) was removed from the original script so it parses and runs; no
> analysis logic was changed.
