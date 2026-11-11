# RhoA paper — figure style guide

Canonical formatting for every figure in the RhoA / membrane-MD paper (target journal: **ACS**, e.g.
JACS). The `plot-polish` agent checks figures against this file and edits the plotting SCRIPT to
conform. Source of truth for the house style: the paper's plotting module
`…/Papers/Orientation Dynamics/RhoA/strided_trajs/hbonds_gr_analysis.py` (`OI_PALETTE`, `style_axes`).

Scope: apply **ad hoc** to figures actually going into the paper/SI — not a bulk pass over the whole
repo.

## 1. Global rcParams / style block (paste into every plotting script)
```python
import matplotlib, matplotlib.pyplot as plt
plt.rcParams["font.family"]   = "Arial"
plt.rcParams["axes.linewidth"] = 2
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["svg.fonttype"] = "none"   # editable text in Illustrator

def style_axes(ax):
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.yaxis.tick_left(); ax.xaxis.tick_bottom(); ax.minorticks_off()
    ax.spines["left"].set_linewidth(2); ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(width=1.5, direction="out")
```
Always save BOTH `.png` (dpi ≥ 300, `bbox_inches="tight"`) and `.svg`.

## 2. Color scheme (match the main figures — never improvise a color)
`OI_PALETTE = ["#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B"]` (blue, orange, green, red,
purple, brown), assigned to lipids **in this order**:

| Lipid | Hex | Color |
|---|---|---|
| POPC  | `#1F77B4` | blue |
| POPE  | `#FF7F0E` | orange |
| POPS  | `#2CA02C` | green |
| PLA18 | `#D62728` | red |
| PSM   | `#9467BD` | purple |

**Structural-region colors** (from main Fig 6/7 — use these wherever regions are colored):
N-patch `#8C564B` (brown), SwI `#C9A227`, interSw `#57A79B`, SwII `#8172B3`, HVR `#4C72B0`.

**Per-system colors** (locked — used wherever a plot is colored by simulation identity; these match the
PyMOL `SIM_COL` in `render_effector_snapshots.py`, so 2D plots and 3D structure renders agree):

| System | Hex |
|---|---|
| RhoA-GDP | `#B279A2` (mauve) |
| RhoA-GTP-1 | `#72B7B2` (teal) |
| RhoA-GTP-2 | `#E45756` (red) |
| RhoA-GTP-2A | `#F58518` (orange) |
| RhoA-GTP-3A | `#4C78A8` (blue) |

Note (accepted): GTP-2's red and GTP-3A's blue coincide with the complex/simple-membrane category
colors, but system-colored and membrane-colored plots are always different figures, so there is no
within-figure clash. Do NOT recolor to avoid this. EXCEPTION: when a figure colors systems by
**membrane type** instead of identity (e.g. S1 equilibration), use the membrane colors (simple
`#1F77B4`, complex `#D62728`), not this per-system table.

**SI category colors** (these categories appear only in SI figures; locked):

| Category | Color |
|---|---|
| simple membrane | `#1F77B4` (blue) |
| complex membrane | `#D62728` (red) |
| simulation (this work) | colored by the system's membrane (e.g. RhoA-GTP-2 → complex red) |
| experimental / reference value | `#7F7F7F` (gray) |
| insertion depth: from bilayer centre | `#9467BD` (purple) |
| insertion depth: from phosphate plane | `#FF7F0E` (orange) |

Rule: any lipid/region/category that appears in more than one figure MUST use the same color
everywhere. Reusing a palette color for a non-lipid category is OK **only** in a panel that shows no
lipids (e.g. simple=blue in a thickness panel, where POPC-blue is not also present).

## 3. Labels & titles
- **No titles on plots.** Describe the figure in the legend/caption. Multi-panel figures use bare bold
  panel letters (A, B, C…) only.
- **Axis labels in Title Case** (capitalize the principal words): "Bilayer Thickness (Å)",
  "Area Per Lipid (Å²)", "Anchor Insertion Depth (Å)", "H-Bond Occupancy", "Time (µs)".
  **Exception:** keep mathematical / variable notation verbatim — `g(r)`, `r (Å)`, `θ`, `φ`, `r_MI`.
- Units in parentheses; use Å, ², µ, θ, φ, Δ directly (Arial has them). For sub/superscripts in text
  use mathtext (`$\sum_i n_i$`), never raw Unicode subscripts (they render as boxes).
- Nomenclature (hard rule): `RhoA-GDP / -GTP-1 / -GTP-2 / -GTP-2A / -GTP-3A`; never bare 3A/2A/GTP.
  Fix display text only (`.replace("-3A","-GTP-3A").replace("-2A","-GTP-2A")`); never edit data keys or
  file paths.

## 4. Fonts (readability is enforced)
Match the paper's visual weight; single-panel figures use large text (ticks ~16, labels ~18, legend
~11). The following are **minimums** — go bigger whenever there is room; only approach the floor when
space is tight, and prefer enlarging the figure over shrinking text:
- tick labels & legend: **≥ 8 pt (prefer 9)**
- axis labels: **≥ 10–11 pt**
- panel letters: ~13 pt bold
Never go below 8 pt.

## 5. Bars / legends / error bars
- Legends: `frameon=False`; never overlapping bars/data — move to clear space, a top strip, or a
  shared legend above the panel row.
- Bars: normal thickness (`width` ~0.5–0.6 + x-padding); value labels where useful, white `bbox` when
  over hatched/colored bars; use hatch for a 2nd categorical dimension.
- Error bars: **default SD**. Use **SEM only when explicitly requested**; for correlated MD data
  compute SEM **block-averaged** (naive SD/√N over frames underestimates). Always state which (SD or
  SEM, and block count) in the caption.

## 6. Page fit & ACS specs
- **Width:** ACS double-column **≤ 7.0 in** (17.78 cm); single-column **3.33 in** (8.25 cm);
  1.5-column ~4.75 in. Max height ~9.5 in. Our multi-panel SI figures target ≤ ~6.9–7.0 in wide.
- **Resolution:** ≥ 300 dpi (halftone); 600 dpi for combination line+text; vector (SVG/PDF/EPS)
  preferred for line art. We always export SVG + a 300+ dpi PNG.
- **Fonts/lines:** Arial/Helvetica; ACS min text ~4.5 pt — but our floor is 8 pt (§4). Line weights
  ≥ 0.5 pt at final size (our `axes.linewidth=2`, tick width 1.5 satisfy this).
- Final-submission files per the ACS graphics/Preparing-Figures appendix on pubs.acs.org (verify
  before submitting, since ACS revises).

## 7. Figure → script catalogue (what to edit for each figure)
*(Seed list — expand as needed. Paths relative to the project root unless noted EXT = OneDrive paper
folder.)*

| Figure | Script(s) |
|---|---|
| Fig 6 (composite) | `rhoa_full_sampled_trajs/compose_figure6.py` ← `fig_A_transition_states.py`, `rhoa_full_sampled_trajs/fig_B_amine_frequency.py`, `rhoa_full_sampled_trajs/fig_cobinding_network.py`, `rhoa_full_sampled_trajs/fig_transition_tiltrot.py` |
| Fig 7 (composite) | `rhoa_full_sampled_trajs/compose_figure7.py` ← `fig_dynetan_pathlength.py`, `switch_rmsd_effector.py`, `dynetan_community_track.py`, `render_community_network.py`, `fig_D_effector_rmsd.py` |
| SI S1 | `SI_replots/figS1/replot_figS1.py` |
| SI S3 | `SI_replots/figS3/plot_figS3.py` |
| SI S5 | `SI_replots/figS5/plot_gr_3d.py` |
| SI S6 | `SI_replots/figS6/replot_figS6.py` |
| SI S7–S10 | `rhoa_full_sampled_trajs/compose_figure6.py` / `compose_figure7.py` (build_* fns) + their panel scripts |
| Membrane APL (SI) | `membrane_apl_compare/fig_apl_psm_check.py` |
| Insertion vs thickness | `insertion_vs_thickness/plot_insertion_thickness.py` |
| Combined membrane SI | `SI_membrane_lipid_combined/si_membrane_lipid_figure.py` |
| Main-text g(r) / H-bonds | EXT `strided_trajs/hbonds_gr_analysis.py` (defines the house style) |

Run figures in the `memb_analysis` conda env:
`source ~/anaconda3/etc/profile.d/conda.sh && conda activate memb_analysis && python <script>`.
