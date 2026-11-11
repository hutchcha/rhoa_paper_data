#!/usr/bin/env python3

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ============================================================
# INPUT
# ============================================================

RDF_CSV = "parallel_rdfs_long.csv"

OUTDIR = Path("hbonds_out/no_stride_rdf_plots")

SIM_NAMES = ["GDP", "GTP-1", "GTP-2"]

# Map names from CSV -> plotting names
SIM_NAME_MAP = {
    "GDP": "GDP",
    "GTP": "GTP-1",
    "Resim": "GTP-2",
}

LIPID_RESNAMES = ["POPC", "POPE", "POPS", "PLA18"]

RESIDUES = [5, 6, 7]

RESID_TO_NAME = {
    5: "Arg5",
    6: "Lys6",
    7: "Lys7",
}

RESID_LABEL_MAP = {
    5: "Arg5",
    6: "Lys6",
    7: "Lys7",
}

# ============================================================
# OI COLOR PALETTE
# ============================================================

LIPID_COLORS = {
    "POPC": "#1F77B4",  # blue
    "POPE": "#FF7F0E",  # orange
    "POPS": "#2CA02C",  # green
    "PLA18": "#D62728", # red
}

FIGSIZE = (5.2, 7.2)
LINEWIDTH = 2.0
SEM_ALPHA = 0.20


# ============================================================
# AXIS STYLE
# ============================================================

def style_axes(ax):

    ax.grid(False)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.spines["left"].set_linewidth(1.25)
    ax.spines["bottom"].set_linewidth(1.25)

    ax.tick_params(
        axis="both",
        labelsize=16,
        width=1.25,
        length=3,
        direction="out",
    )


# ============================================================
# PLOTTING
# ============================================================

def plot_rdf_single_residue_all_sims_blockavg(
    rdf_df,
    resid,
    outpath,
):

    resid_name = RESID_TO_NAME[resid]

    n_sims = len(SIM_NAMES)

    fig, axes = plt.subplots(
        nrows=n_sims,
        ncols=1,
        figsize=FIGSIZE,
        sharex=True,
        sharey=True,
    )

    if n_sims == 1:
        axes = [axes]

    # --------------------------------------------------------
    # global ymax
    # --------------------------------------------------------

    global_ymax = 0
    found_any = False

    for sim in SIM_NAMES:

        for lip in LIPID_RESNAMES:

            sub = rdf_df[
                (rdf_df.sim_name == sim) &
                (rdf_df.residue_name == resid_name) &
                (rdf_df.lipid_name == lip)
            ]

            if sub.empty:
                continue

            y = sub.rdf_mean.values

            if np.any(np.isfinite(y)):

                global_ymax = max(global_ymax, np.nanmax(y))
                found_any = True

    if not found_any:
        raise ValueError("No RDF data found")

    # --------------------------------------------------------
    # plot
    # --------------------------------------------------------

    for ax, sim in zip(axes, SIM_NAMES):

        sim_found = False

        for lip in LIPID_RESNAMES:

            sub = rdf_df[
                (rdf_df.sim_name == sim) &
                (rdf_df.residue_name == resid_name) &
                (rdf_df.lipid_name == lip)
            ]

            if sub.empty:
                continue

            sub = sub.sort_values("r_A")

            r = sub.r_A.values
            g = sub.rdf_mean.values
            sem = sub.rdf_sem.values

            color = LIPID_COLORS[lip]

            ax.plot(
                r,
                g,
                color=color,
                linewidth=LINEWIDTH,
                label=lip if ax is axes[0] else None
            )

            ax.fill_between(
                r,
                g - sem,
                g + sem,
                color=color,
                alpha=SEM_ALPHA
            )

            sim_found = True

        style_axes(ax)

        ax.set_ylabel("g(r)", fontsize=18)

        ax.set_title(sim, fontsize=15)

        ax.set_ylim(0, global_ymax * 1.05)

        if not sim_found:

            ax.text(
                0.5,
                0.5,
                "No data",
                transform=ax.transAxes,
                ha="center",
                va="center"
            )

    axes[-1].set_xlabel("r (Å)", fontsize=18)

    fig.suptitle(resid_name, fontsize=16, y=0.995)

    legend_handles = [
        plt.Line2D([0], [0], color=LIPID_COLORS[lip], lw=3, label=lip)
        for lip in LIPID_RESNAMES
    ]

    axes[0].legend(
        handles=legend_handles,
        frameon=False,
        fontsize=10,
        loc="upper right",
        title="Lipid"
    )

    fig.tight_layout()

    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    fig.savefig(outpath, bbox_inches="tight")

    plt.close(fig)

    print("Wrote", outpath)

    return fig


# ============================================================
# MAIN
# ============================================================

def main():

    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(RDF_CSV)

    # ------------------------------------------
    # Fix simulation names
    # ------------------------------------------

    df["sim_name"] = df["sim_name"].map(SIM_NAME_MAP)

    # ------------------------------------------
    # individual residue plots
    # ------------------------------------------

    for resid in RESIDUES:

        resid_name = RESID_TO_NAME[resid]

        outpath = OUTDIR / f"rdf_{resid_name}_all_sims.pdf"

        plot_rdf_single_residue_all_sims_blockavg(
            df,
            resid,
            outpath,
        )

    # ------------------------------------------
    # combined PDF
    # ------------------------------------------

    combined_pdf = OUTDIR / "all_residue_rdfs.pdf"

    with PdfPages(combined_pdf) as pdf:

        for resid in RESIDUES:

            fig = plot_rdf_single_residue_all_sims_blockavg(
                df,
                resid,
                OUTDIR / f"_tmp_{resid}.pdf",
            )

            pdf.savefig(fig, bbox_inches="tight")

    print("Wrote", combined_pdf)


if __name__ == "__main__":
    main()
