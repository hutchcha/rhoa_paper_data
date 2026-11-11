#!/usr/bin/env python3
"""Mutual-information generalized-correlation matrices (Lange-Grubmuller) of Cα motion,
all sims, over the frame-matched full-resolution (0.1 ns) post-state excerpts.

Generalized correlation r_MI[i,j] = sqrt(1 - exp(-2 I[i,j] / 3)), with I the (linearized,
Gaussian) MUTUAL INFORMATION between the 3-D Cα displacement vectors of residues i and j:
    I = 0.5 * ln( det(Cii) det(Cjj) / det(Cij) )
where Cii/Cjj are 3x3 covariances and Cij the 6x6 joint covariance. Unlike the linear
DCCM (which only sees collinear motion and can read ~0 for real but orthogonal coupling),
r_MI captures coupling in ANY relative direction, in [0,1]. This is the substrate for the
allosteric network analysis.

Self-contained: uses only stable_state_config.py (SIMS/paths) from THIS bundle. Robust to
the path layout — falls back to <bundle>/<name>_last<POST_US>us.dcd if the config path is a
dangling symlink. Frame-matches every sim to the common minimum (MI is sample-size-biased).

  python dccm_mi.py [--sims 2A,3A,GTP-2,GTP-1,GDP] [--max-frames 0] [--stride 1]

Outputs (transition_states/dccm_mi/):
  mi_<sim>.npy            NxN generalized-correlation matrix
  mi_resids.npy           residue ids for the axes
  mi_all.png/.svg         all sims side by side
  mi_diff.png/.svg        difference vs WT (GTP-2), if present  -> what the mutation rewires
  npatch_mi_profile.*     N-patch(5-7) ↔ residue generalized correlation, regions highlighted
"""
import warnings; warnings.filterwarnings("ignore")
import os, argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
import stable_state_config as C

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "transition_states", "dccm_mi")

# invariant alignment core = G-domain (4-179) minus the mobile elements
SWITCH1 = (28, 42); SWITCH2 = (56, 78); RHO_INSERT = (120, 137); NPATCH = (5, 7)
CORE_SEL = (f"protein and name CA and resid 4-179 and not "
            f"(resid {SWITCH1[0]}-{SWITCH1[1]} or resid {SWITCH2[0]}-{SWITCH2[1]} "
            f"or resid {RHO_INSERT[0]}-{RHO_INSERT[1]})")
ALLCA = "protein and name CA and resid 4-189"
ORDER = ["GTP-2", "2A", "3A", "GTP-1", "GDP"]
# regions to shade / annotate (name, lo, hi)
HL = [("N-patch", 5, 7), ("Switch I", 28, 42), ("Switch II", 56, 78),
      ("Rho ins", 120, 137), ("α5", 164, 177), ("HVR", 180, 189)]


def resolve_traj(sim):
    """config traj path, or the bundle-local fallback if that path is missing."""
    t = sim["traj"]
    if os.path.exists(t):
        return t
    alt = os.path.join(HERE, f"{sim['name']}_last{C.POST_US}us.dcd")
    return alt if os.path.exists(alt) else t


def aligned_displacement(sim, nframes, stride):
    """core-aligned all-Cα displacement (T,N,3), mean-removed, + resids."""
    u = mda.Universe(sim["top"], resolve_traj(sim))
    core = u.select_atoms(CORE_SEL); allca = u.select_atoms(ALLCA)
    T = len(u.trajectory)
    sel = np.arange(0, T, stride)
    if nframes and len(sel) > nframes:                       # even coverage of the excerpt
        sel = np.linspace(0, len(sel) - 1, nframes).round().astype(int)
        sel = np.arange(0, T, stride)[sel]
    u.trajectory[int(sel[0])]
    ref = core.positions - core.positions.mean(0)
    P = np.empty((len(sel), len(allca), 3), np.float32)
    for k, f in enumerate(sel):
        u.trajectory[int(f)]
        c = core.positions; com = c.mean(0)
        R, _ = rotation_matrix(c - com, ref)
        P[k] = (allca.positions - com) @ R.T
    return P - P.mean(0), allca.resids.copy()


def mi_matrix(D):
    """Linearized mutual-information generalized correlation r_MI (N,N) from displacement D."""
    T, N, _ = D.shape
    Dflat = D.reshape(T, 3 * N).astype(np.float64)
    Sigma = (Dflat.T @ Dflat) / T                                       # (3N,3N) covariance
    Sigma += np.eye(3 * N) * 1e-8                                        # numerical floor
    Bmat = Sigma.reshape(N, 3, N, 3).transpose(0, 2, 1, 3)              # (N,N,3,3) blocks
    detblk = np.linalg.det(np.einsum('iikl->ikl', Bmat))               # (N,) det of Cii diagonal blocks
    iu = np.triu_indices(N, 1); ii, jj = iu
    A = np.empty((len(ii), 6, 6))
    A[:, :3, :3] = Bmat[ii, ii]; A[:, :3, 3:] = Bmat[ii, jj]
    A[:, 3:, :3] = Bmat[jj, ii]; A[:, 3:, 3:] = Bmat[jj, jj]
    detij = np.linalg.det(A)
    ratio = np.clip(detblk[ii] * detblk[jj] / np.clip(detij, 1e-30, None), 1.0, 1e12)
    I = 0.5 * np.log(ratio)
    r = np.sqrt(np.clip(1.0 - np.exp(-2.0 * I / 3.0), 0.0, 1.0))
    M = np.eye(N)
    M[ii, jj] = r; M[jj, ii] = r
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", default="2A,3A,GTP-2,GTP-1,GDP")
    ap.add_argument("--max-frames", type=int, default=0, help="0 = use common minimum across sims")
    ap.add_argument("--stride", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    want = [s for s in ORDER if s in args.sims.split(",")]
    present = [s for s in [x for x in C.SIMS if x["name"] in want]
               if os.path.exists(resolve_traj(s))]
    missing = [n for n in want if n not in [s["name"] for s in present]]
    if missing:
        print(f"!! missing trajectories, skipping: {missing}")
    if not present:
        print("no trajectories found — nothing to do."); return

    # frame-match: common minimum (after stride)
    counts = {}
    for s in present:
        counts[s["name"]] = len(mda.Universe(s["top"], resolve_traj(s)).trajectory) // args.stride
    nframes = args.max_frames or min(counts.values())
    print(f"frame-matching all sims to {nframes} frames (per-sim available: {counts})", flush=True)

    mats = {}; resids = None
    for s in present:
        D, resids = aligned_displacement(s, nframes, args.stride)
        M = mi_matrix(D)
        mats[s["name"]] = M
        np.save(f"{OUT}/mi_{s['name']}.npy", M)
        off = M[~np.eye(len(M), dtype=bool)]
        print(f"  [{s['name']}] r_MI {M.shape} from {D.shape[0]} frames  "
              f"(mean off-diag {off.mean():.3f})", flush=True)
    np.save(f"{OUT}/mi_resids.npy", resids)

    plot(mats, resids)


def _region_lines(ax, resids):
    r2i = {int(r): i for i, r in enumerate(resids)}
    for _, lo, hi in HL:
        idx = [r2i[r] for r in range(lo, hi + 1) if r in r2i]
        if idx:
            for a in (ax.axhline, ax.axvline):
                a(min(idx) - 0.5, color="k", lw=0.3, alpha=0.3)
                a(max(idx) + 0.5, color="k", lw=0.3, alpha=0.3)


def plot(mats, resids):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sims = [s for s in ORDER if s in mats]
    r2i = {int(r): i for i, r in enumerate(resids)}
    tickres = [10, 30, 50, 70, 90, 110, 130, 150, 170, 189]
    ticks = [r2i[r] for r in tickres if r in r2i]; tlab = [r for r in tickres if r in r2i]

    # --- all MI heatmaps ---
    fig, axes = plt.subplots(1, len(sims), figsize=(4.2 * len(sims), 4.6), constrained_layout=True)
    for ax, name in zip(np.atleast_1d(axes), sims):
        Mv = mats[name].copy(); np.fill_diagonal(Mv, np.nan)
        im = ax.imshow(Mv, cmap="magma", vmin=0, vmax=0.7, origin="lower", interpolation="nearest")
        _region_lines(ax, resids)
        ax.set_title(f"RhoA-{name}", fontsize=12)
        ax.set_xticks(ticks); ax.set_xticklabels(tlab, fontsize=7)
        ax.set_yticks(ticks); ax.set_yticklabels(tlab, fontsize=7); ax.set_xlabel("residue", fontsize=9)
    fig.colorbar(im, ax=list(np.atleast_1d(axes)), shrink=0.7, label="MI generalized correlation r$_{MI}$")
    fig.suptitle("Mutual-information generalized correlation (Cα), full-res post state", fontsize=13, fontweight="bold")
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/mi_all.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- difference vs WT ---
    if "GTP-2" in mats:
        diffs = [(n, mats[n] - mats["GTP-2"]) for n in ["2A", "3A", "GTP-1", "GDP"] if n in mats]
        if diffs:
            fig, axes = plt.subplots(1, len(diffs), figsize=(4.4 * len(diffs), 4.6), constrained_layout=True)
            for ax, (name, dm) in zip(np.atleast_1d(axes), diffs):
                dv = dm.copy(); np.fill_diagonal(dv, np.nan)
                im = ax.imshow(dv, cmap="PuOr_r", vmin=-0.4, vmax=0.4, origin="lower", interpolation="nearest")
                _region_lines(ax, resids)
                ax.set_title(f"RhoA-{name} − WT", fontsize=12)
                ax.set_xticks(ticks); ax.set_xticklabels(tlab, fontsize=7)
                ax.set_yticks(ticks); ax.set_yticklabels(tlab, fontsize=7); ax.set_xlabel("residue", fontsize=9)
            fig.colorbar(im, ax=list(np.atleast_1d(axes)), shrink=0.7, label="Δ r$_{MI}$ (sim − WT)")
            fig.suptitle("MI correlation difference vs WT — what the mutation rewires", fontsize=13, fontweight="bold")
            for e in ("png", "svg"):
                fig.savefig(f"{OUT}/mi_diff.{e}", dpi=150, bbox_inches="tight")
            plt.close(fig)

    # --- N-patch highlight profile ---
    npi = [r2i[r] for r in range(NPATCH[0], NPATCH[1] + 1) if r in r2i]
    COL = {"3A": "#4C78A8", "2A": "#F58518", "GTP-2": "#E45756", "GTP-1": "#72B7B2", "GDP": "#B279A2"}
    fig, ax = plt.subplots(figsize=(15, 4.8))
    for name in sims:
        M = mats[name]
        prof = M[npi].mean(0).copy()
        for i in npi:
            prof[i] = np.nan                                     # blank self
        ax.plot(resids, prof, lw=1.6, color=COL.get(name), label=f"RhoA-{name}")
    ymax = 0.6
    for nm, lo, hi in HL:
        if nm == "N-patch":
            ax.axvspan(lo, hi, color="#888", alpha=0.25)
        else:
            ax.axvspan(lo, hi, color="#ffe08a" if "Switch" in nm else "#d9d9d9", alpha=0.35)
        ax.text((lo + hi) / 2, ymax * 0.98, nm, ha="center", va="top", fontsize=8, rotation=0)
    ax.set_xlabel("residue", fontsize=12)
    ax.set_ylabel("N-patch (5-7) ↔ residue  r$_{MI}$", fontsize=12)
    ax.set_ylim(0, ymax)
    ax.set_title("Where the N-terminal patch shows correlated motion (MI generalized correlation)",
                 fontsize=13, loc="left")
    ax.legend(frameon=False, ncol=5, fontsize=10)
    ax.spines[["right", "top"]].set_visible(False)
    fig.tight_layout()
    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/npatch_mi_profile.{e}", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # print N-patch top distal partners per sim
    print("\n=== N-patch top distal generalized-correlation partners (resid: r_MI) ===")
    for name in sims:
        prof = mats[name][npi].mean(0)
        cand = [(int(resids[k]), prof[k]) for k in range(len(resids))
                if not (NPATCH[0] - 2 <= resids[k] <= NPATCH[1] + 2)]
        top = sorted(cand, key=lambda x: -x[1])[:8]
        print(f"  {name:6}: " + "  ".join(f"{r}:{v:.2f}" for r, v in top))
    print(f"\nwrote {OUT}/mi_all.png, mi_diff.png, npatch_mi_profile.png")


if __name__ == "__main__":
    main()
