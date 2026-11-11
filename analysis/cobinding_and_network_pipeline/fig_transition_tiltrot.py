#!/usr/bin/env python3
"""Fig 6 Panel A: the RhoA-3A orientational transition. Left = tilt (theta) and rotation (phi)
time series with the transition marked by a dotted line; right = polar density of the settled
(post-transition) orientation (radial = tilt theta, angular = rotation phi), matching the polar
maps used earlier in the paper. Only 3A is shown here (the GTP-2 transition appears earlier in the
paper). Establishes that RhoA-3A transitions into the same stable membrane configuration.

Tilt/rotation traces from orientation_out*/tilt_rotation/RhoA-3A_{tilt,rotation}.txt; the crossing
window (transition) from stable_state_config.STATES.

  python fig_transition_tiltrot.py -> paper_figures/figure6/fig_transition_tiltrot.{png,svg}
"""
import warnings; warnings.filterwarnings("ignore")
import os, sys
import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import stable_state_config as C

plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(REPO, "paper_figures", "figure6")
SIM = "3A"
TILT_COL = "#333333"; ROT_COL = "#B279A2"
SMOOTH_NS = 100.0


def movavg(x, n):
    if n <= 1:
        return np.asarray(x, float)
    x = np.asarray(x, float); k = np.ones(n) / n; pad = n // 2
    return np.convolve(np.pad(x, (pad, pad), mode="edge"), k, mode="same")[pad:pad + len(x)]


def load(name):
    for root in ("orientation_out_mut", "orientation_out"):
        t = f"{REPO}/{root}/tilt_rotation/RhoA-{name}_tilt.txt"
        r = f"{REPO}/{root}/tilt_rotation/RhoA-{name}_rotation.txt"
        if os.path.exists(t) and os.path.exists(r):
            return np.loadtxt(t), np.loadtxt(r)
    return None, None


def polar_density(ax, rot_deg, tilt_deg, cmap="cividis", tilt_max=180):
    """Smoothed polar density (radial=tilt, angular=rotation) with circular wrap in phi."""
    tb = np.linspace(-np.pi, np.pi, 361)
    rb = np.linspace(0, np.deg2rad(tilt_max), 101)
    H, _, _ = np.histogram2d(np.deg2rad(rot_deg), np.deg2rad(tilt_deg), bins=[tb, rb])
    pad = 20
    Hs = gaussian_filter(np.pad(H, ((pad, pad), (0, 0)), mode="wrap"), sigma=(5, 2), mode="wrap")
    Hs = Hs[pad:-pad, :].T
    tc = (tb[:-1] + tb[1:]) / 2; rc = (rb[:-1] + rb[1:]) / 2
    T, R = np.meshgrid(tc, rc)
    Hs /= Hs.max() if Hs.max() > 0 else 1.0
    ax.set_theta_zero_location("N"); ax.set_theta_direction("clockwise")
    ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(["0°", "45°", "90°", "135°", "180°", "-135°", "-90°", "-45°"], fontsize=8)
    ax.set_yticks(np.deg2rad([30, 60, 90, 120, 150]))
    ax.set_yticklabels(["30", "", "90", "", "150"], fontsize=8)
    ax.grid(linestyle="-", color="black", alpha=0.4, zorder=1)
    ax.contourf(T, R, Hs, levels=np.linspace(0.02, 1.0, 11), cmap=cmap, zorder=0)


def main():
    os.makedirs(OUT, exist_ok=True)
    dt = next(s["dt_ns"] for s in C.SIMS if s["name"] == SIM)
    tilt, rot = load(SIM)
    x = np.arange(len(tilt)) * dt / 1000.0
    sm = max(1, int(SMOOTH_NS / dt))
    cross = C.STATES.get(SIM, {}).get("cross")
    post = C.STATES.get(SIM, {}).get("post")
    tstar = 0.5 * (cross[0] + cross[1]) if cross else None
    f0 = int(post[0] * 1000 / dt); f1 = int(min(post[1] * 1000 / dt, len(tilt)))

    # exact 130 x 71 mm box; constrained_layout fits all content inside (no tight-bbox trimming)
    fig = plt.figure(figsize=(130 / 25.4, 71 / 25.4), layout="constrained")
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.9, 1.0])
    axt = fig.add_subplot(gs[0, 0])
    axt.plot(x, movavg(tilt, sm), color=TILT_COL, lw=1.4, label="tilt θ")
    axt.plot(x, movavg(rot, sm), color=ROT_COL, lw=1.2, alpha=0.85, label="rotation φ")
    axt.axhline(0, color="0.85", lw=0.7, zorder=0)
    if tstar is not None:
        axt.axvline(tstar, color="#E45756", ls=":", lw=1.8)
        axt.text(tstar, 0.98, "  transition", transform=axt.get_xaxis_transform(),
                 ha="left", va="top", fontsize=8.5, color="#E45756")
    axt.set_xlim(0, x[-1]); axt.set_ylim(-90, 120)
    axt.set_xlabel("time (µs)", fontsize=11); axt.set_ylabel("angle (deg)", fontsize=11)
    axt.set_title(f"RhoA-{SIM}".replace("-3A", "-GTP-3A").replace("-2A", "-GTP-2A"),
                  fontsize=13, fontweight="bold", loc="left")
    axt.tick_params(labelsize=10)
    axt.spines[["right", "top"]].set_visible(False)
    axt.legend(frameon=False, fontsize=9, ncol=2, loc="upper right")

    axp = fig.add_subplot(gs[0, 1], projection="polar")
    polar_density(axp, rot[f0:f1], tilt[f0:f1])
    axp.set_title("settled orientation", fontsize=9.5, color="0.35", pad=10)

    for e in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_transition_tiltrot.{e}", dpi=300)   # no tight bbox -> output == 130x71 mm
    plt.close(fig)
    print(f"wrote {OUT}/fig_transition_tiltrot.png")


if __name__ == "__main__":
    main()
