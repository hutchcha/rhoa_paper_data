#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# Your line plot function
# =========================

def line_plot(xvalues, yvalues, ax, title, color, z, label, ma_window=5, ma_color='red', ma_z=2):

    xvalues = np.array(xvalues)
    yvalues = np.array(yvalues)

    y_series = pd.Series(yvalues)
    moving_avg = y_series.rolling(window=ma_window, center=True).mean()

    ax.plot(xvalues, yvalues, alpha=.25, color=color, zorder=z, label=None)
    ax.plot(xvalues, moving_avg, alpha=1.0, color=ma_color, zorder=ma_z)

    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', labelsize=12)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.minorticks_off()
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(width=1.5, direction='out')


# =========================
# USER INPUT
# =========================

# Tilt files
tilt_files = [
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GDP-tilt.txt",
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GTP-tilt.txt",
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GTP-2-tilt.txt"
]

# Rotation files
rot_files = [
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GDP-rotation.txt",
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GTP-rotation.txt",
    "/mnt/data/orientation_dynamics_figures/analysis/Tilt and Rotation/RhoA-GTP-2-rotation.txt"
]

# labels and colors
labels = ["GDP", "GTP-1", "GTP-2"]
colors = ["orange", "blue", "green"]

# total simulation time (μs)
total_time = 20

# moving average window
ma_window = 100


# =========================
# Helper function
# =========================

def load_last_half(file):

    data = np.loadtxt(file)
    half = len(data) // 2
    print(len(data))
    return data


# =========================
# Load data
# =========================

tilts = [load_last_half(f) for f in tilt_files]
rots = [load_last_half(f) for f in rot_files]

# time axes
time_axes = [
    np.linspace(0, total_time, len(r)) for r in rots
]


# =========================
# Plot
# =========================

fig, ax = plt.subplots(2,1, figsize=(6,3.5))

# ----- rotation -----
for i in range(3):

    line_plot(
        time_axes[i],
        rots[i],
        ax=ax[0],
        title=None,
        color=colors[i],
        z=1,
        label=labels[i],
        ma_window=ma_window,
        ma_color=colors[i],
        ma_z=2
    )

ax[0].set_ylabel("φ (°)", fontsize=15)
ax[0].set_xticks([])
ax[0].set_ylim(-200,200)
# ----- tilt -----
for i in range(3):

    line_plot(
        time_axes[i],
        tilts[i],
        ax=ax[1],
        title=None,
        color=colors[i],
        z=1,
        label=labels[i],
        ma_window=ma_window,
        ma_color=colors[i],
        ma_z=2
    )

ax[1].set_ylabel("θ (°)", fontsize=15)
ax[1].set_ylim(0,200)
fig.supxlabel("Time (μs)", fontsize=15)

plt.tight_layout()

fig.savefig("RhoA/strided_trajs/RhoA_time_evo_tiltrot_last_half.png", dpi=600)

plt.show()