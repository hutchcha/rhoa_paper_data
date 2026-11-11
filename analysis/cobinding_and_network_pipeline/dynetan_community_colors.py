#!/usr/bin/env python3
"""Shared WT-matched dynamical-community coloring so the sequence track and the 3-D structure
renders use IDENTICAL colors. Raw Louvain community indices are arbitrary per sim; we color every
residue by the reference (WT = GTP-2) community its own community overlaps most.

  from dynetan_community_colors import community_colors
  refcol, per_sim_resid_rgb = community_colors(D, ORDER)   # D = dynetan dir with <sim>_nodes.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REF = "GTP-2"
PALETTE = [tuple(c) for c in (list(plt.cm.tab10.colors) + list(plt.cm.Set3.colors))]


def community_colors(D, order):
    nd = {s: pd.read_csv(f"{D}/{s}_nodes.csv")[["resid", "community"]] for s in order}
    lo = min(int(d.resid.min()) for d in nd.values())
    hi = max(int(d.resid.max()) for d in nd.values())
    refd = nd[REF].set_index("resid")["community"]
    # reference communities ordered by first appearance along the sequence -> stable colors
    seen, ref_order = set(), []
    for r in range(lo, hi + 1):
        if r in refd.index:
            c = int(refd[r])
            if c not in seen:
                seen.add(c); ref_order.append(c)
    refcol = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(ref_order)}
    ref_members = {c: set(refd.index[refd == c]) for c in ref_order}

    def map_to_ref(simd):
        m = {}
        for c in simd.community.unique():
            mem = set(simd.resid[simd.community == c])
            m[int(c)] = max(ref_order, key=lambda rc: len(mem & ref_members[rc]))
        return m

    per_sim = {}
    for s in order:
        d = nd[s]
        smap = {int(c): int(c) for c in d.community.unique()} if s == REF else map_to_ref(d)
        per_sim[s] = {int(r): refcol[smap[int(c)]]
                      for r, c in zip(d.resid, d.community)}
    return refcol, ref_order, per_sim, (lo, hi)
