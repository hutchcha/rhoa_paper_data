#!/usr/bin/env python3
"""Organize the supplementary figures into one directory per planned SI figure (S6-S10).

Each figSN_* directory gets the assembled composite plus the component plots that make it up
(PNG + SVG where available), copied from the flat supplementary/ (and figure7/) folders so the
originals stay in place for the plotting scripts. Re-runnable.

  python organize_si.py
"""
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
PF = os.path.abspath(os.path.join(HERE, "..", "paper_figures"))
SUPP = os.path.join(PF, "supplementary")
F7 = os.path.join(PF, "figure7")

# dir -> (composite stem[, dir], [component stems...])   stems resolved in SUPP unless (stem, F7)
PLAN = {
    "figS6_structural_stability": {
        "composite": "supp_2A_global_similarity",
        "components": [],
    },
    "figS7_refs_cobinding": {
        "composite": "supp_S7_refs_cobinding",
        "components": ["supp_hbond_occupancy_nontransition",
                       "supp_cobinding_frequency_nontransition",
                       "supp_cobinding_network_nontransition"],
    },
    "figS8_swII_network_quant": {
        "composite": "supp_S8_swII_network_quant",
        "components": ["supp_effector_competence_swII",
                       ("fig_dynetan_pathlength", F7),
                       ("fig_dynetan_betweenness", F7)],
    },
    "figS9_refs_switch_competence": {
        "composite": "supp_S9_refs_switch_competence",
        "components": ["supp_swI_competence_refs", "supp_swII_competence_refs"],
    },
    "figS10_refs_network": {
        "composite": "supp_S10_refs_network",
        "components": ["supp_dynetan_community_track", "supp_dynetan_community_network",
                       "supp_dynetan_pathlength_ref", "supp_dynetan_betweenness_ref"],
    },
}


def resolve(stem):
    return (stem, SUPP) if isinstance(stem, str) else stem


def copy_stem(stem, srcdir, dstdir, tag):
    got = []
    for ext in ("png", "svg"):
        src = os.path.join(srcdir, f"{stem}.{ext}")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dstdir, f"{tag}_{stem}.{ext}"))
            got.append(ext)
    return got


def main():
    for d, spec in PLAN.items():
        dst = os.path.join(SUPP, d)
        os.makedirs(dst, exist_ok=True)
        cstem, csrc = resolve(spec["composite"])
        got = copy_stem(cstem, csrc, dst, "COMPOSITE")
        print(f"{d}/")
        print(f"    composite  {cstem}  [{','.join(got) or 'MISSING'}]")
        for comp in spec["components"]:
            stem, src = resolve(comp)
            got = copy_stem(stem, src, dst, "part")
            print(f"    part       {stem}  [{','.join(got) or 'MISSING'}]")
    print(f"\norganized under {SUPP}")


if __name__ == "__main__":
    main()
