import numpy as np
import pandas as pd
import MDAnalysis as mda
import lipyphilic
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle



# ——— your universes & names ———
rhebgdp = mda.Universe("/home/chutchins@uthouston.edu/projects/rheb/rhebgdpanalysis/20us/sys.psf",
                       "/home/chutchins@uthouston.edu/projects/rheb/rhebgdpanalysis/20us/cenwrapprod.dcd")
rhebgtp = mda.Universe("/home/chutchins@uthouston.edu/projects/rheb/rhebgtpanalysis/20us/sys.psf",
                       "/home/chutchins@uthouston.edu/projects/rheb/rhebgtpanalysis/20us/cenwrapprod.dcd")
rhogdp  = mda.Universe("/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/sys.psf",
                       "/home/chutchins@uthouston.edu/projects/rhoa/rhoagdpanalysis/20us/cenwrapprod.dcd")
rhogtp  = mda.Universe("/mnt/nas/chase_projects/rhoa/rhoagtpanalysis/20us/sys.psf",
                       "/mnt/nas/chase_projects/rhoa/rhoagtpanalysis/20us/cenwrapprod.dcd")
np.printoptions(threshold=np.inf)
system_names = ["Rheb-GDP", "Rheb-GTP", "RhoA-GDP", "RhoA-GTP"]
universes     = [rhebgdp, rhebgtp, rhogdp, rhogtp]

# APL settings
lips    = ["POPC","POPE","POPS","CHL1", "SAPI"]
step    = 10   # stride
out_pkl = "APL/all_sims_apl.pkl"

def compute_apl_all(unis, names, step=1000, out_pkl="APL/all_sims_apl.pkl"):
    import numpy as np
    import pandas as pd
    import lipyphilic
    from tqdm import tqdm
    import pickle

    results = {}

    for name, u in zip(names, unis):
        # 1) choose selection
        if name in ("Rheb-GDP", "Rheb-GTP"):
            selection = "resname POPC SAPI CHL1 POPE and name P O3"
        elif name in ("RhoA-GDP", "RhoA-GTP"):
            selection = "name P"
        else:
            raise ValueError(f"Unknown system {name}")

        print(f"Running APL for {name} with selection: {selection}")
        leaf = lipyphilic.AssignLeaflets(universe=u, lipid_sel=selection)
        leaf.run(step=step, verbose=True)

        areas = lipyphilic.analysis.area_per_lipid.AreaPerLipid(
            universe=u,
            lipid_sel=selection,
            leaflets=leaf.leaflets
        )
        areas.run(step=step, verbose=True)

        # 2) build the residue→row map
        lipid_ag    = u.select_atoms(selection)
        residue_ids = lipid_ag.residues.resids  # e.g. array([101,105,110,...])
        resid_to_row = { resid: idx for idx, resid in enumerate(residue_ids) }

        # 3) start building the DataFrame
        apl_df = pd.DataFrame()
        apl_df["Total"] = areas.areas.mean(axis=0)

        # 4) per-lipid columns using the correct row indices
        for lip in lips:
            # pick the right atom selection for this lipid
            if lip == "CHL1":
                ag_lip = u.select_atoms("resname CHL1 and name O3")
            else:
                ag_lip = u.select_atoms(f"resname {lip} and name P")

            lip_resids = ag_lip.residues.resids
            if lip_resids.size == 0:
                print(f"  • {lip} not in {name}, skipping.")
                continue

            # map to row indices, skipping any that aren’t in the map
            row_idxs = [resid_to_row[r] for r in lip_resids if r in resid_to_row]
            if not row_idxs:
                print(f"  • no matching rows for {lip} in {name}, skipping.")
                continue

            # slice and average
            lip_areas = areas.areas[row_idxs, :].mean(axis=0)
            apl_df[lip] = lip_areas

        # 5) upper/lower POPC the same way
        memcom_z    = u.select_atoms("resname POPC POPS POPE SAPI and name P").center_of_mass()[2]
        popcup_ids  = u.select_atoms(f"(resname POPC and name P) and prop z > {memcom_z}").residues.resids
        popclow_ids = u.select_atoms(f"(resname POPC and name P) and prop z < {memcom_z}").residues.resids

        upper_rows = [resid_to_row[r] for r in popcup_ids if r in resid_to_row]
        lower_rows = [resid_to_row[r] for r in popclow_ids if r in resid_to_row]

        if upper_rows:
            apl_df["Upper POPC"] = areas.areas[upper_rows, :].mean(axis=0)
        if lower_rows:
            apl_df["Lower POPC"] = areas.areas[lower_rows, :].mean(axis=0)

        results[name] = apl_df

    # pickle and return
    with open(out_pkl, "wb") as f:
        pickle.dump(results, f)

    return results


# ——— run & pickle ———
apl_results = compute_apl_all(universes, system_names,
                              step=1,
                              out_pkl=out_pkl)
