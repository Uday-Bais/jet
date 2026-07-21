import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

file_path = "full_history_pythia_default_ptHat.root"
tree_name = "ParticleTree"
event_id = 0
target_idx = 5  # the hard-process gluon

branches = ["EventID", "Index", "ParticleID", "Status",
            "Mother1", "Mother2", "Daughter1", "Daughter2",
            "IsFinal", "Px", "Py", "Pz"]


def load_event(path, tree_name, event_id):
    tree = uproot.open(path)[tree_name]
    p = tree.arrays(branches)
    return p[p["EventID"] == event_id]


def add_trees(ev):
    px, py, pz = ev["Px"], ev["Py"], ev["Pz"]
    pt = np.sqrt(px**2 + py**2)
    pt_safe = ak.where(pt == 0, 1e-10, pt)

    ev = ak.with_field(ev, pt, "pT_calc")
    ev = ak.with_field(ev, np.arcsinh(pz / pt_safe), "eta_calc")
    ev = ak.with_field(ev, np.arctan2(py, px), "phi_calc")
    return ev


def get_descendants(root_idx, d1_map, d2_map):
    
    seen = set()
    stack = [root_idx]
    while stack:
        idx = stack.pop()
        d1 = d1_map.get(idx)
        d2 = d2_map.get(idx)
        if not d1 or d1 <= 0:
            continue
        if d2 > d1:
            kids = range(d1, d2 + 1)
        elif d2 in (0, d1):
            kids = (d1,)
        else:
            kids = (d1, d2)
        for k in kids:
            if k in d1_map and k not in seen:
                seen.add(k)
                stack.append(k)
    return seen


event0 = load_event(file_path, tree_name, event_id)
event0 = add_trees(event0)

indices = ak.to_numpy(event0["Index"])
d1_map = dict(zip(indices, ak.to_numpy(event0["Daughter1"])))
d2_map = dict(zip(indices, ak.to_numpy(event0["Daughter2"])))

desc = get_descendants(target_idx, d1_map, d2_map)
is_desc = ak.Array(np.isin(indices, list(desc)))
is_final = event0["IsFinal"] == True

hard_gluon = event0[event0["Index"] == target_idx]
target_pts = event0[is_final & is_desc]
other_pts = event0[is_final & ~is_desc]

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(ak.to_numpy(other_pts["eta_calc"]), ak.to_numpy(other_pts["phi_calc"]),
           color="lightgrey", s=40, alpha=0.5, label="Other final state particles")

sc = ax.scatter(ak.to_numpy(target_pts["eta_calc"]), ak.to_numpy(target_pts["phi_calc"]),
                 c=ak.to_numpy(target_pts["pT_calc"]), cmap="viridis",
                 s=60, alpha=0.9, edgecolor="black", linewidth=0.5,
                 label=f"Final state constituents of gluon {target_idx}")

ax.scatter(ak.to_numpy(hard_gluon["eta_calc"]), ak.to_numpy(hard_gluon["phi_calc"]),
           c="red", marker="X", s=200, label=f"Original gluon (index {target_idx})")

cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Transverse Momentum $p_T$ (GeV)")

ax.set_xlim(-10, 10)
ax.set_xlabel("Pseudorapidity ($\\eta$)")
ax.set_ylabel("Azimuthal Angle ($\\phi$)")
ax.set_title(f"Event {event_id}: Final State Particles")
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=True)

plot_path = "single_gluon_scatter_uproot.png"
plt.savefig(plot_path, bbox_inches="tight", dpi=300)
print(f"Plot saved")