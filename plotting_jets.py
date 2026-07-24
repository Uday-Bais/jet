import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

file_path = "full_history_pythia_20gev_ptHat.root"
tree_name = "ParticleTree"
event_id = 0
pt_min = 1 


branches = ["EventID", "ParticleID", "IsFinal", "Px", "Py", "Pz", "Energy"]

def load_event(path, tree_name, event_id):
    tree = uproot.open(path)[tree_name]
    p = tree.arrays(branches)
    return p[p["EventID"] == event_id]

def add_trees(ev):
    px, py, pz, e = ev["Px"], ev["Py"], ev["Pz"], ev["Energy"]
    pt = np.sqrt(px**2 + py**2)
    
    # Calculating Rapidity
    e_plus = np.maximum(e + pz, 1e-10)
    e_minus = np.maximum(e - pz, 1e-10)
    y = 0.5 * np.log(e_plus / e_minus)

 
    phi = np.arctan2(py, px)
    phi_wrapped = ak.where(phi < 0, phi + 2 * np.pi, phi)

    ev = ak.with_field(ev, pt, "pT_calc")
    ev = ak.with_field(ev, y, "y_calc")
    ev = ak.with_field(ev, phi_wrapped, "phi_calc")
    return ev


event0 = load_event(file_path, tree_name, event_id)
event0 = add_trees(event0)

is_final = event0["IsFinal"] == True
particle_ids = np.abs(ak.to_numpy(event0["ParticleID"]))
is_visible = ak.Array(~np.isin(particle_ids, [12, 14, 16]))
is_above_pt = event0["pT_calc"] >= pt_min

final_pts = event0[is_final & is_visible & is_above_pt]

\
fig, ax = plt.subplots(figsize=(10, 8))

sc = ax.scatter(ak.to_numpy(final_pts["y_calc"]), ak.to_numpy(final_pts["phi_calc"]),
                c=ak.to_numpy(final_pts["pT_calc"]), cmap="viridis",
                s=60, alpha=0.9, edgecolor="black", linewidth=0.5,
                label="Final state particles")

cbar = plt.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label("Transverse Momentum $p_T$ (GeV)")

#Manually adding the jets properties.
jets = [
    {"label": "Jet 0", "pt": 20.7728, "y": -3.07519, "phi": 2.22168},
    {"label": "Jet 1", "pt": 21.3747, "y": -1.46637, "phi": 5.45769} 
]

jet_y = [j["y"] for j in jets]
jet_phi = [j["phi"] for j in jets]


ax.scatter(jet_y, jet_phi,
           s=400, facecolors='none', edgecolors='red', linewidth=3, marker='s',
           label="Jet centers")


for j in jets:
    ax.text(j["y"], j["phi"] + 0.3, f'{j["label"]}\n$p_T$={j["pt"]:.1f} GeV', 
            color='red', fontsize=11, ha='center', va='bottom', weight='bold')

ax.set_xlim(-10, 10)
ax.set_ylim(-0.5, 2 * np.pi + 0.5) 
ax.set_xlabel("Rapidity ($y$)")
ax.set_ylabel("Azimuthal Angle ($\\phi$)")

pi_ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
pi_labels = ['0', '$\\pi/2$', '$\\pi$', '$3\\pi/2$', '$2\\pi$']
ax.set_yticks(pi_ticks)
ax.set_yticklabels(pi_labels)

ax.set_title(f"Event {event_id}: Phase_Space_pTmin of 20 GeV ($p_T \geq {pt_min}$ GeV of final particles)")
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend(loc="upper right", frameon=True)

plot_path = "event_with_jets_rapidity_scatter.png"
plt.savefig(plot_path, bbox_inches="tight", dpi=300)
print(f"Plot saved to {plot_path}")
