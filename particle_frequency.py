import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from particle import Particle, ParticleNotFound


file = uproot.open("full_history_pythia_20gev_ptHat.root")
tree = file["ParticleTree"]
arrays = tree.arrays()

is_event_0 = arrays.EventID == 0
is_final = arrays.IsFinal == True
is_neutrino = (abs(arrays.ParticleID) == 12) | (abs(arrays.ParticleID) == 14) | (abs(arrays.ParticleID) == 16)

detected_event_0 = arrays[is_event_0 & is_final & ~is_neutrino]

pdg_ids = detected_event_0.ParticleID.to_numpy()


unique_ids, counts = np.unique(pdg_ids, return_counts=True)


sort_indices = np.argsort(counts)[::-1]
unique_ids = unique_ids[sort_indices]
counts = counts[sort_indices]


labels = []
for pid in unique_ids:
    try:
        p = Particle.from_pdgid(pid)
        sym = f"${p.latex_name}$"
    except ParticleNotFound:
        sym = ""
    
    label = f"{sym} ({pid})".strip() if sym else f"({pid})"
    labels.append(label)


plt.figure(figsize=(14, 7))
bars = plt.bar(labels, counts, color="teal", alpha=0.8, edgecolor="black")

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.5, int(yval), ha='center', va='bottom', fontsize=12)

plt.title("Detected Final-State Particle Frequencies (Event 0)", fontsize=22, fontweight='bold', pad=15)
plt.xlabel("Particle Type (PDG ID)", fontsize=22, labelpad=10)
plt.ylabel("Frequency", fontsize=22, labelpad=10)

plt.xticks(rotation=45, ha='right', fontsize=22)
plt.yticks(fontsize=22)

plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

print(f"Total detected final state particles in Event 0: {len(pdg_ids)}")
print("\nTop Particle Frequencies:")
for label, count in zip(labels, counts):
    print(f"  {label}: {count}")
