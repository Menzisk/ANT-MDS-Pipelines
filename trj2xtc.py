#  Converstion Script: "Name of System"
#  Section: Converts -trj folder into xtc for downstream analysis 
#  Author : Menzi Sikakane

import mdtraj as md
Conversion = "_trj to .xtc"

print(f"   — {Conversion}")
traj = md.load(
    "desmond_md_trj/clickme.dtr", # Use the .dtr file from the -trj folder as input
    top=".pdb" # Use the same PDB/GRO file as topology reference
)

print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
traj.save_xtc("_trajectory.xtc") # Save as .xtc for use in DCCM, dPCA, and FEL scripts
print("Saved: _trajectory.xtc")