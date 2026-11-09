#  DCCM Pipeline "Name of System"
#  Author : Menzi Sikakane
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import align
import warnings

warnings.filterwarnings('ignore')
TOPOLOGY = '.pdb'
TRAJECTORY = '_trajectory.xtc'
SYSTEM_LABEL = ' System'  # e.g. " SaANT Apo"
SELECTION = 'protein and name CA' # CA only, one point per residue
REF_FRAME = 0 # frame everything gets aligned to
CMAP = 'viridis'
DPI = 800
OUTPUT_DIR = '.'

print(f'  DCCM PIPELINE  |  System: {SYSTEM_LABEL}')

print('\n[STEP 1] Loading trajectory')
u = mda.Universe(TOPOLOGY, TRAJECTORY)
print(f'         Topology : {TOPOLOGY}')
print(f'         Trajectory: {TRAJECTORY}')
print(f'         Frames loaded : {u.trajectory.n_frames}')

ag = u.select_atoms(SELECTION)
n_atoms = ag.n_atoms
resids = ag.resids
resnames = ag.resnames # not used below, keeping it for labelling later

print(f"Selection   : '{SELECTION}'")
print(f'Residues    : {n_atoms}  (resid {resids[0]}–{resids[-1]})')

print(f'\n[STEP 2] Aligning all frames to frame {REF_FRAME}')
# need to remove overall rotation/translation first, otherwise the
# correlations just reflect the whole protein tumbling around
ref = mda.Universe(TOPOLOGY, TRAJECTORY) # separate copy so the reference isn't moved by the alignment
aligner = align.AlignTraj(u, ref, select=SELECTION, in_memory=True, ref_frame=REF_FRAME).run()
print(f'Alignment complete.')

print(f'\n[STEP 3] Extracting Cα coordinates from {u.trajectory.n_frames} frames')
n_frames = u.trajectory.n_frames
coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32) # float32 to keep memory down on long runs
for i, ts in enumerate(u.trajectory):
    coords[i] = ag.positions.copy()
print(f'coords array shape: {coords.shape}')
print(f' = ({n_frames} frames) × ({n_atoms} atoms) × (3 xyz)')

print(f'\n[STEP 4] Computing DCCM')
# C_ij = <dr_i . dr_j> / sqrt(<dr_i^2> * <dr_j^2>)
mean_pos = coords.mean(axis=0) # average structure
print(f' Mean position array shape: {mean_pos.shape}')

delta = coords - mean_pos # how far each CA sits from its average, per frame
print(f'Displacement array shape : {delta.shape}')
# t = frame, i/j = residues, a = xyz. dot product over xyz, averaged over frames
numerator = np.einsum('tia,tja->ij', delta, delta) / n_frames
print(f'Numerator matrix shape   : {numerator.shape}')
# same thing but i with itself, gives the mean squared fluctuation per residue
variance_i = np.einsum('tia,tia->i', delta, delta) / n_frames
denominator = np.sqrt(np.outer(variance_i, variance_i))
print(f'Denominator matrix shape : {denominator.shape}')

dccm = numerator / (denominator + 1e-12) # tiny offset so a frozen residue doesn't give nan
np.fill_diagonal(dccm, 1.0) # self-correlation is 1, just cleaning up rounding error
print(f'\n   DCCM computed.')
print(f' Shape : {dccm.shape}  ({n_atoms} × {n_atoms})')
print(f' Min   : {dccm.min():.4f}')
print(f' Max   : {dccm.max():.4f}')
print(f' Mean  : {dccm.mean():.4f}')

print(f'\n[STEP 5] Plotting DCCM heatmap')
fig, ax = plt.subplots(figsize=(8, 7))
# colour scale pinned to -1..1 with 0 in the middle so systems can be compared side by side
norm = mcolors.TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
im = ax.imshow(dccm, cmap=CMAP, norm=norm, origin='lower', aspect='equal', interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Cross-correlation coefficient  $C_{ij}$', fontsize=14, fontweight='bold')
cbar.ax.tick_params(labelsize=10)
n = dccm.shape[0]
step = max(1, n // 15) # roughly 15 ticks regardless of protein size
tick_positions = np.arange(0, n, step)
ax.set_xticks(tick_positions)
ax.set_yticks(tick_positions)
# tick labels show the real resid numbers, not the 0-based matrix index
ax.set_xticklabels(resids[tick_positions], rotation=45, ha='right', fontsize=15)
ax.set_yticklabels(resids[tick_positions], fontsize=15)
ax.set_xlabel('Residue index', fontsize=14, fontweight='bold')
ax.set_ylabel('Residue index', fontsize=14, fontweight='bold')
ax.set_title(f'Dynamical Cross-Correlation Matrix\n{SYSTEM_LABEL}', fontsize=14, fontweight='bold', pad=12)
ax.plot([0, n - 1], [0, n - 1], color='black', lw=0.5, ls='--', alpha=0.3) # faint guide along the diagonal
fig.tight_layout()
png_path = f'{OUTPUT_DIR}/{SYSTEM_LABEL}_DCCM.png'
plt.savefig(png_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f'PNG saved : {png_path}')

print(f'\n[STEP 6] Saving DCCM matrix to CSV')
# resids as row/column labels so it can be reloaded for other plots
dccm_df = pd.DataFrame(dccm, index=resids, columns=resids)
csv_path = f'{OUTPUT_DIR}/{SYSTEM_LABEL}_DCCM.csv'
dccm_df.to_csv(csv_path)
print(f'CSV saved : {csv_path}')

print(f'\n[STEP 7] Summary statistic')
mask_offdiag = ~np.eye(n_atoms, dtype=bool) # drop the diagonal, it's always 1 and skews the stats
off_vals = dccm[mask_offdiag]
# matrix is symmetric so every pair is counted twice, hence the /2
strong_pos = (off_vals > 0.5).sum() / 2
strong_neg = (off_vals < -0.5).sum() / 2
print(f' Off-diagonal correlations:')
print(f' Mean         : {off_vals.mean():.4f}')
print(f' Std dev      : {off_vals.std():.4f}')
print(f' Strong pos (>0.5)  pairs : {int(strong_pos)}')
print(f' Strong neg (<-0.5) pairs : {int(strong_neg)}')
print(f"\n{'=' * 65}")
print(f' DCCM pipeline complete for {SYSTEM_LABEL}')
print(f' Outputs in: {OUTPUT_DIR}/')
