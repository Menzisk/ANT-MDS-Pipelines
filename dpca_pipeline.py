"""
dPCA analysis pipeline for protein backbone dihedral angles.

Extracts phi/psi angles from an MD trajectory, applies sin/cos encoding,
then performs PCA, clustering, FEL estimation, and convergence analysis.

Author: Menzi Sikakane
"""

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from datetime import datetime
import sys

# Loggging setup: how the pipeline records itself 
class Tee:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()


# User-editable configuration
TOPOLOGY = ".pdb"            # Topology file for the system
TRAJECTORY = "_trajectory.xtc"  # MD trajectory file
CONDITION = ""               # e.g. "AbANT_Apigenin"; used for output naming
N_PCS = 10                    # Number of principal components to retain
TEMP_K = 300.0                # Temperature used for FEL calculations
K_RANGE = range(1, 9)         # Range of cluster values for K in silhouette analysis


# Redirect stdout to a log file so all analysis output is recorded
sys.stdout = Tee(f"{CONDITION}_log.txt")

print("=" * 60)
print(f"  dPCA PIPELINE — {CONDITION}")
print(f"  Run date/time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

print("\n[SECTION 1] Loading trajectory")
u = mda.Universe(TOPOLOGY, TRAJECTORY) #Creates the Universal object u (complete simulation record)

n_frames = len(u.trajectory)  # Counts number frames loaded 
n_atoms = u.atoms.n_atoms     #Total atom count in the system (protein+ligand+solvent+ions)
n_residues_total = u.residues.n_residues  # Total residue count

print(f"  Topology      : {TOPOLOGY}")
print(f"  Trajectory    : {TRAJECTORY}")
print(f"  Frames loaded : {n_frames}")
print(f"  Total atoms   : {n_atoms}")
print(f"  Total residues: {n_residues_total}")

print("\n[SECTION 2] Extracting phi/psi dihedral angles")
protein = u.select_atoms("protein") #isolate protein atoms only
rama = Ramachandran(protein).run() #Measures phi/psi for every residue 

angles = rama.results.angles 
n_frames, n_residues, _ = angles.shape
angles_flat = angles.reshape(n_frames, -1)

print(f"  Angle array shape     : {angles.shape}  (frames, residues, 2)")

angle_cols = [
    f"{angle}_{i+1}"
    for i in range(n_residues)
    for angle in ("phi", "psi")
]

# Store raw dihedral angles before any transformation or feature engineering
# to allow later inspection and reproducibility checks.
df_angles = pd.DataFrame(angles_flat, columns=angle_cols)
df_angles.insert(0, "frame", np.arange(n_frames))
df_angles.to_csv(f"{CONDITION}_raw_angles.csv", index=False)
print(f"  Residues with phi/psi : {n_residues}")
print(f"  Saved: {CONDITION}_raw_angles.csv")

print("\n[SECTION 3] Applying sin/cos transformation")
sincos_cols = []

for i in range(n_residues):
    phi_col = f"phi_{i+1}"
    psi_col = f"psi_{i+1}"

    phi_rad = np.deg2rad(df_angles[phi_col].values)
    psi_rad = np.deg2rad(df_angles[psi_col].values)

    # Convert circular angle values to continuous features; this avoids
    # artificial discontinuities near ±180° in the PCA analysis.
    df_angles[f"sin_phi_{i+1}"] = np.sin(phi_rad)
    df_angles[f"cos_phi_{i+1}"] = np.cos(phi_rad)
    df_angles[f"sin_psi_{i+1}"] = np.sin(psi_rad)
    df_angles[f"cos_psi_{i+1}"] = np.cos(psi_rad)

    sincos_cols.extend([
        f"sin_phi_{i+1}",
        f"cos_phi_{i+1}",
        f"sin_psi_{i+1}",
        f"cos_psi_{i+1}",
    ])

df_angles.to_csv(f"{CONDITION}_sincos_features.csv", index=False)

print(f"  Sin/cos features added : {len(sincos_cols)} columns")
print(f"  Saved: {CONDITION}_sincos_features.csv")

print("\n[SECTION 4] Standardising features and running PCA")
X = df_angles[sincos_cols].values  # shape: (n_frames, 4*n_rama_res)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X) # zero mean, unit variance

pca = PCA(n_components=N_PCS, random_state=42, svd_solver="full")
X_pca = pca.fit_transform(X_scaled) # shape: (n_frames, N_PCS)

pc_cols = [f"PC{i+1}" for i in range(N_PCS)]
for i, col in enumerate(pc_cols):
    df_angles[col] = X_pca[:, i]

var_df = pd.DataFrame({
    "PC": pc_cols,
    "explained_variance_ratio": pca.explained_variance_ratio_,
    "cumulative_variance": np.cumsum(pca.explained_variance_ratio_),
})

eig_df = pd.DataFrame({
    "PC": pc_cols,
    "eigenvalue": pca.explained_variance_,
    "explained_variance_ratio": pca.explained_variance_ratio_,
})

scree_df = pd.DataFrame({
    "PC": pc_cols,
    "eigenvalue": pca.explained_variance_,
})

loadings_df = pd.DataFrame(
    pca.components_.T,
    index=sincos_cols,
    columns=pc_cols,
)

pc_proj_df = pd.DataFrame(X_pca, columns=pc_cols)

scaler_df = pd.DataFrame({
    "feature": sincos_cols,
    "mean": scaler.mean_,
    "std": scaler.scale_,
})

var_df.to_csv(f"{CONDITION}_dPCA_explained_var.csv", index=False)
eig_df.to_csv(f"{CONDITION}_dPCA_eigenvalues.csv", index=False)
scree_df.to_csv(f"{CONDITION}_dPCA_scree.csv", index=False)
loadings_df.to_csv(f"{CONDITION}_dPCA_loadings.csv")
pc_proj_df.to_csv(f"{CONDITION}_dPCA_projections.csv", index=False)
scaler_df.to_csv(f"{CONDITION}_dPCA_scaler.csv", index=False)

print(f"  PC1 variance            : {pca.explained_variance_ratio_[0]*100:.2f}%")
print(f"  PC2 variance            : {pca.explained_variance_ratio_[1]*100:.2f}%")
print(f"  PC1+PC2                 : {np.sum(pca.explained_variance_ratio_[:2])*100:.2f}%")
print(f"  Cumulative PC1-PC{N_PCS}: {np.sum(pca.explained_variance_ratio_)*100:.2f}%")
print(f"  PC1 eigenvalue          : {pca.explained_variance_[0]:.4f}")
print(f"  PC2 eigenvalue          : {pca.explained_variance_[1]:.4f}")
print(f"  Saved: explained_var, eigenvalues, scree, loadings, projections, scaler")

print("\n[SECTION 5a] Silhouette analysis")
X_cluster = X_pca[:, :2]
silhouette_avgs = []
valid_k = []

for k in K_RANGE:
    km_test = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels_k = km_test.fit_predict(X_cluster)
    n_labels = len(np.unique(labels_k))

    if n_labels < 2:
        print(f"  k={k}: skipped (only {n_labels} cluster found)")
        continue

    avg_score = silhouette_score(X_cluster, labels_k)
    silhouette_avgs.append(avg_score)
    valid_k.append(k)
    print(f"  k={k}: average silhouette score = {avg_score:.4f}")

BEST_K = valid_k[np.argmax(silhouette_avgs)]
print(f"\n  Best K: {BEST_K}")

sil_df = pd.DataFrame({
    "k": valid_k,
    "silhouette_score": silhouette_avgs,
})
sil_df.to_csv(f"{CONDITION}_silhouette_scores.csv", index=False)
print(f"  Saved: {CONDITION}_silhouette_scores.csv")

# Override BEST_K here if needed after inspecting silhouette scores
# FORCE_K = 4
# if FORCE_K is not None:
#     print(f"  [INFO] Overriding BEST_K → using k={FORCE_K}")
#     BEST_K = FORCE_K

print(f"\n[SECTION 5b] K-Means clustering (k={BEST_K})")
kmeans = KMeans(n_clusters=BEST_K, random_state=42, n_init=10)
labels_final = kmeans.fit_predict(X_cluster)

# Rebuild the frame table before merging labels back to trajectories so the
# frame index stays aligned with the original trajectory order.
df_angles = pd.DataFrame(angles_flat, columns=angle_cols)
df_angles.insert(0, "frame", np.arange(n_frames))
frame_ids = df_angles["frame"].to_numpy()

cluster_df = pd.DataFrame({
    "frame": frame_ids[:len(labels_final)],
    "PC1": X_cluster[:, 0],
    "PC2": X_cluster[:, 1],
    "cluster": labels_final,
}) # PC1 + PC2 only

# Guard against silent label/frame misalignment that would otherwise produce
# incorrect downstream cluster assignments.
assert len(cluster_df) == len(labels_final), (
    f"Mismatch: frame_ids={len(frame_ids)} vs labels_final={len(labels_final)}"
)

cluster_df.to_csv(f"{CONDITION}_dPCA_clusters.csv", index=False)

centers_df = pd.DataFrame(
    kmeans.cluster_centers_,
    columns=["PC1", "PC2"],
)
centers_df.to_csv(f"{CONDITION}_dPCA_cluster_centers.csv", index=False)

n_frames = len(cluster_df)

for c in range(BEST_K):
    count = (cluster_df["cluster"] == c).sum()
    pct = count / n_frames * 100
    print(f"  Cluster {c}: {count} frames ({pct:.1f}%)")

print("\n[SECTION 6] Computing Free Energy Landscape FEL")
kB = 8.314e-3  # kJ mol^-1 K^-1
pc1 = cluster_df["PC1"].values
pc2 = cluster_df["PC2"].values

hist, xedges, yedges = np.histogram2d(pc1, pc2, bins=100, density=True)
# Replace zeros to avoid log(0)
hist[hist == 0] = np.nan

# ΔG = −kBT·ln(P/Pmax)
with np.errstate(invalid="ignore"):
    delta_G = -kB * TEMP_K * np.log(hist / np.nanmax(hist))

# Build FEL DataFrame, one row per grid cell
x_centres = 0.5 * (xedges[:-1] + xedges[1:])
y_centres = 0.5 * (yedges[:-1] + yedges[1:])
xx, yy = np.meshgrid(x_centres, y_centres, indexing="ij")

fel_df = pd.DataFrame({
    "PC1": xx.ravel(),
    "PC2": yy.ravel(),
    "delta_G": delta_G.ravel(),
}).dropna() # drop zero-probability bins

fel_df.to_csv(f"{CONDITION}_FEL_dPCA.csv", index=False)

hist_df = pd.DataFrame(hist)
hist_df.to_csv(f"{CONDITION}_dPCA_hist2d.csv", index=False)
print(f"  Saved: {CONDITION}_FEL_dPCA.csv")
print(f"  Saved: {CONDITION}_dPCA_hist2d.csv")

print("\n[SECTION 7] Extracting representative structures")
rep_frames = []

for i, center in enumerate(kmeans.cluster_centers_):
    distances = np.linalg.norm(X_cluster - center, axis=1)
    closest_frame = np.argmin(distances)
    rep_frames.append(closest_frame)
    print(f"  Cluster {i}: representative frame = {closest_frame}")

rep_df = pd.DataFrame({
    "cluster": np.arange(len(rep_frames)),
    "representative_frame": rep_frames,
})
rep_df.to_csv(f"{CONDITION}_representative_frames.csv", index=False)

for i, frame_idx in enumerate(rep_frames):
    u.trajectory[frame_idx]
    filename = f"{CONDITION}_cluster_{i}_rep.pdb"
    protein.write(filename)
    print(f"  Saved: {filename}")

# Merge the reduced coordinates with the original frame information for a final
# master table that connects structural descriptors to cluster labels.
master_df = df_angles.merge(cluster_df, on="frame", how="inner")
master_csv = f"{CONDITION}_dPCA_master.csv"
master_df.to_csv(master_csv, index=False)
print(f" Saved: {master_csv}")
print(f"  Master shape: {master_df.shape}  (raw angles + sincos + PC1/PC2/cluster)")

print("\n[SECTION 8] PCA Convergence (Subspace Overlap)")

def subspace_overlap(A, B):
    return np.linalg.norm(A @ B.T, ord="fro") ** 2 / A.shape[0]


fractions = np.arange(0.05, 1.01, 0.05).tolist()
subspaces = []

for frac in fractions:
    n_cut = int(len(X_scaled) * frac)
    X_sub = X_scaled[:n_cut]
    pca_sub = PCA(n_components=2, svd_solver="full")
    pca_sub.fit(X_sub)
    subspaces.append(pca_sub.components_)

ref = subspaces[-1]
overlaps = [subspace_overlap(s, ref) for s in subspaces]

conv_df = pd.DataFrame({
    "fraction": fractions,
    "subspace_overlap": overlaps,
})
conv_df.to_csv(f"{CONDITION}_PCA_convergence.csv", index=False)

for f, o in zip(fractions, overlaps):
    print(f"  {int(f*100)}% trajectory → overlap = {o:.3f}")

CONVERGENCE_THRESHOLD = 0.90
crossover_frac = None
crossover_overlap = None

for frac, overlap in zip(fractions, overlaps):
    if overlap >= CONVERGENCE_THRESHOLD:
        crossover_frac = frac
        crossover_overlap = overlap
        break

if crossover_frac is not None and crossover_frac <= 0.75:
    verdict = "CONVERGED"
elif crossover_frac is not None and crossover_frac <= 0.90:
    verdict = "CONVERGED"
elif crossover_frac is not None:
    verdict = "MARGINAL"
else:
    verdict = "NOT CONVERGED"

if crossover_frac is not None:
    verdict_str = (
        f"{verdict} "
        f"(overlap ≥ 0.90 first reached at "
        f"{int(crossover_frac*100)}% of trajectory, "
        f"overlap = {crossover_overlap:.3f})"
    )
else:
    verdict_str = f"NOT CONVERGED (overlap never reached 0.90)"

print(f"\n  Convergence verdict : {verdict_str}")
print(f"  Saved: {CONDITION}_PCA_convergence.csv")

print("\n" + "=" * 60)
print("  dPCA PIPELINE COMPLETE")
print("=" * 60)
print(f"  Frames analysed         : {n_frames}")
print(f"  Residues (φ/ψ)          : {n_residues}")
print(f"  Sin/cos features        : {len(sincos_cols)}")
print(f"  PCs retained            : {N_PCS}")
print(f"  Variance PC1+PC2        : {np.sum(pca.explained_variance_ratio_[:2])*100:.2f}%")
print(f"  Cumulative PC1-PC{N_PCS}: {np.sum(pca.explained_variance_ratio_)*100:.2f}%")
print(f"  Best K (silhouette)     : {BEST_K}")
print(f"  Convergence verdict     : {verdict_str}")
print("=" * 60)

sys.stdout.log.close()
sys.stdout = sys.stdout.terminal
