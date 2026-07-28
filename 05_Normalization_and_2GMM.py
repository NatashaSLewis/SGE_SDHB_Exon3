"""
Two-component Gaussian mixture model for functional variant classification.

This script performs the following steps:

1. Extracts functional scores from three biological replicate sets.
2. Robust-normalizes each replicate independently using:
       z_robust = (x - median) / IQR
   This places replicate sets on comparable scales while reducing sensitivity
   to extreme values.
3. Fits a two-component Gaussian mixture model (GMM) jointly to the three
   normalized scores. Each variant is therefore represented as a point in
   three-dimensional functional-score space rather than by its mean score.
4. Designates the component with the lower centroid across the three
   normalized dimensions as the functionally-abnormal-like component.
5. Uses GMM posterior component-membership probabilities to classify variants:
       p_abnormal >= 0.95  -> functionally-abnormal
       p_abnormal <= 0.05  -> functionally-normal
       otherwise           -> uncertain
6. Writes posterior probabilities, classifications, normalized scores, and
   component labels back to the original dataset.

The reported probabilities represent membership in the fitted GMM components;
they are not direct probabilities of clinical functionally_abnormality.

The model is unsupervised and assumes that:
- the dataset contains two underlying functional distributions;
- lower functional scores correspond to reduced function;
- the three replicate scores are sufficiently comparable after normalization;
- each component can be reasonably approximated by a multivariate Gaussian.

Model fit and stability should be evaluated using component centroids,
covariance matrices, convergence status, repeated initializations, and
sensitivity analyses using alternative covariance structures.
"""

import numpy as np
from sklearn.mixture import GaussianMixture
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# Input
# ----------------------------

'''
file_path = "path/to/file/file_with_rawFunctional_scores.xlsx"
output_path = "path/to/file/file_with_normalizedscores_and_probabilityClass.xlsx"
'''
file_path = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/FL/TableS3_MAVE_mainTableFinal.xlsx"
output_path = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/FL/temp7.xlsx"


df = pd.read_excel(file_path, sheet_name="Sheet1", engine="openpyxl")
df.columns = df.columns.astype(str).str.strip()

# ----------------------------
# Columns to use
# ----------------------------
#cols = ["Functional_score_haploid_set1",  "Functional_score_haploid_set2", "Functional_score_haploid_set3"] # Raw scores haploid
cols = ["Functional_score_hemizygous_set1",  "Functional_score_hemizygous_set2", "Functional_score_hemizygous_set3"] # Raw scores hemizygous


# check columns exist
missing = [c for c in cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# keep only rows where all 3 scores exist
df_sub = df[cols].dropna().copy()

# remember original row index so results can be written back correctly
valid_idx = df_sub.index

# ----------------------------
# Robust normalization
# z_robust = (x - median) / IQR
# ----------------------------
X_norm = df_sub.copy()

norm_params = {}
for c in cols:
    med = X_norm[c].median()
    iqr = X_norm[c].quantile(0.75) - X_norm[c].quantile(0.25)

    if iqr == 0:
        raise ValueError(f"IQR is zero for column '{c}', cannot robust-normalize.")

    X_norm[c] = (X_norm[c] - med) / iqr
    norm_params[c] = {"median": med, "iqr": iqr}

# convert to numpy for GMM
X = X_norm.to_numpy()

# ----------------------------
# Fit 2-component GMM in 3D
# ----------------------------
gmm = GaussianMixture(n_components=2, random_state=42, covariance_type="full")
gmm.fit(X)

# component with LOWER mean across normalized dimensions = more functionally_abnormal-like
component_mean_score = gmm.means_.mean(axis=1)
functionally_abnormal_k = int(np.argmin(component_mean_score))
functionally_normal_k = 1 - functionally_abnormal_k

# posterior probabilities
proba = gmm.predict_proba(X)
p_functionally_abnormal = proba[:, functionally_abnormal_k]
p_functionally_normal = proba[:, functionally_normal_k]



# ----------------------------
# Classification thresholds (adjust threshold accordingly)
# ----------------------------
p_hi = 0.95   # confident functionally_abnormal-like
p_lo = 0.05   # confident functionally_normal-like

call = np.where(
    p_functionally_abnormal >= p_hi, "functionally_abnormal_like",
    np.where(p_functionally_abnormal <= p_lo, "functionally_normal_like", "uncertain")
)

# ----------------------------
# Write back to original dataframe
# ----------------------------
df["p_functionally_abnormal"] = np.nan
df["p_functionally_normal"] = np.nan
df["call"] = np.nan

df.loc[valid_idx, "p_functionally_abnormal"] = p_functionally_abnormal
df.loc[valid_idx, "p_functionally_normal"] = p_functionally_normal
df.loc[valid_idx, "call"] = call

# optional: save normalized values too
for c in cols:
    df[c + "_robust_norm"] = np.nan
    df.loc[valid_idx, c + "_robust_norm"] = X_norm[c].values

# optional: most likely cluster label
raw_cluster = gmm.predict(X)

# Remap to biological labels:
# 1 = functionally abnormal
# 0 = functionally normal
gmm_cluster = np.where(raw_cluster == functionally_abnormal_k, 1, 0)

df["gmm_cluster"] = np.nan
df.loc[valid_idx, "gmm_cluster"] = gmm_cluster

# ----------------------------
# Print summary
# ----------------------------
print("Columns used:", cols)
print("\nRobust normalization parameters:")
for c in cols:
    print(f"{c}: median={norm_params[c]['median']:.4f}, IQR={norm_params[c]['iqr']:.4f}")

print("\nGMM component means (normalized space):")
print(gmm.means_)

print("\nAverage component means:")
print(component_mean_score)


print("\nCall counts:")
print(pd.Series(call).value_counts(dropna=False))

print("\nPreview:")
print(df[[*cols, "p_functionally_abnormal", "p_functionally_normal", "call"]].head(20))

# ----------------------------
# Save to Excel
# ----------------------------
df.to_excel(output_path, index=False, engine="openpyxl")
print(f"\nSaved output to: {output_path}")