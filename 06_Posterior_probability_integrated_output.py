"""
Integration of haploid and hemizygous functional probabilities.

The input file contains GMM-derived posterior probabilities of membership in
the functionally abnormal component for the haploid and hemizygous assays.

For variants with results in both assays, the two probabilities are combined
using an odds-based formula:

    p_combined = (p_haploid * p_hemizygous) /
                 [(p_haploid * p_hemizygous) +
                  ((1 - p_haploid) * (1 - p_hemizygous))]

This formulation strengthens concordant evidence:
    - two high probabilities produce a higher combined probability;
    - two low probabilities produce a lower combined probability;
    - opposing probabilities move the result toward 0.5.

Probabilities are clipped slightly away from 0 and 1 to prevent numerical
problems caused by exact boundary values.

A confidence rule is then applied. A result is considered high confidence when:
    1. at least one assay produces a probability >= 0.95 or <= 0.05; and
    2. the difference between the assay probabilities is < 0.9.

Combined probabilities that do not satisfy this rule are conservatively
shrunk toward 0.5:

    p_adjusted = 0.5 + (p_combined - 0.5) * 0.25

Thus, only 25% of the original deviation from 0.5 is retained for results
with insufficient or strongly discordant evidence. This adjustment is a
prespecified conservative decision rule rather than a Bayesian update.

The resulting probabilities summarize integrated functional evidence and
should not be interpreted directly as clinical probabilities of pathogenicity.
This value is subsequently used for functional classification.
"""
# ----------------------------
# Input and output files
# ----------------------------
import pandas as pd
import numpy as np


# ==========================
# Files
# ==========================

df = pd.read_excel('path/to/haploid_and_hemizygous_probability.xlsx', 
                   sheet_name="Sheet1",engine="openpyxl")

output_file = "path/to/output.xlsx"



# ==========================
# Columns
# ==========================

p_haploid_col = "p_damaging_haploid"
p_hemi_col = "p_damaging_hemizygous"

# ==========================
# Clean probabilities
# ==========================

df[p_haploid_col] = pd.to_numeric(
    df[p_haploid_col],
    errors="coerce"
)

df[p_hemi_col] = pd.to_numeric(
    df[p_hemi_col],
    errors="coerce"
)


# remove missing
df = df.dropna(
    subset=[p_haploid_col, p_hemi_col]
)


# avoid division problems
eps = 1e-10

p1 = df[p_haploid_col].clip(eps, 1-eps)
p2 = df[p_hemi_col].clip(eps, 1-eps)


# ==========================
# Calculate combined probability
# ==========================

numerator = p1 * p2

denominator = (
    (p1 * p2) +
    ((1-p1) * (1-p2))
)

df["new_combined_p"] = numerator / denominator


# ==========================
# Apply confidence adjustment
# ==========================

confidence_rule = (
    (
        (p1 >= 0.95) |
        (p2 >= 0.95) |
        (p1 <= 0.05) |
        (p2 <= 0.05)
    )
    &
     (df["delta_nosign"] < 0.9)

)


df["new_combined_p_adjusted"] = df["new_combined_p"]



# push uncertain variants toward 0.5
df.loc[
    ~confidence_rule,
    "new_combined_p_adjusted"
] = 0.5 + (
    df.loc[
        ~confidence_rule,
        "new_combined_p"
    ] - 0.5
) * 0.25


# ==========================
# Optional: label reason
# ==========================

df["confidence_status"] = np.where(
    confidence_rule,
    "high_confidence",
    "adjusted_to_0.5"
)


# ==========================
# Save
# ==========================

df.to_excel(
    output_file,
    index=False
)

print("Saved:", output_file)
print(
    "Adjusted variants:",
    (~confidence_rule).sum()
)