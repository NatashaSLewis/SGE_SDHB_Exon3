"""
Functional-score calculation from variant allele frequencies at two timepoints.

This script performs the following steps:

1. Reads variant data from the `TP` worksheets of two Excel files representing
   an earlier and a later experimental timepoint.
2. Merges the datasets using genomic position, reference allele, and alternate
   allele:
   POS, REF, and ALT
3. Uses an outer merge so that variants present at only one timepoint are
   retained.
4. Replaces missing or very low allele frequencies with a pseudocount:
   epsilon = 6 × 10⁻⁵
   This value corresponds approximately to 10 alternate reads at the expected
   sequencing depth and prevents division by zero or undefined logarithms.
5. Calculates the functional score for each variant using:
   Functional_score = log2(AF_later / AF_earlier)
6. Writes the variant identifiers and calculated functional scores to an Excel
   file.

A positive functional score indicates an increase in variant allele frequency
between the two timepoints, whereas a negative score indicates a decrease.

"""
import pandas as pd
import numpy as np
import statsmodels.api as sm

# -----------------------------
# ---- Load Excel files ----
# -----------------------------

day4_df = pd.read_excel(
    "path/to/NGS_folder/timepoint1.xlsx",
    sheet_name='TP',
    engine='openpyxl'
)

day6_df = pd.read_excel(
    "path/to/NGS_folder/timepoint2.xlsx",
    sheet_name='TP',
    engine='openpyxl'
)

# -----------------------------
# ---- Step 3a: Functional scores (day6 live, alternative method) ----
# -----------------------------
def calculate_functional_scores_alt(day6_df, day4_df):
    merged = pd.merge(day4_df, day6_df, on=['POS','ALT','REF'], how='outer', suffixes=('_day4','_day6'))
    eps = 6e-5 #this is around 10 ALT reads
   
    merged["AF_day6"] = merged["AF_day6"].fillna(eps).clip(lower=eps)
    merged["AF_day4"] = merged["AF_day4"].fillna(eps).clip(lower=eps)
    
    merged['Functional_score'] = np.log2(merged['AF_day6'] / merged['AF_day4'])
    return merged[['POS','ALT','REF','Functional_score']]
       
# -----------------------------
# ---- Step 4: Run workflow ----
# -----------------------------

final_df = calculate_functional_scores_alt(day6_df, day4_df)


# Save
output_path = "/path/to/output_file/output.xlsx"
final_df.to_excel(output_path, index=False, sheet_name='Scores')

# -----------------------------
