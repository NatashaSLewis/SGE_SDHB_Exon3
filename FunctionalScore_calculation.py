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
