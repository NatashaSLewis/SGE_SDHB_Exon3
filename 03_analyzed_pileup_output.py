"""
This script compares a reference set of expected variants with variants detected from pileup output. Variants with ALT_COUNT < 10 are filtered out before comparison.

The reference and detected variants are matched using REF, ALT, and POS columns to identify True Positives (TP), False Negatives (FN), and False Positives (FP). 

Theresults are saved to an Excel workbook with separate worksheets for TP, FN, and FP.
"""

import pandas as pd

# Read input Excel files (change filename/path if needed)
df1 = pd.read_excel("path/to/Exon3_PosRefAlt.xlsx",sheet_name='FL',engine='openpyxl') #this file has POS, REF, ALT columns for the variants of interest
df2 = pd.read_csv("path/to/pileup_output.csv") #change the plasmid POS to plas_POS in csv file. # add the hg38 position manually
df2['ALT_COUNT'] = pd.to_numeric(df2['ALT_COUNT'], errors='coerce')
df2 = df2[df2['ALT_COUNT'] >= 10] #low ALT_COUNTs are removed

output_path="path/to/output.xlsx"


merged_df = pd.merge(
    df1,
    df2,
    left_on=['REF', 'ALT', 'POS'],
    right_on=['REF', 'ALT', 'POS'],
    how='inner'
)

FN = pd.merge(
    df1,
    df2,
    left_on=['REF', 'ALT', 'POS'],
    right_on=['REF', 'ALT', 'POS'],
    how='left',
    indicator=True
).query('_merge == "left_only"').drop(columns=['_merge'])

FP = pd.merge(
    df1,
    df2,
    left_on=['REF', 'ALT', 'POS'],
    right_on=['REF', 'ALT', 'POS'],
    how='right',
    indicator=True
).query('_merge == "right_only"').drop(columns=['_merge'])


with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
    merged_df.to_excel(writer, index=False, sheet_name='TP')            # matched rows
    FN.to_excel(writer, index=False, sheet_name='FN')  
    FP.to_excel(writer, index=False, sheet_name='FP') 

print("Merged Excel file saved successfully!")
