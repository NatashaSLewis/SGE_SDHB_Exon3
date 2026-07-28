"""
Pileup parsing and basic SNV allele-frequency calculation.

This script performs the following steps:

1. Reads pileup output generated with minimum base-quality and mapping-quality
   thresholds of 30.
2. Removes pileup read-start and read-end markers.
3. Counts reference matches represented by '.' and ',' and counts explicit
   A, C, G, T, and N observations at each genomic position.
4. Calculates alternate allele frequency using:
       AF = ALT_COUNT / REPORTED_PILEUP_DEPTH
5. Reports non-reference alleles with at least one counted observation.
6. Writes chromosome, position, reference allele, alternate allele,
   alternate-allele count, counted nucleotide depth, and allele frequency
   to a CSV file.

The pileup was generated using snakemake output with base quality below 30 and reads with mapping quality below 30 were
excluded during pileup generation.

"""
import pandas as pd
from collections import Counter
import re

pileup_file = "path/to/output_from_snakemake.pileup" #changes
output_csv = "path/to/output.csv"  #change


data = []

def parse_pileup_line(line):
    cols = line.strip().split('\t')
    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[2].upper()
    depth = int(cols[3])
    bases = cols[4]

    # remove read start (^.) and end ($) markers
    bases = re.sub(r'\^.', '', bases)
    bases = bases.replace('$', '')

    # map pileup symbols to canonical bases
    counts = Counter()
    for b in bases:
        b = b.upper()
        if b == '.' or b == ',':
            counts[ref] += 1
        elif b in ['A','C','G','T','N']:
            counts[b] += 1
        else:
            pass  # ignore other symbols

    # compute AF for each alt base
    af_dict = {}
    for base in ['A','C','G','T']:
        if base != ref:
            af_dict[base] = counts[base] / depth if depth > 0 else 0

    return chrom, pos, ref, counts, af_dict

with open(pileup_file) as f:
    for line in f:
        chrom, pos, ref, counts, af_dict = parse_pileup_line(line)
        for alt, af in af_dict.items():
            if af > 0:
                data.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "ALT_COUNT": counts[alt],
                    "DEPTH": sum(counts.values()),
                    "AF": af
                })

# Convert to DataFrame
df = pd.DataFrame(data)


# Save as CSV
df.to_csv(output_csv, index=False)

# Optional: save as Excel
#df.to_excel("sample_vaf_table.xlsx", index=False)
