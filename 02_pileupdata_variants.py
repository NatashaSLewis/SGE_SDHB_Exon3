
"""
Python 3.9.0

Pileup parsing and SNV allele-frequency calculation.

This script performs the following steps:

1. Reads a tab-separated samtools mpileup file. The pileup with minimum base-quality and mapping-quality thresholds of 30.

2. Parses the read-bases column character by character.

3. Removes or ignores pileup control symbols, including:

   * read-start markers: '^' and the following mapping-quality character
   * read-end markers: '$'
   * insertion and deletion sequence annotations following +<length> or -<length>
   * deletion placeholders: '*' and '#'
   * reference-skip symbols: '<' and '>'

4. Counts nucleotide observations at each genomic position:

   * '.' and ',' are counted as matches to the reference allele
   * explicit A, C, G, T, and N symbols are counted case-insensitively
   * N observations are reported but excluded from SNV depth calculations

5. Records counts of non-SNV pileup events, including insertions, deletions,
   deletion placeholders, reference skips, and unrecognized symbols.

6. Calculates callable nucleotide depth as:

   ```
   CALLABLE_DEPTH = A_COUNT + C_COUNT + G_COUNT + T_COUNT
   ```

   N bases, indel annotations, deletion placeholders, reference skips, and
   unrecognized symbols are excluded from this depth.

7. Reports each non-reference A, C, G, or T allele with at least one counted
   observation.

8. Calculates alternate allele frequency as:

   ```
   AF = ALT_COUNT / CALLABLE_DEPTH
   ```

9. Writes the following fields to a CSV file:

   ```
   CHROM
   POS
   REF
   ALT
   REF_COUNT
   ALT_COUNT
   CALLABLE_DEPTH
   MPILEUP_DEPTH
   N_COUNT
   INSERTION_EVENTS
   DELETION_EVENTS
   DELETION_PLACEHOLDERS
   REFERENCE_SKIPS
   OTHER_SYMBOLS
   AF
   ```

The reported mpileup depth is retained as MPILEUP_DEPTH for comparison, but it
is not used as the allele-frequency denominator.
"""


import pandas as pd
from collections import Counter
from pathlib import Path


pileup_file = Path("path/to/output_from_snakemake.pileup" ) 
output_csv = Path("path/to/output.csv")

data = []


def parse_pileup_bases(bases: str, ref: str) -> tuple[Counter, Counter]:
    """
    Parse the read-bases column from samtools mpileup.

    Returns
    -------
    base_counts:
        Counts of callable A/C/G/T/N bases at the current reference position.

    event_counts:
        Counts of insertions, deletions, deletion placeholders,
        reference skips, and other ignored symbols.

    Important
    ---------
    Bases written after +<length> or -<length> describe an indel sequence.
    They must be skipped and must not be counted as SNVs.
    """
    base_counts = Counter()
    event_counts = Counter()

    i = 0

    while i < len(bases):
        symbol = bases[i]

        # Start of a read segment:
        # "^" is followed by one mapping-quality character.
        if symbol == "^":
            if i + 1 >= len(bases):
                raise ValueError("Malformed pileup: '^' lacks mapping-quality character")
            i += 2
            continue

        # End of a read segment
        if symbol == "$":
            i += 1
            continue

        # Reference match:
        # "." = forward strand, "," = reverse strand
        if symbol in ".,":
            base_counts[ref] += 1
            i += 1
            continue

        # Explicit nucleotide mismatch
        if symbol.upper() in {"A", "C", "G", "T", "N"}:
            base_counts[symbol.upper()] += 1
            i += 1
            continue

        # Insertion or deletion annotation:
        # Example: +3ACT or -11TCTTCGGAAGG
        if symbol in "+-":
            event_type = "INSERTION" if symbol == "+" else "DELETION"
            i += 1

            number_start = i
            while i < len(bases) and bases[i].isdigit():
                i += 1

            if number_start == i:
                raise ValueError(
                    f"Malformed pileup: '{symbol}' is not followed by an indel length"
                )

            indel_length = int(bases[number_start:i])
            indel_end = i + indel_length

            if indel_end > len(bases):
                raise ValueError(
                    f"Malformed pileup: declared indel length {indel_length} "
                    "extends beyond the read-bases string"
                )

            event_counts[event_type] += 1
            i = indel_end
            continue

        # Placeholder for a base deleted in a preceding pileup column.
        # "*" is normally forward strand; "#" may occur for reverse strand
        # in newer samtools output.
        if symbol in "*#":
            event_counts["DELETION_PLACEHOLDER"] += 1
            i += 1
            continue

        # Reference skip, commonly associated with spliced alignments
        if symbol in "<>":
            event_counts["REFERENCE_SKIP"] += 1
            i += 1
            continue

        # Unexpected symbols are recorded rather than silently treated as bases
        event_counts["OTHER_SYMBOL"] += 1
        i += 1

    return base_counts, event_counts


def parse_pileup_line(line: str, line_number: int):
    cols = line.rstrip("\n").split("\t")

    if len(cols) < 5:
        raise ValueError(
            f"Line {line_number}: expected at least 5 tab-separated columns, "
            f"found {len(cols)}"
        )

    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[2].upper()
    mpileup_depth = int(cols[3])
    bases = cols[4]

    if ref not in {"A", "C", "G", "T", "N"}:
        raise ValueError(
            f"Line {line_number}: unsupported reference base '{ref}'"
        )

    base_counts, event_counts = parse_pileup_bases(bases, ref)

    # Only A/C/G/T observations contribute to the SNV denominator.
    # N bases, deletions, reference skips, and indel annotations are excluded.
    callable_depth = sum(base_counts[base] for base in "ACGT")

    records = []

    for alt in "ACGT":
        if alt == ref:
            continue

        alt_count = base_counts[alt]

        if alt_count == 0:
            continue

        af = alt_count / callable_depth if callable_depth > 0 else 0.0

        records.append(
            {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "REF_COUNT": base_counts[ref],
                "ALT_COUNT": alt_count,
                "CALLABLE_DEPTH": callable_depth,
                "MPILEUP_DEPTH": mpileup_depth,
                "N_COUNT": base_counts["N"],
                "INSERTION_EVENTS": event_counts["INSERTION"],
                "DELETION_EVENTS": event_counts["DELETION"],
                "DELETION_PLACEHOLDERS": event_counts["DELETION_PLACEHOLDER"],
                "REFERENCE_SKIPS": event_counts["REFERENCE_SKIP"],
                "OTHER_SYMBOLS": event_counts["OTHER_SYMBOL"],
                "AF": af,
            }
        )

    return records


if not pileup_file.exists():
    raise FileNotFoundError(f"Pileup file not found: {pileup_file}")

with pileup_file.open("r", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, start=1):
        if not line.strip():
            continue

        try:
            data.extend(parse_pileup_line(line, line_number))
        except ValueError as error:
            raise ValueError(
                f"Failed to parse {pileup_file} at line {line_number}: {error}"
            ) from error


columns = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "REF_COUNT",
    "ALT_COUNT",
    "CALLABLE_DEPTH",
    "MPILEUP_DEPTH",
    "N_COUNT",
    "INSERTION_EVENTS",
    "DELETION_EVENTS",
    "DELETION_PLACEHOLDERS",
    "REFERENCE_SKIPS",
    "OTHER_SYMBOLS",
    "AF",
]

df = pd.DataFrame(data, columns=columns)

output_csv.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_csv, index=False)

print(f"Saved {len(df):,} SNV records to {output_csv}")
