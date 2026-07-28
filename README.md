The analysis pipeline consists of 6 separate codes to be run in sequential steps
Step1: Snakemake pipeline guide 
Step2: pileup file analyses 
Step3: True positives from pileup 
Step4: Calculate Functional scores 
Step5: Normalize functional scores, 2GMM and probability of pathogenicity
Step6: Integration of ploidy data

Each code has explanation of workflow. This docucument is for specifying the sequence of analysis.

SNAKEMAKE:
Step1: Snakemake pipeline guide:
## Requirements
-Snakemake 7.24.2

Paired-end sequencing quality control, plasmid alignment, duplicate marking,
and pileup generation using Snakemake.

Before running the workflow:
Activate the Snakemake environment:
conda activate snakemake
Place paired-end FASTQ files in the following directory:
01_Fastq/                        ###adjust according to your file format
Files must follow this naming convention:
<sample>_R1_001.fastq.gz
<sample>_R2_001.fastq.gz
Add each sample name to the sample list without the read suffixes:
sample = ["Sample1","Sample2"]

Provide the plasmid reference files:
11_Sequences/plasmid/SDHB_Ex3_HDR.fasta       ###adjust according to your file format
11_Sequences/plasmid/SDHB_Ex3_HDR.bed        ###adjust according to your file format
    
The FASTA reference must be indexed for BWA before running the workflow:
bwa index 11_Sequences/plasmid/SDHB_Ex3_HDR.fasta
A FASTA index may also be required by samtools:
samtools faidx 11_Sequences/plasmid/SDHB_Ex3_HDR.fasta

Ensure that the Conda environment definitions used by the rules are
available:
1. fastqc
2. bwa
3. samtools
    
Run the workflow using:
time nice -n 5 snakemake \
    -s code_name.smk \  ###adjust according to the snakemake code saved
    --use-conda \
    --cores 8 \
    --conda-prefix /software/tmp \
    --latency-wait 30
    
Principal output files are written to: #adjust according to your file path
10_Fastqc/       FastQC reports
02_Mapping/      Initial BAM alignments
03_MarkedDups/   Sorted, duplicate-marked, and indexed BAM files
05_Variants/     Quality-filtered pileup files

Important considerations:
•	The region name SDHB_Ex3_HDR must exactly match the sequence identifier in
the reference FASTA.
•	The requested region SDHB_Ex3_HDR:294-465 uses one-based genomic
coordinates.
•	Duplicate reads are marked in the BAM but are retained in the pileup.
•	Because --ff 0 removes the default flag-exclusion mask, secondary and
QC-failed alignments may also be included unless they were removed earlier.
•	The maximum pileup depth is set to 200,000 reads per position (to override default setting)


PYTHON:
for all following steps 
## Requirements
- Python 3.6
  
Step2: pileup file analyses 
Relevant packages before running code
import pandas as pd
from collections import Counter
import re

command:
python 02_pileupdata_variants.py



Step3: True positives from pileup 
Relevant packages before running code
import pandas as pd

Additional requirement: 
xlsx with POS, REF, ALT of your SNV library 

command:
python 03_analyzed_pileup_output.py



Step4: Calculate Functional scores 
Relevant packages before running code
import pandas as pd
import numpy as np
import statsmodels.api as sm

Additional requirement: 
separate xlsx for timepoint 1 and 2 

command:
python FunctionalScore_calculation.py



Step5: Normalize functional scores, 2GMM and probability of pathogenicity
Relevant packages before running code
import numpy as np
from sklearn.mixture import GaussianMixture
import pandas as pd
import matplotlib.pyplot as plt

command:
python 05_Normalization_and_2GMM.py



Step6: Integration of ploidy scores
Relevant packages before running code
import numpy as np
import pandas as pd

command:
python 06_Posterior_probability_integrated_output.py



