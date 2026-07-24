# time nice -5 snakemake -s 128.smk --use-conda --cores 8 --conda-prefix /software/tmp  --latency-wait 30
#Library
#conda activate snakemake7_TH  
import yaml

#Paths
path_ref_plas = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/11_Sequences/plasmid/SDHB_Ex3_HDR.fasta"
path_bed_plas = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/11_Sequences/plasmid/SDHB_Ex3_HDR.bed"


#sample = ['AN232938_S89', 'AN232941_S90', 'AN232942_S91', 'AN232943_S92', 'AN232944_S93', 'AN232945_S94', 'AN232948_S95', 'AN232949_S96']
sample = [ 'AN232942_S91']


#'AN218971_S62',
rule all:
    input:
        o1 = expand("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/10_Fastqc/{sample}_R{direction}_fastqc/fastqc_data.txt", sample = sample, direction = ["1","2"]),
        o2 = expand ("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/10_Fastqc/{sample}_R{direction}_fastqc.html", sample = sample, direction = ["1","2"]),
        mapping = expand("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/02_Mapping/{sample}.bam", sample=sample),
        allbams = expand("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam", sample = sample),
        allbamindices = expand("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai", sample = sample),
        pileup_quality_more = expand("/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}_quality_more.pileup.no_bq", sample = sample),
        
        
# fastqc rule
rule fastqc1:
    input:
        r = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/01_Fastq/{sample}_R{direction}_001.fastq.gz"
    threads: 1
    priority: 50
    output: 
        o1 = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/10_Fastqc/{sample}_R{direction}_fastqc/fastqc_data.txt",
        o2 = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/10_Fastqc/{sample}_R{direction}_fastqc.html"
    conda:
        "fastqc"
    shell:
        "fastqc "\
        f"-o /mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/10_Fastqc/ "\
        "-t {threads} "\
        "--extract {input.r}"

rule bwa:
    input:
        r1 = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/01_Fastq/{sample}_R1_001.fastq.gz",
        r2 = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/01_Fastq/{sample}_R2_001.fastq.gz",
        ref = path_ref_plas
    params:
        threads = 8
    output: 
        mapping = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/02_Mapping/{sample}.bam"
    conda:
        "bwa"
    shell:
        """
        bwa mem \
        -t {threads} \
        -c 250 \
        -M \
        -R '@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tPU:{wildcards.sample}\\tSM:{wildcards.sample}' \
        {input.ref} \
        {input.r1} \
        {input.r2} \
        | samtools view -b > {output.mapping}
        """

# Mark duplicate alignments
rule MarkDup:
    input:
        bam = lambda wildcards: "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/02_Mapping/{sample}.bam".format(sample=wildcards.sample)
    output:
        sortbam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_sorted.bam",
        fixbam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_fixmate.bam",
        dupbam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
    conda:
        "bwa"
    shell:
        """
        samtools fixmate -m {input.bam} {output.fixbam}; \
        samtools sort -m 4G -o {output.sortbam} {output.fixbam}; \
        samtools markdup {output.sortbam} {output.dupbam}
        """

# Add index files for marked BAM files
rule IndexBam:
    input:
        dupbam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam"
    output:
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai"
    conda:
        "bwa"
    shell:
        "samtools index {input.dupbam}"


rule GeneratePileup:
    input:
        bam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}.pileup"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup})
        samtools mpileup -d 100000 -f {input.ref} {input.bam} > {output.pileup}
        """

rule GeneratePileupNoBQ:
    input:
        bam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup_no = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}.pileup.no_bq"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup_no})
        samtools mpileup -d 200000 -f {input.ref} --ff 0 {input.bam} > {output.pileup_no}
        """

rule GeneratePileupNoBQrange:
    input:
        bam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup_no = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}_range.pileup.no_bq"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup_no})
        samtools mpileup -d 200000 -f {input.ref} --ff 0 -r SDHB_Ex3_HDR:294-465 {input.bam} > {output.pileup_no}
        """

rule GeneratePileupNoBQquality:
    input:
        bam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup_no = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}_quality.pileup.no_bq"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup_no})
        samtools mpileup -d 200000 -f {input.ref} --ff 0 -Q 30 -r SDHB_Ex3_HDR:294-465 {input.bam} > {output.pileup_no}
        """

rule GeneratePileupNoBQquality_more:
    input:
        bam = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam",
        bamindex = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup_no = "/mnt/nct-zfs/CMTD/01_Projects/999-2025/00_Users/lewisnata/05_Variants/{sample}_quality_more.pileup.no_bq"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup_no})
        samtools mpileup -d 200000 -f {input.ref} --ff 0 -Q 30 -q 30 -r SDHB_Ex3_HDR:294-465 {input.bam} > {output.pileup_no}
        """