# conda activate snakemake 
# time nice -5 snakemake -s code_name.smk --use-conda --cores 8 --conda-prefix /software/tmp  --latency-wait 30

#Paths
path_ref_plas = "11_Sequences/plasmid/SDHB_Ex3_HDR.fasta"
path_bed_plas = "11_Sequences/plasmid/SDHB_Ex3_HDR.bed"


sample = [ 'Name_of_sample_file_string']

rule all:
    input:
        o1 = expand("10_Fastqc/{sample}_R{direction}_fastqc/fastqc_data.txt", sample = sample, direction = ["1","2"]),
        o2 = expand ("10_Fastqc/{sample}_R{direction}_fastqc.html", sample = sample, direction = ["1","2"]),
        mapping = expand("02_Mapping/{sample}.bam", sample=sample),
        allbams = expand("03_MarkedDups/{sample}_marked.bam", sample = sample),
        allbamindices = expand("03_MarkedDups/{sample}_marked.bam.bai", sample = sample),
        pileup_quality_more = expand("05_Variants/{sample}_quality_more.pileup.no_bq", sample = sample),
        
        
# fastqc rule
rule fastqc1:
    input:
        r = "01_Fastq/{sample}_R{direction}_001.fastq.gz"
    threads: 1
    priority: 50
    output: 
        o1 = "10_Fastqc/{sample}_R{direction}_fastqc/fastqc_data.txt",
        o2 = "10_Fastqc/{sample}_R{direction}_fastqc.html"
    conda:
        "fastqc"
    shell:
        "fastqc "\
        f"-o 10_Fastqc/ "\
        "-t {threads} "\
        "--extract {input.r}"

rule bwa:
    input:
        r1 = "01_Fastq/{sample}_R1_001.fastq.gz",
        r2 = "01_Fastq/{sample}_R2_001.fastq.gz",
        ref = path_ref_plas
    params:
        threads = 8
    output: 
        mapping = "02_Mapping/{sample}.bam"
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
        bam = lambda wildcards: "02_Mapping/{sample}.bam".format(sample=wildcards.sample)
    output:
        sortbam = "03_MarkedDups/{sample}_sorted.bam",
        fixbam = "03_MarkedDups/{sample}_fixmate.bam",
        dupbam = "03_MarkedDups/{sample}_marked.bam",
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
        dupbam = "03_MarkedDups/{sample}_marked.bam"
    output:
        bamindex = "03_MarkedDups/{sample}_marked.bam.bai"
    conda:
        "bwa"
    shell:
        "samtools index {input.dupbam}"



rule GeneratePileupNoBQquality_more:
    input:
        bam = "03_MarkedDups/{sample}_marked.bam",
        bamindex = "03_MarkedDups/{sample}_marked.bam.bai",
        ref = path_ref_plas
    output:
        pileup_no = "05_Variants/{sample}_quality_more.pileup.no_bq"
    conda: "samtools_115"  
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.pileup_no})
        samtools mpileup -d 200000 -f {input.ref} --ff 0 -Q 30 -q 30 -r SDHB_Ex3_HDR:294-465 {input.bam} > {output.pileup_no}
        """
