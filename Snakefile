import pandas as pd

# Загрузка конфига и образцов
configfile: "config/config.yaml"
samples = pd.read_table("config/samples.tsv").set_index("sample", drop=False)

# Правило ALL
rule all:
    input:
        # QC
        expand("{output_dir}/qc/fastqc/{sample}_R1_fastqc.html", sample=samples.index),
        expand("{output_dir}/qc/trimmed/{sample}_R1.trimmed.fastq", sample=samples.index),
        # Сборка
        "{output_dir}/assembly/final_assembly.fasta",
        # Биннинг
        "{output_dir}/binning/metabat2/bins",
        # Аннотация
        expand("{output_dir}/annotation/{sample}/prokka.gff", sample=samples.index),
        # BUSCO
        "{output_dir}/busco/short_summary.txt"

# --- Контроль качества ---
rule fastqc:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "R1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "R2"]
    output:
        html1 = "{output_dir}/qc/fastqc/{sample}_R1_fastqc.html",
        html2 = "{output_dir}/qc/fastqc/{sample}_R2_fastqc.html"
    conda:
        "envs/qc.yaml"
    shell:
        "fastqc {input.r1} {input.r2} -o {output_dir}/qc/fastqc/"

rule trim_adapters:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "R1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "R2"]
    output:
        r1 = "{output_dir}/qc/trimmed/{sample}_R1.trimmed.fastq",
        r2 = "{output_dir}/qc/trimmed/{sample}_R2.trimmed.fastq"
    conda:
        "envs/qc.yaml"
    shell:
        "trimmomatic PE -threads {config[threads]} "
        "{input.r1} {input.r2} "
        "{output.r1} {output_dir}/qc/trimmed/{wildcards.sample}_R1.unpaired.fastq "
        "{output.r2} {output_dir}/qc/trimmed/{wildcards.sample}_R2.unpaired.fastq "
        "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50"

# --- Сборка генома ---
rule assemble_megahit:
    input:
        r1 = expand("{output_dir}/qc/trimmed/{sample}_R1.trimmed.fastq", sample=samples.index),
        r2 = expand("{output_dir}/qc/trimmed/{sample}_R2.trimmed.fastq", sample=samples.index)
    output:
        "{output_dir}/assembly/final_assembly.fasta"
    conda:
        "envs/assembly.yaml"
    shell:
        "megahit -1 {input.r1} -2 {input.r2} "
        "-o {output_dir}/assembly/megahit "
        "-t {config[threads]} "
        "--min-contig-len 1000 && "
        "cp {output_dir}/assembly/megahit/final.contigs.fa {output}"

# --- Биннинг ---
rule binning_metabat2:
    input:
        assembly = "{output_dir}/assembly/final_assembly.fasta",
        bams = expand("{output_dir}/mapping/{sample}.sorted.bam", sample=samples.index)
    output:
        directory("{output_dir}/binning/metabat2/bins")
    conda:
        "envs/binning.yaml"
    shell:
        "jgi_summarize_bam_contig_depths "
        "--outputDepth {output_dir}/binning/metabat2/depth.txt "
        "{input.bams} && "
        "metabat2 -i {input.assembly} "
        "-a {output_dir}/binning/metabat2/depth.txt "
        "-o {output_dir}/binning/metabat2/bins/bin "
        "-t {config[threads]}"

# --- Аннотация ---
rule annotate_prokka:
    input:
        fasta = "{output_dir}/binning/metabat2/bins/bin.{bin_id}.fa"
    output:
        "{output_dir}/annotation/{bin_id}/prokka.gff"
    conda:
        "envs/annotation.yaml"
    shell:
        "prokka --outdir {output_dir}/annotation/{wildcards.bin_id} "
        "--prefix prokka {input.fasta}"

# --- Оценка качества ---
rule run_busco:
    input:
        "{output_dir}/assembly/final_assembly.fasta"
    output:
        "{output_dir}/busco/short_summary.txt"
    conda:
        "envs/busco.yaml"
    shell:
        "busco -i {input} -l {config[busco_db]} "
        "-o {output_dir}/busco "
        "-m genome --cpu {config[threads]} && "
        "cp {output_dir}/busco/run_{config[busco_db]}/short_summary.txt {output}"