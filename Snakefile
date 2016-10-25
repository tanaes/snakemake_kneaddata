import tempfile
import os
import glob


# configfile: "config.yaml"

TMP_DIR_ROOT = config["TMP_DIR_ROOT"]
RUN = config["RUN"]
SAMPLES_PE = config["samples_pe"] if "samples_pe" in config else []
SAMPLES_SE = config["samples_se"] if "samples_se" in config else []

# Path to programs (or element on path)
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]


localrules: raw_make_links_pe, raw_make_links_se, multiQC_run, multiQC_all


#### Top-level rules: rules to execute a subset of the pipeline

rule all:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    """
    input:
        expand( # fastqc zip and html for raw PE data
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            run = RUN,
            end = "R1 R2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_SE,
            run = RUN,
            end = "SE".split(),
            extension = "zip html".split()
        ),
        expand( # trimmomatic output for PE data
            "data/{sample}/{run}/trimmed/{sample}_{end}.trimmed.fq.gz",
            sample = SAMPLES_PE,
            run = RUN,
            end = "R1 R2 up".split()
        ) + expand( # fastqc zip and html for raw SE data
            "data/{sample}/{run}/trimmed/{sample}_{end}.trimmed.fq.gz",
            sample = SAMPLES_SE,
            run = RUN,
            end = "SE".split()
        ),
        expand(
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand(
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        ),
        expand(
            "data/multiQC/{run}/multiqc_report.html",
            run = RUN
        ),
        "data/multiQC/all/multiqc_report.html"


rule raw_fastqc:
    """
    Rule to do just QC on raw reads:
        - raw_fastqc
    """
    input:
        expand(
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand(
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        )


rule raw_make_links_pe:
    input:
        forward = lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        reverse = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    threads:
        1
    log:
        "logs/{run}/raw/make_links_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/raw/make_links_pe_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.forward}) {output.forward} 2> {log}
        ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}
        """ 


rule raw_make_links_se:
    input:
        single = lambda wildcards: config["samples_se"][wildcards.sample]["forward"],
    output:
        single = "data/{sample}/{run}/raw/{sample}_SE.fq.gz"
    threads:
        1
    log:
        "logs/{run}/raw/make_links_se_{sample}.log"
    benchmark:
        "benchmarks/{run}/raw/make_links_se_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.single}) {output.single} 2>  {log}
        """ 


rule raw_fastqc_sample:
    input:
        fastq = "data/{sample}/{run}/raw/{sample}_{end}.fq.gz"
    output:
        html = "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip"
    threads:
        1
    params:
        html = "data/{sample}/{run}/raw/{sample}_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/raw/{sample}_{end}_fastqc.zip"
    log:
        "logs/{run}/raw/fastqc_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/raw/fastqc_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_raw \
            {input.fastq} \
        2> {log} 1>&2
        """


rule qc_kneaddata_pe:
    """
    Run kneaddata on paired end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    """
    input:
        forward = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        reverse = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    output:
        paired_f  = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_R1.fq.gz",
        paired_r  = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_R2.fq.gz",
        unpaired_f = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_unmatched_R1.fq.gz",
        unpaired_r = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_unmatched_R2.fq.gz",
        all_f = temp("data/{sample}/{run}/kneaddata/{sample}_kneaddata_R1.fq.gz"),
        all_r = temp("data/{sample}/{run}/kneaddata/{sample}_kneaddata_R2.fq.gz")
    params:
        db       = config["kneaddata_db"],
        output_prefix = "{sample}_kneaddata",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/{run}/qc/kneaddata_pe_{sample}.json"
    log:
        "logs/{run}/qc/kneaddata_pe_{sample}.log" 
    threads:
        8
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  kneaddata \
                    --input {input.forward} \
                    --input {input.reverse} \
                    --output %s \
                    --output-prefix {params.output_prefix} \
                    --reference_db {params.db} \
                    --quality-scores {params.phred} \
                    --threads {threads} \
                    --trimmomatic {trimmomatic} \
                    --trimmomatic-options "--ILLUMINACLIP:{params.adaptor}:2:30:10 {params.trimmomatic_params}" \
                  2> {log}

                  {gzip} %s/*
                  
                  scp %s/{sample}_kneaddata_paired_1.fastq {output.paired_f}
                  scp %s/{sample}_kneaddata_paired_2.fastq {output.paired_r}
                  scp %s/{sample}_kneaddata_unmatched_1.fastq {output.unpaired_f}
                  scp %s/{sample}_kneaddata_unmatched_2.fastq {output.unpaired_r}
                  zcat {output.paired_f} {output.unpaired_f} > {output.all_f}
                  zcat {output.paired_r} {output.unpaired_r} > {output.all_r}
                  """ % (temp_dir, temp_dir, temp_dir, temp_dir, temp_dir,
                         temp_dir))



rule qc_kneaddata_se:
    """
    Run kneaddata on single end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    """
    input:
        single = "data/{sample}/{run}/raw/{sample}_SE.fq.gz",
    output:
        single  = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_SE.fq.gz"
    params:
        db       = config["kneaddata_db"],
        output_prefix = "{sample}_kneaddata",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/{run}/qc/kneaddata_pe_{sample}.json"
    log:
        "logs/{run}/qc/kneaddata_pe_{sample}.log" 
    threads:
        8
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  kneaddata \
                    --input {input.forward} \
                    --output %s \
                    --output-prefix {params.output_prefix} \
                    --reference_db {params.db} \
                    --quality-scores {params.phred} \
                    --threads {threads} \
                    --trimmomatic {trimmomatic} \
                    --trimmomatic-options "--ILLUMINACLIP:{params.adaptor}:2:30:10 {params.trimmomatic_params}" \
                  2> {log}

                  {gzip} %s/*
                  
                  scp %s/{sample}_kneaddata.fastq {output.single}
                  """ % (temp_dir, temp_dir, temp_dir))




rule qc_fastqc:
    """
    Do FASTQC reports
    One thread per fastq.gz file
    """
    input:
        fastq = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_{end}.fq.gz"
    output:
        html = "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.zip"
    threads:
        1
    log:
        "logs/{run}/qc/fastqc_kneaddata_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/qc/fastqc_kneaddata_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_trimmed \
            {input.fastq} 
        2> {log} 1>&2
        """


rule multiQC_run:
    input: 
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split())
    output:
        "data/multiQC/{run}/multiqc_report.html"
    threads:
        1
    log:
        "logs/{run}/qc/multiqc.log"
    benchmark:
        "benchmarks/{run}/qc/multiqc.json"
    shell:
        """
        set +u; source activate multiqc; set -u
        multiqc -f -o data/multiQC/{wildcards.run} data/*/{wildcards.run} 2> {log} 1>&2
        """


rule multiQC_all:
    input:
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split())
    output:
        "data/multiQC/all/multiqc_report.html"
    threads:
        1
    log:
        "logs/multiqc_all.log"
    benchmark:
        "benchmarks/multiqc_all.json"
    shell:
         """
         set +u; source activate multiqc; set -u
         multiqc -f -o data/multiQC/all data 2> {log} 1>&2
         """
