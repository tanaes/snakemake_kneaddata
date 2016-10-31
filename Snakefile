import tempfile
import os
import glob


# Config parameters from config.yaml file
TMP_DIR_ROOT = config["TMP_DIR_ROOT"]
RUN = config["RUN"]
SAMPLES_PE = config["samples_pe"] if "samples_pe" in config else []
SAMPLES_SE = config["samples_se"] if "samples_se" in config else []

# Path to programs
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]

# These rules are not processor intensive, and will execute on the head node
# without being allocated to compute nodes
localrules: raw_make_links_pe, raw_make_links_se, multiQC_run, multiQC_all

# This specifices the environment setup command for running compute jobs
# (for example, source activate kneaddata) and is specified in config.yaml
ENV_KNEAD = config["KNEAD_ENV"]
shell.prefix(ENV_KNEAD + '; ')


#### Top-level rules: rules to execute a subset of the pipeline

rule all:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_kneaddata_pe
        - qc_kneaddata_se
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
            "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_{end}.fq.gz",
            sample = SAMPLES_PE,
            run = RUN,
            end = "R1 R2".split()
        ) + expand( # fastqc zip and html for raw SE data
            "data/{sample}/{run}/kneaddata/{sample}_kneaddata_{end}.fq.gz",
            sample = SAMPLES_SE,
            run = RUN,
            end = "SE".split()
        ),
        expand(
            "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_paired_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand(
            "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_single_{end}_fastqc.{extension}",
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
    """
    Makes symbolic links to original read files.
    We do this so we can be sure of file naming in downstream steps. 

    Note that readlink -f does not have the same behavior on Mac OSX, so this
    step will fail. If running on MacOSX, can ensure that config.yaml specifies
    path to original file and replace with:

    ln -s {input.forward} {output.forward} 2> {log}
    ln -s {input.reverse} {output.reverse} 2> {log}
    """
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
    """
    Makes symbolic links to original read files.
    We do this so we can be sure of file naming in downstream steps. 

    Note that readlink -f does not have the same behavior on Mac OSX, so this
    step will fail. If running on MacOSX, can ensure that config.yaml specifies
    path to original file and replace with:

    ln -s {input.single} {output.single} 2> {log}
    """
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
    Run kneaddata on paired end mode to filter host reads, eliminate Illumina
    adaptors and remove low quality regions and reads.
    """
    input:
        forward = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        reverse = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    output:
        paired_f  = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_R1.fq.gz",
        paired_r  = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_R2.fq.gz",
        unpaired_f = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_unmatched_R1.fq.gz",
        unpaired_r = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_unmatched_R2.fq.gz",
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
                    --reference-db {params.db} \
                    --quality-scores {params.phred} \
                    --threads {threads} \
                    --trimmomatic {trimmomatic} \
                    --trimmomatic-options 'ILLUMINACLIP:{params.adaptor}:2:30:10 {params.trimmomatic_params}' \
                  2> {log}

                  {gzip} %s/*

                  touch %s/{wildcards.sample}_kneaddata_paired_1.fastq.gz
                  touch %s/{wildcards.sample}_kneaddata_paired_2.fastq.gz
                  touch %s/{wildcards.sample}_kneaddata_unmatched_1.fastq.gz
                  touch %s/{wildcards.sample}_kneaddata_unmatched_2.fastq.gz
                  
                  scp %s/{wildcards.sample}_kneaddata_paired_1.fastq.gz {output.paired_f}
                  scp %s/{wildcards.sample}_kneaddata_paired_2.fastq.gz {output.paired_r}
                  scp %s/{wildcards.sample}_kneaddata_unmatched_1.fastq.gz {output.unpaired_f}
                  scp %s/{wildcards.sample}_kneaddata_unmatched_2.fastq.gz {output.unpaired_r}

                  """ % (temp_dir, temp_dir, temp_dir, temp_dir, temp_dir,
                         temp_dir, temp_dir, temp_dir, temp_dir, temp_dir))



rule qc_kneaddata_se:
    """
    Run kneaddata on single end mode to filter host reads, eliminate Illumina
    adaptors and remove low quality regions and reads.
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
                    --reference-db {params.db} \
                    --quality-scores {params.phred} \
                    --threads {threads} \
                    --trimmomatic {trimmomatic} \
                    --trimmomatic-options 'ILLUMINACLIP:{params.adaptor}:2:30:10 {params.trimmomatic_params}' \
                  2> {log}

                  {gzip} %s/*
                  
                  scp %s/{wildcards.sample}_kneaddata.fastq.gz {output.single}
                  """ % (temp_dir, temp_dir, temp_dir))


rule qc_fastqc_pe:
    """
    Do FASTQC reports
    One thread per fastq.gz file
    """
    input:
        fastq = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_paired_{end}.fq.gz",
   output:
        html = "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_paired_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_paired_{end}_fastqc.zip"
    threads:
        1
    log:
        "logs/{run}/qc/fastqc_kneaddata_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/qc/fastqc_kneaddata_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_kneaddata \
            {input.fastq} 
        2> {log} 1>&2
        """


rule qc_fastqc_se:
    """
    Do FASTQC reports
    One thread per fastq.gz file
    """
    input:
        fastq = "data/{sample}/{run}/kneaddata/{sample}_kneaddata_{end}.fq.gz"
    output:
        html = "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_single_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_kneaddata/{sample}_kneaddata_single_{end}_fastqc.zip"
    threads:
        1
    log:
        "logs/{run}/qc/fastqc_kneaddata_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/qc/fastqc_kneaddata_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_kneaddata \
            {input.fastq} 
        2> {log} 1>&2
        """

rule multiQC_run:
    """
    Run MultiQC summarization on all of the FastQC datafiles produced by the
    pipeline for a given run. 
    """
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
    """
    When multiple runs have been done, this will summarize all of the outputs
    for all of the runs. 

    Note that this will just search the entire ./data folder for FastQC output;
    the inputs are all provided to ensure that it doesn't execute until all the
    files for the current run have executed.
    """
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
