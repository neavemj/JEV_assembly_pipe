"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned Illumina sequencing reads
"""


# get the sample and file names from the config file
def getFastq(wildcards):
    return config['samples'][wildcards.sample]


# combine paired reads together and unzip for downstream analysis
rule combine_pairs:
    input:
        reads = getFastq
    output:
        "01_raw/{sample}_combined.fastq"
    shell:
        # if we only have single reads they will just be copied and unzipped
        """
        cat {input} | gunzip -c > {output}
        """ 
        

# remove Illumina adapters and trim for quality using trimmomatic
rule trimmomatic_PE:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        reads = getFastq,
    output:
        R1_P = "01_trim/{sample}_1P.fastq.gz",
        R1_U = "01_trim/{sample}_1U.fastq.gz",
        R2_P = "01_trim/{sample}_2P.fastq.gz",
        R2_U = "01_trim/{sample}_2U.fastq.gz"
    params:
        qual = config["trimmomatic_quality"],
        adapters = config["program_dir"] + config["trimmomatic_adapters"],
        minlen = config["trimmomatic_minlen"]
    threads: 8
    log:
        "logs/trimmomatic_PE/{sample}.log"
    benchmark:
        "benchmarks/trimmomatic_PE/{sample}.txt"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.reads} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """

