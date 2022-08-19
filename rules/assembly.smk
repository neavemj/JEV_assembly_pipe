"""
These rules will assemble a consensus genome
from reads mapped to a reference

The inputs are:
    - cleaned, sorted bam file
The outputs are:
    - consensus genome in fasta format
"""


rule ivar_consensus:
    message:
        """
        ** assembly **
        Assembling a consensus genome from {wildcards.sample} with iVar
        """
    input:
        "03_mapping/{sample}_reference.primer_trim.sorted.bam"
    output:
        # ivar only takes a prefix for a file name
        # then produces *fa and *qual.txt files
        # putting this name here but can't use it in the shell
        "04_genomes/{sample}_consensus.fa"
    params:
        # putting the actual prefix name here to avoid the
        # snakemake 'waiting for output file' error
        output_prefix = "04_genomes/{sample}_consensus",
        frequency_threshold = config["ivar_frequency_threshold"],
        min_depth = config["ivar_min_depth"],
    log:
        "logs/ivar_consensus/{sample}.log"
    threads: 8
    shell:
        # the samtools flags are
        # -a, output all positions (even those with 0 depth)
        # -d, maximum depth to report (to save memory)
        # -A, do not discard anomalous read pairs
        # ivar flags are explained in the config.yaml file     
        """
        samtools mpileup \
            -a \
            -d 1000 \
            -A \
            {input} | \
        ivar consensus \
            -t {params.frequency_threshold} \
            -m {params.min_depth} \
            -p {params.output_prefix} > {log}
        """
