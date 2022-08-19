"""
These rules will trim sequencing reads
by removing Illumina adapters and tiling primers

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned bam file
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
        # if the reads are already unzipped, gunzip will just pass the file through
        """
        cat {input} | gunzip -c > {output}
        """ 
        

# remove Illumina adapters and trim for quality using trimmomatic
# using single end mode as pairs are now be combined
rule trimmomatic:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        "01_raw/{sample}_combined.fastq"
    output:
        "02_trim/{sample}_trim.fastq"
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
        trimmomatic SE \
            -threads {threads} \
            {input} {output} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """


# to trim the tiling primers by position, we need a sorted bam file
# first create a bowtie index of the JEV reference sequence
rule build_bowtiedb:
    message:
        """
        ** preprocessing **
        Building a bowtie2 database of the reference genome
        """
    input:
        config["program_dir"] + config["JEV_reference_genome"]
    output:
        # bowtie2-build needs a basename for the database
        # often just the same name as the input
        # and it appends several *bt2 files
        # will trick snakemake by using this as an output even though
        # I won't use it in the shell command
        config["program_dir"] + config["JEV_reference_genome"] + ".1.bt2"
    threads: 8
    shell:
        # use the same name for basename reference database
        """
        bowtie2-build \
            --threads {threads} \
            {input} \
            {input} > /dev/null
        """

rule bowtie_to_reference:
    message:
        """
        ** preprocessing **
        Mapping {wildcards.sample} reads to reference genome
        """
    input:
        reads = "02_trim/{sample}_trim.fastq",
        db_trick = config["program_dir"] + config["JEV_reference_genome"] + ".1.bt2"
    output:
        sam_fl = "03_mapping/{sample}_reference.sam"
    params:
        assembly_db = config["program_dir"] + config["JEV_reference_genome"] 
    log:
        "logs/bowtie/{sample}.log"
    benchmark:
        "benchmarks/bowtie/{sample}.txt"
    threads: 16
    shell:
        """
        bowtie2 \
            -x {params.assembly_db} \
            -U {input.reads} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule sam_to_sorted_bam:
    message:
        """
        ** preprocessing **
        Converting {wildcards.sample} sam file to bam
        """
    input:
        "03_mapping/{sample}_reference.sam"
    output:
        "03_mapping/{sample}_reference.sorted.bam"
    threads: 8
    shell:
        """
        samtools view \
            -@ {threads} \
            -S -b \
            {input} | \
        samtools sort \
            - \
            -o {output}
        """

rule trim_primers_ivar:
    message:
        """
        ** preprocessing **
        Triming tiling primers from {wildcards.sample} with iVar
        """
    input:
        "03_mapping/{sample}_reference.sorted.bam"
    output:
        "03_mapping/{sample}_reference.primer_trim.bam"
    params:
        primers_bed = config["program_dir"] + config["primal_scheme_bed"],
        ivar_min_length = config["ivar_min_length"]
    log:
        "logs/ivar_primers/{sample}.log"
    threads: 8
    shell:
        # the -e flag will keep reads even if they don't contain primers
        # this is possible if people use tagmentation - ie. Nextera
        """
        ivar trim \
            -i {input} \
            -b {params.primers_bed} \
            -m {params.ivar_min_length} \
            -e \
            -p {output} > {log}
        """

# need to sort the bam file again after trimming of primers
rule trimmed_bam_to_sorted:
    message:
        """
        ** preprocessing **
        Re-sorting trimmed bam for {wildcards.sample}
        """
    input:
        "03_mapping/{sample}_reference.primer_trim.bam"
    output:
        "03_mapping/{sample}_reference.primer_trim.sorted.bam"
    threads: 8
    shell:
        """
        samtools sort \
            {input} \
            -o {output}
        """


