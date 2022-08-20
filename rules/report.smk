"""
These rules will provide some stats and short report

The inputs are:
    - log files and bam mapping files
The outputs are:
    - stats table and coverage plots
"""


rule summarise_logs:
    message:
        """
        ** reporting **
        Summarising results from log files
        """
    input:
        trim = "logs/trimmomatic_PE/{sample}.log",
        bowtie = "logs/bowtie/{sample}.log",
        primers = "logs/ivar_primers/{sample}.log",
        consensus = "logs/ivar_consensus/{sample}.log",
    output:
        "05_report/{sample}_stats.txt"
    threads: 8
    run:
        # will use python to go through the fairly inconsistent
        # log files produced by the different tools

        # first grab useful info from the trimmomatic report        
        with open(input.trim) as t:
            for line in t:
                line = line.strip()
                if line.startswith("Input"):
                    cols = line.split(" ")
                    input_reads = cols[2]
                    surviving_reads = cols[4]
                    stats_list = [input_reads, surviving_reads]

        # now the bowtie mapping results
        with open(input.bowtie) as b:
            for line in b:
                line = line.strip()
                if line.endswith("time"):
                    aligned_once = line.split(" ")[0]
            # note, could also collect multi-mappers
            # but that should be rare with the JEV genome
            stats_list = stats_list + [aligned_once]
        
        # get primer trimming information
        with open(input.primers) as p:
            for line in p:
                line = line.strip()
                if line.startswith("Trimmed"):
                    primers_trimmed = line.split("(")[1].split(")")[0]
           # decided not to report this for now
           # its in the logs anyway
           # stats_list = stats_list + [primers_trimmed]
        
        # get consensus genome info
        with open(input.consensus) as c:
            for line in c:
                line = line.strip()
                if line.startswith("Reference"):
                    ref_len = line.split(" ")[2]
                elif line.startswith("Positions with 0"):
                    pos_0 = line.split(" ")[4]
                elif line.startswith("Positions with depth"):
                    pos_low = line.split(" ")[5]
            # to figure out genome coverage,
            # divide the regions below threshold coverage by the reference length
            cov = round((100 - ((int(pos_low) / int(ref_len)) * 100)), 2)
            stats_list = stats_list + [str(cov)]                  

        # write that information into a file
        with open(output[0], "w") as out:
            out.write(wildcards.sample + "\t" + "\t".join(stats_list))


rule make_table:
    input:
        expand("05_report/{sample}_stats.txt", sample=config["samples"])
    output:
        "05_report/summary_stats.txt"
    run:
        # now combine the individual samples into a single table
        # write to file and also do a 'pretty print' to screen
        from tabulate import tabulate
        
        with open(output[0], "w") as fl:
            header = ["Sample_name", "Input_reads", "Cleaned_reads", "Aligned_reads", "Ref_coverage_(%)"]
            fl.write("\t".join(header) + "\n")
            table_data = []
            for f in input:
                with open(f) as samp:
                    sample_line = samp.readlines()
                    sample_list = sample_line[0].split("\t")
                    table_data.append(sample_list)
                    fl.write(sample_line[0] + "\n")
            # also print the stats to the terminal
            print(tabulate(table_data, headers=header))
            

rule get_coverage:
    message:
        """
        ** reporting **
        Calculating genome coverage for {wildcards.sample}
        """
    input:
        "03_mapping/{sample}_reference.primer_trim.sorted.bam"
    output:
        "03_mapping/{sample}_reference.primer_trim.sorted.depth"
    threads: 8
    shell:
        # -a means output all positions even those with 0 depth
        # that will make it easier to plot over the genome
        """
        samtools depth \
            -@ {threads} \
            -a \
            {input} > {output}
        """


rule plot_coverage:
    input:
        depth = "03_mapping/{sample}_reference.primer_trim.sorted.depth",  
    params:
        cut_off = config["ivar_min_depth"],
        primer_bed = config["program_dir"] + config["primal_scheme_bed"]        
    output:
        pdf = "05_report/{sample}_coverage.pdf",
        png = "05_report/{sample}_coverage.png",
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_genome_coverage.R \
            {wildcards.sample} {input} {params.cut_off} {params.primer_bed} {output.pdf} {output.png}
        """












