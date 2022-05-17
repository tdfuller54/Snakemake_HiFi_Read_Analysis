## Create snakemake workflow to bin HiFi contigs using illumina reads for
## alignment generation

configfile: "config.yaml"

include: "rules/taxonomy.smk"
#include: "rules/eukdetect.smk"

rule all:
    input:
        expand("sorted_reads/{assembly}/{sample}.bam.bai",
            assembly=config["assembly"], sample=config["samples"]),
        expand("binning/{assembly}/metabat2/bins",
            assembly=config["assembly"]),
        expand("binning/{assembly}/concoct/bins",
            assembly=config["assembly"]),
        expand("binning/{assembly}/DASTool/output"
            "/DASTool_results_DASTool_scaffolds2bin.tsv",
            assembly=config["assembly"]),
        expand("taxonomy/gtdb/{assembly}_classify.done",
            assembly=config["assembly"]),
        expand("binning/{assembly}/DASTool/output/checkm/result.tsv",
            assembly=config["assembly"]),
        #"taxonomy/EukDetect/finished_eukdetect"

#Use bwa and samtools to generate and sort BAM files for downstream analysis
rule bwa_index:
    input:
        "assembly/{assembly}.fa"
    output:
        "assembly/{assembly}.fa.amb",
        "assembly/{assembly}.fa.ann",
        "assembly/{assembly}.fa.bwt",
        "assembly/{assembly}.fa.pac",
        "assembly/{assembly}.fa.sa"
    log:
        "logs/bwa_index/{assembly}.log"
    conda:
        "envs/bwa_sam.yaml"
    shell:
        "bwa index {input} 2> {log}"

rule bwa_sort:
    input:
        ref = "assembly/{assembly}.fa",
        amb = "assembly/{assembly}.fa.amb",
        ann = "assembly/{assembly}.fa.ann",
        bwt = "assembly/{assembly}.fa.bwt",
        pac = "assembly/{assembly}.fa.pac",
        sa = "assembly/{assembly}.fa.sa",
        r1 = lambda wildcards: config["samples"][wildcards.sample][0],
        r2 = lambda wildcards: config["samples"][wildcards.sample][1]
    output:
        "sorted_reads/{assembly}/{sample}.bam"
    conda:
        "envs/bwa_sam.yaml"
    log:
        "logs/bwa_sort/{assembly}_{sample}.log"
    threads: config["threads"]
    shell:
        "(bwa mem  -t {threads} {input.ref} "
        "{input.r1} {input.r2} | "
        "samtools sort -o {output} -) 2> {log}"

#Create index files of sorted BAM files
rule samtools_index:
    input:
        "sorted_reads/{assembly}/{sample}.bam"
    output:
        "sorted_reads/{assembly}/{sample}.bam.bai"
    conda:
        "envs/bwa_sam.yaml"
    shell:
        "samtools index {input}"

#get depth file for Metabat2
rule get_depth_file:
    input:
        expand("sorted_reads/{assembly}/{sample}.bam",
            assembly=config["assembly"], sample=config["samples"])
    output:
        "sorted_reads/{assembly}/metabat2_depth.txt"
    log:
        "logs/metabat2/{assembly}_depth.log"
    threads:  config["threads"]
    conda:
        "envs/metabat2.yaml"
    shell:
        "(jgi_summarize_bam_contig_depths "
        "--outputDepth {output} {input}) 2> {log}"

#Run MEtabat2 with assembly file and depth file
rule metabat2:
    input:
        assembly="assembly/{assembly}.fa",
        depth="sorted_reads/{assembly}/metabat2_depth.txt"
    output:
        directory("binning/{assembly}/metabat2/bins")
    log:
        "logs/metabat2/{assembly}_metabat2.log"
    threads: config["threads"]
    conda:
        "envs/metabat2.yaml"
    shell:
        """
        (mkdir {output}
        metabat2 -i {input.assembly} -a {input.depth} \
        -t {threads} -o {output}/metabat2) \
        2> {log}
        """

#Slice contigs and get concoct coverage
rule concoct_cut:
    input:
        "assembly/{assembly}.fa"
    output:
        contigs="binning/{assembly}/concoct/contigs_10k.fa",
        bed="binning/{assembly}/concoct/contigs_10k.bed"
    log:
        "logs/concoct/{assembly}_cut.log"
    threads: 10
    conda:
        "envs/concoct.yaml"
    shell:
        "(cut_up_fasta.py {input} -c 10000 "
        "-o 0 --merge_last -b {output.bed} > {output.contigs}) 2> {log}"

rule concoct_coverage:
    input:
        mapped=expand("sorted_reads/{assembly}/{sample}.bam",
            assembly=config["assembly"], sample=config["samples"]),
        contigs="binning/{assembly}/concoct/contigs_10k.bed"
    output:
        "binning/{assembly}/concoct/coverage_table.tsv"
    log:
        "logs/concoct/{assembly}_coverage.log"
    threads: 10
    conda:
        "envs/concoct.yaml"
    shell:
        "(concoct_coverage_table.py {input.contigs} "
        "{input.mapped} > {output}) 2> {log}"

#Run concoct, merge clusters, and extract bins
rule run_concoct:
    input:
        contigs="binning/{assembly}/concoct/contigs_10k.fa",
        coverage="binning/{assembly}/concoct/coverage_table.tsv"
    output:
        dir=directory("binning/{assembly}/concoct/output/"),
        cluster="binning/{assembly}/concoct/output/clustering_gt1000.csv"
    log:
        "logs/concoct/{assembly}_concoct.log"
    threads: 10
    conda:
        "envs/concoct.yaml"
    shell:
        "(concoct --composition_file {input.contigs} "
        "--coverage_file {input.coverage} -t {threads} "
        "-b {output.dir}) 2> {log}"

rule concoct_merge:
    input:
        "binning/{assembly}/concoct/output/clustering_gt1000.csv"
    output:
        "binning/{assembly}/concoct/output/clustering_merged.csv"
    log:
        "logs/concoct/{assembly}_merge.log"
    threads: 10
    conda:
        "envs/concoct.yaml"
    shell:
        "(merge_cutup_clustering.py {input} > {output}) 2> {log}"

rule concoct_extract_bins:
    input:
        assembly="assembly/{assembly}.fa",
        cluster="binning/{assembly}/concoct/output/clustering_merged.csv"

    output:
        directory("binning/{assembly}/concoct/bins")
    log:
        "logs/concoct/{assembly}_extract_bins.log"
    threads: 10
    conda:
        "envs/concoct.yaml"
    shell:
        """
        (mkdir {output}
        extract_fasta_bins.py {input.assembly} \
        {input.cluster} --output_path {output}) \
        2> {log}
        """

#Reformat concoct and metabat2 output for DASTool input
rule csv_to_tsv:
    input:
        rules.concoct_merge.output
    output:
        "binning/{assembly}/DASTool/concoct.scaffolds2bin.tsv"
    conda:
        "envs/DASTool.yaml"
    shell:
        """
        perl -pe "s/,/\tconcoct./g;" {input} > {output}
        """

rule bins_to_scaffold:
    input:
        rules.metabat2.output
    output:
        "binning/{assembly}/DASTool/metabat2.scaffolds2bin.tsv"
    conda:
        "envs/DASTool.yaml"
    shell:
        "Fasta_to_Scaffolds2Bin.sh -i {input} "
        "-e fa > {output}"

rule DASTool:
    input:
        bins=expand("binning/{assembly}/DASTool/{binner}.scaffolds2bin.tsv",
            assembly=config["assembly"], binner=config["binners"]),
        ref="assembly/{assembly}.fa"
    output:
        scaffolds2bin="binning/{assembly}/DASTool/output"
        "/DASTool_results_DASTool_scaffolds2bin.tsv",
        bins=directory("binning/{assembly}/DASTool/output"
        "/DASTool_results_DASTool_bins")
    params:
        binners=",".join(config["binners"]),
        scaffolds2bin=lambda wildcards, input: ",".join(input.bins),
        output_prefix="binning/{assembly}/DASTool/output/DASTool_results"
    log:
        "logs/DASTool_{assembly}.log"
    threads: 10
    conda:
        "envs/DASTool.yaml"
    shell:
        """
        DAS_Tool --outputbasename {params.output_prefix} \
        --bins {params.scaffolds2bin} --labels {params.binners} \
        --contigs {input.ref} --search_engine diamond --write_bin_evals 1 \
        --write_bins 1 --create_plots 1 --threads {threads} --debug &> {log};
        mv {params.output_prefix}_DASTool_scaffolds2bin.txt {output.scaffolds2bin} &>> {log}
        """
rule download_checkm:
    output:
        touch("checkm_data/download.done")
    log:
        "logs/checkm_download.log"
    conda:
        "envs/checkm.yaml"
    shell:
        """
        (mkdir checkm_data
        pip3 install pysam
        pip3 install checkm-genome
        cd checkm_data
        wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
        tar -xzvf checkm_data_2015_01_16.tar.gz
        rm -r checkm_data_2015_01_16.tar.gz) 2> {log}
        """

rule checkm_lineage:
    input:
        rules.DASTool.output.scaffolds2bin,
        "checkm_data/download.done"
    output:
        dir=directory("binning/{assembly}/DASTool/output/checkm"),
        "binning/{assembly}/DASTool/output/checkm/lineage.ms",
    log:
        "logs/checkm_lineage_{assembly}.log"
    threads: 16
    conda:
        "envs/checkm.yaml"
    shell:
        """
        (mkdir {output.dir}
        checkm data setRoot {workflow.basedir}/checkm_data
        cd binning/{assembly}/DASTool/output
        checkm lineage_wf -t {threads} -x fa DASTool_results_DASTool_bins checkm) \
        2> {log}
        """

rule checkm_qa:
    input:
        "binning/{assembly}/DASTool/output/checkm/lineage.ms"
    output:
        "binning/{assembly}/DASTool/output/checkm/result.tsv"
    log:
        "logs/checkm_qa_{assembly}.log"
    threads: 16
    conda:
        "envs/checkm.yaml"
    shell:
        """
        (cd binning/{assembly}/DASTool/output
        checkm qa -t 16 --out_format 1 --tab_table \
        -f checkm/result.tsv checkm/lineage.ms checkm) \
        2> {log}
        """
