
rule download_gtdb:
    output:
        touch("taxonomy/gtdb/gtdb_download.done")
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        "mkdir taxonomy/gtdb ;"
        "download-db.sh ;"

#Need to find out where bins are written to from DASTool
#based on my code should be in binning/{assembly}/DASTool/output/DASTool_results_DASTool_bins/
rule identify:
    input:
        dir="binning/{assembly}/DASTool/output"
        "/DASTool_results_DASTool_bins",
        db=rules.download_gtdb.output
    output:
        directory("taxonomy/gtdb/{assembly}/identify")
    log:
        "logs/gtdb_{assembly}/identify.log"
    threads: config["threads"]
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        mkdir {output};
        gtdbtk identify --genome_dir {input.dir} \
        --out_dir {output} --extension fa \
        --cpus {threads} &> {log}
        """

rule align:
    input:
        rules.identify.output
    output:
        directory("taxonomy/gtdb/{assembly}/align")
    log:
        "logs/gtdb_{assembly}/align.log"
    threads: config["threads"]
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        mkdir {output};
        gtdbtk align --identify_dir {input} \
        --out_dir {output} --cpus {threads} &> {log}
        """

rule classify:
    input:
        dir="binning/{assembly}/DASTool/output"
        "/DASTool_results_DASTool_bins",
        align=rules.align.output
    output:
        directory("taxonomy/gtdb/{assembly}/classify"),
        touch("taxonomy/gtdb/{assembly}_classify.done")
    log:
        "logs/gtdb_{assembly}/classify.log"
    threads: config["threads"]
    resources:
        mem_mb=250000,
        disk_mb=250000
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        mkdir {output[0]};
        gtdbtk classify --genome_dir {input.dir} \
        --align_dir {input.align} --out_dir {output[0]} \
        --extension fa --cpus {threads} &> {log}
        """
