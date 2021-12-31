
#download eukdetect database and install EukDetect
rule install_eukdetect:
    output:
        touch("taxonomy/EukDetect/eukdb_download.done")
    conda:
        "../envs/eukdetect.yaml"
    shell:
        "cd taxonomy ;"
        "git clone https://github.com/allind/EukDetect.git ;"
        "cd EukDetect ;"
        "mkdir eukdb ;"
        "cd eukdb ;"
        "wget https://ndownloader.figshare.com/files/26173346 ;"
        "tar -zxvf 26173346 ;"
        "rm 26173346 ;"
        "cd .. ;"
        "python setup.py install ;"

rule runall_eukdetect:
    input:
        "taxonomy/EukDetect/eukdb_download.done",
    output:
        touch("taxonomy/EukDetect/finished_eukdetect")
    conda:
        "../envs/eukdetect.yaml"
    log:
        "/lustre/project/rumen_longread_metagenome_assembly/tfuller/Sheep_2_HiFi_binning/logs/EukDetect2.log"
    threads: 8
    resources:
        mem_mb=250000,
        disk_mb=250000
    shell:
        "cd taxonomy/EukDetect ;"
        "mkdir output ;"
        "(eukdetect --mode runall --configfile euk_config.yaml) &> {log}"

#rule runtest_eukdetect:
#    input:
#        "taxonomy/EukDetect/eukdb_download.done",
#    output:
#        touch("taxonomy/EukDetect/finished_eukdetect")
#    conda:
#        "../envs/eukdetect.yaml"
#    log:
#        "logs/EukDetect_test.log"
#    threads: 8
#    resources:
#        mem_mb=250000,
#        disk_mb=250000
#    shell:
#        """
#        cd taxonomy/EukDetect ;
#        (python tests/test_eukdetect.py) 2> {log}
#        """
