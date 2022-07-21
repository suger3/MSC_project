###################################################################
#pipeline for mapping to pan-genome  and calculate depth
#Yizhou zhang
#09/07/2022
###################################################################
READ_DIR="/datastore/yzhang/filter/"
MAPP2PANGNEOM="/datastore/yzhang/"+"mapp2pan_genome/"
#PAN_GENOME_FASTA=
#PAN_GENOME_INDEX=BASE_DIR+"E28_assembly_index/"+
GENOME_FASTA="pan_genome51.fa"
GENOME_INDEX="/localdisk/home/s2224743/data/E28_assembly_index/"+GENOME_FASTA
DEPTH="/datastore/yzhang/"+"depth/"
DEPTH_PY="/localdisk/home/s2224743/data/snakemake_pipeline/scripts/get_depth_single.py"
E_nums, = glob_wildcards(READ_DIR+"E{E_num}_f.fq.gz")
rule all:
    message:
        print(E_nums)
    input:
        expand(DEPTH+"E{E_num}.tsv",E_num=E_nums)
###################################################################
# Do bwa_mem2 to map the short read to genome
# input raw short read |  output mapping sam file
# ###################################################################
rule BWA_MEM2:
    input:
        READ_DIR+"E{E_num,[0-9]+}_f.fq.gz",
        READ_DIR+"E{E_num,[0-9]+}_r.fq.gz"
    output:
        DEPTH+"E{E_num,[0-9]+}.tsv"
    threads:40
    shell:
        "bwa-mem2 mem -t 40 {GENOME_INDEX} {input[0]} {input[1]} | samtools view  -@ 40 -h -F 1796 | python {DEPTH_PY} - > {output}"

