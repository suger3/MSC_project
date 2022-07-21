###################################################################
#pipeline for short reads mapping to reference genome  and do kmc
#Yizhou zhang
#12/05/2022
###################################################################

BASE_DIR = "/datastore/yzhang/"
READ_DIR=BASE_DIR+"filter/"
GENOME_FASTA="E28_assembly.fa"
GENOME_INDEX="/localdisk/home/s2224743/data/E28_assembly_index/"+GENOME_FASTA
#GET_POOR_MAP_PY="/localdisk/home/s2224743/snakemake_test/filter_bwa_add_dep_kmer.py"
GET_POOR_MAP_PY="/localdisk/home/s2224743/data/snakemake_pipeline/scripts/get_unmap_read.py"
SEQTK_PATH="/localdisk/home/s2224743/software/seqtk/seqtk"
UNIVEC_DATABASE_PATH="/localdisk/home/s2224743/UniVec_database/UniVec"
ASSEMBLY_RESULT=BASE_DIR+"assembly/"
KMER_KMC_DIR = BASE_DIR + "KMC_dataset/"
LOG_DIR = BASE_DIR + "log/"
###################################################################
# The list of samples to be processed
###################################################################
E_nums, = glob_wildcards(READ_DIR+"E{E_num}_f.fq.gz")
#E_nums=["031"]
paired=["f","r"]

import os
if os.path.exists("tmp"):
    pass
else:
    os.makedirs("tmp")


rule all:
    input:
        expand(KMER_KMC_DIR+"hist/E{E_num}_21er.histogram",E_num=E_nums)
        #expand(KMER_DIR+"E{E_num}_{KMERS}er_{MINI_K}.kmer_num",E_num=E_nums,KMERS=KMERS,MINI_K=MIN_K_NUM),
        #expand(PICTURE_DIR+"E{E_num}_{KMERS}er_{MINI_K}_genomescope2",E_num=E_nums,KMERS=KMERS,MINI_K=MIN_K_NUM)

###################################################################
# Prepare for KMC
###################################################################
rule PRE_KMC:
    input:
        READ_DIR+"E{E_num,[0-9]+}_f.fq.gz",
        READ_DIR+"E{E_num,[0-9]+}_r.fq.gz"
    output:
        temporary(READ_DIR+"E{E_num,[0-9]+}.txt")
    shell:
        "ls {READ_DIR}E{wildcards.E_num}*.fq.gz > {output}"

###################################################################
# Do KMC to get kmc database
# input raw short read or fasta |  output kmc.pre kmc.sub
###################################################################
rule KMC_dataset:
    input:
        READ_DIR+"E{E_num,[0-9]+}.txt"
    output:
        KMER_KMC_DIR+"E{E_num}_21er.kmc_suf",
        KMER_KMC_DIR+"E{E_num}_21er.kmc_pre"
    params:
        dbname = KMER_KMC_DIR+"E{E_num}_21er"
    log:
        LOG_DIR+"E{E_num}_21er_kmc.log"
    threads:25
    shell:
        "kmc -k21 -t25 -m200  @{input} {params.dbname}  tmp/ > {log}  2>&1"

###################################################################
# deal wtih KMC dataset
# input kmc kmer file  | output sorted kmer txt file
###################################################################
rule KMC_transfer:
    input:
        KMER_KMC_DIR+"E{E_num}_21er.kmc_suf",
        KMER_KMC_DIR+"E{E_num}_21er.kmc_pre"      
    params:
        dbname = KMER_KMC_DIR+"E{E_num}_21er"
    output:
        KMER_KMC_DIR+"hist/E{E_num}_21er.histogram"
    log:
        LOG_DIR+"E{E_num}_21er_transfer.log"
    threads:5
    shell:
        "kmc_tools -t{threads} transform {params.dbname} histogram {output[0]} > {log}  2>&1"
