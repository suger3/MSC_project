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
Pan_genome=BASE_DIR+"pan_genome/"
iterate_reference="/localdisk/home/s2224743/data/E28_assembly_index/E28_assembly_over200.fa"
cd_hit_iterate_reference_sh="/localdisk/home/s2224743/data/snakemake_pipeline/scripts/cd_hit_iterate_reference.sh"
minimap2_iterate_reference_sh="/localdisk/home/s2224743/data/snakemake_pipeline/scripts/minimap2_iterate_reference.sh"
REPORT=Pan_genome+"all_report.csv"
write_report_py="/localdisk/home/s2224743/data/snakemake_pipeline/scripts/write_report.py"
###################################################################
# The list of samples to be processed
###################################################################
#E_nums, = glob_wildcards(READ_DIR+"E{E_num}_f.fq.gz")

#use for test
#E_nums=["025"]
paired=["f","r"]

import os
if os.path.exists("tmp"):
    pass
else:
    os.makedirs("tmp")

# def process_lines(file_name):
#     """generates id/run, ignoring non-numeric lines"""
#     with open("data_list", "r") as f:
#         for line in f:
#             E_nums, double_copy, *_ = line.split()
#             if E_nums.isnumeric() and double_copy.isnumeric():
#                 E_nums = E_nums.zfill(8)
#                 single_copy_min=int(double_copy/4)
#                 single_copy_max=int(double_copy*5)
#                 single_copy_min = single_copy_min.zfill(8)
#                 single_copy_max = single_copy_max.zfill(8)
#                 yield E_nums, single_copy_max , single_copy_min
import pandas as pd

def add_new_attribute(pd_name,list,col_name):
    copy_time=len(list)
    data_num=pd_name.shape[0]
    tmp_pd=pd_name.copy()
    for i in range(copy_time-1):
        tmp_pd=pd.concat([tmp_pd,pd_name],ignore_index=True)
    col_data=[]
    for i in list:
        col_data+=[i for j in range(data_num)]
    col_data=pd.Series(col_data)
    tmp_pd[col_name]=col_data
    #all=pd.concat([tmp_pd,col_data],axis=1)
    return tmp_pd

run_list = pd.read_csv("data_list",converters={"Enum":str,"double_copy":float})
run_list["single_copy_min"]=(run_list["double_copy"]/4).astype("int32")
run_list["single_copy_max"]=(run_list["double_copy"]*5).astype("int32")
#run_list=add_new_attribute(run_list,paired,"pair")
print(run_list)


rule all:
    #message:
        #print(E_nums)
    input:
        #expand(ASSEMBLY_RESULT+"E{E_num}_spades_ci{min}_cx{max}",zip,Pair=run_list.pair,\
        #E_num=run_list.Enum,min=run_list.single_copy_min,max=run_list.single_copy_max),
        #Pan_genome+"cd_hit_finall.fa",
        Pan_genome+"minimap2_finall.fa",
        Pan_genome+"all_report.csv"
###################################################################
# fastp filter and QC
# input raw short read |  output fastp file
###################################################################
# rule FASTP_FILTER:
#     input:
#         READ_DIR+"E{E_num}_f.fq.gz",
#         READ_DIR+"E{E_num}_r.fq.gz"
#     output:
#         temp(BASE_DIR+"filter/E{E_num}_f.fq.gz"),
#         temp(BASE_DIR+"filter/E{E_num}_r.fq.gz"),
#         BASE_DIR+"fastp_report/E{E_num}.json",
#         BASE_DIR+"fastp_report/E{E_num}.html"
#     log:
#         LOG_DIR+"E{E_num}_fastp.log"
#     threads:8
#     shell:
#         "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -z 5 -q 20 -u 30 --thread {threads}\
#              -j {output[2]} -h {output[3]} >{log} 2>&1 "

###################################################################
# Do single copy filter
# input raw short read |  output single copy reads
# ###################################################################
rule KMC_filter:
    input:
        READ_DIR+"E{E_num,[0-9]+}_{Pair,[fr]}.fq.gz",
        KMER_KMC_DIR+"E{E_num}_21er.kmc_suf",
        KMER_KMC_DIR+"E{E_num}_21er.kmc_pre"
    output:
        temp(READ_DIR+"E{E_num,[0-9]+}_{Pair,[fr]}_ci{min}_cx{max}.fq")
    log:
        LOG_DIR+"E{E_num}_{Pair,[fr]}_kmc_filter_ci{min}_cx{max}.log"
    threads:15
    params:
        dbname = KMER_KMC_DIR+"E{E_num}_21er"
    shell:
        "kmc_tools -t{threads} filter {params.dbname} -ci{wildcards.min} -cx{wildcards.max}  {input[0]} -ci0.9 {output} >{log} 2>&1 "

###################################################################
# Prepare for KMC
###################################################################
rule PRE_KMC:
    input:
        READ_DIR+"E{E_num,[0-9]+}_f.fq",
        READ_DIR+"E{E_num,[0-9]+}_r.fq"
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
    threads:30
    shell:
        "kmc -k21 -t30  @{input} {params.dbname}  tmp/ > {log}  2>&1"

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
        "kmc_tools -t{threads} transform {params.dbname} histogram {output} > {log}  2>&1"


###################################################################
# Do bwa_mem2 to map the short read to genome
# input raw short read |  output mapping sam file
# ###################################################################
rule BWA_MEM2:
    input:
        READ_DIR+"E{E_num,[0-9]+}_{Pair,[fr]}_ci{min}_cx{max}.fq"
    output:
        temporary(BASE_DIR+"aligned/E{E_num,[0-9]+}_{Pair,[fr]}_ci{min}_cx{max}.bam")
    log:
        LOG_DIR+"E{E_num}_{Pair,[fr]}_BWA_ci{min}_cx{max}.log"
    threads:15
    shell:
        "bwa-mem2 mem -t {threads} {GENOME_INDEX} {input[0]} | samtools view -bS -o {output} >{log} 2>&1 "

# ###################################################################
# # convert sam to bam
# # input sam file  | output sorted bam file
# ###################################################################
rule SORT_BAM:
    input:
        BASE_DIR+"aligned/E{E_num,[0-9]+}_{Pair,[fr]}_ci{min}_cx{max}.bam"
    output:
        temp(BASE_DIR+"aligned/E{E_num,[0-9]+}_{Pair,[fr]}_sorted_ci{min}_cx{max}.bam")
    threads:10
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp -@ {threads} {input}"
#
rule BAM_INDEX:
    input:
        BASE_DIR+"aligned/E{E_num,[0-9]+}_{Pair,[fr]}_sorted_ci{min}_cx{max}.bam"
    output:
        temp(BASE_DIR+"aligned/E{E_num,[0-9]+}_{Pair,[fr]}_sorted_ci{min}_cx{max}.bam.bai")
    threads:10
    shell:
        "samtools index -m {threads} {input} {output}"

###################################################################
# extract the indetifier of single copy
###################################################################
rule POOR_MAP:
    input:
        BASE_DIR+"aligned/E{E_num,[0-9]+}_f_sorted_ci{min}_cx{max}.bam",
        BASE_DIR+"aligned/E{E_num,[0-9]+}_r_sorted_ci{min}_cx{max}.bam",
        BASE_DIR+"aligned/E{E_num,[0-9]+}_f_sorted_ci{min}_cx{max}.bam.bai",
        BASE_DIR+"aligned/E{E_num,[0-9]+}_r_sorted_ci{min}_cx{max}.bam.bai"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_ci{min}_cx{max}.namelist"),
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_r_poormap_ci{min}_cx{max}.namelist")
    #params:
    #    dbname = KMER_KMC_DIR+"E{E_num}_21er.kmer_num"
    threads:2
    shell:
        "python3 {GET_POOR_MAP_PY} {input[0]} {input[1]} {output[0]} {output[1]}"

# ###################################################################
# # remove duplicate identifier
# ###################################################################
rule RM_dup_Name:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.namelist"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.del_dup_namelist")
    shell:
        "sort {input} | uniq > {output}"

# ###################################################################
# # extract the identifier fastq for blast against Univec
# ###################################################################
rule EXTRACT_DEL_DUP_FASTQ:
    input:
        READ_DIR+"E{E_num,[0-9]+}_{Pair,[fr]}.fq.gz",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.del_dup_namelist"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.fq.gz")
    shell:
        "{SEQTK_PATH} subseq {input[0]} {input[1]} | pigz -5c > {output}"

# ###################################################################
# # convert fastq to fasta
# ###################################################################
rule FASTQ2FASTA:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.fq.gz"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.fa")
    shell:
        "{SEQTK_PATH} seq -a {input} > {output}"

# ###################################################################
# # blast against Univec
# ###################################################################
rule UNIVEC:
    input: 
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.fa"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.result")
    threads:4
    shell:
        "blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true \
        -evalue 700 -searchsp 1750000000000 -db {UNIVEC_DATABASE_PATH} -outfmt 6 -num_threads {threads} -query {input} > {output}"

# ###################################################################
# # Get the identifiers of  contamination reads
# ###################################################################
rule DEL_UNIVEC_LIST:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.result"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_univec_ci{min}_cx{max}.namelist")
    shell:
        "cut -f1 {input} | sort | uniq > {output}"


# ###################################################################
# # delete the identifiers of  contamination reads from single copy list
# ###################################################################
rule DEL_UNIVEC_FQ:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_ci{min}_cx{max}.del_dup_namelist",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_univec_ci{min}_cx{max}.namelist"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_del_univec_dup_namelist_ci{min}_cx{max}")
    shell:
        "sort {input[0]} {input[1]} | uniq -u > {output} "


# ###################################################################
# # conbine the f and r identifier list
# ###################################################################
rule COMBINE_FR:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_del_univec_dup_namelist_ci{min}_cx{max}",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_r_poormap_del_univec_dup_namelist_ci{min}_cx{max}"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_combine_poormap_del_univec_dup_namelist_ci{min}_cx{max}")
    shell:
        "sort {input[0]} {input[1]} | uniq > {output}"

# ###################################################################
# # extract single copy paired reads
# ###################################################################
rule EXTRACT_DEL_DUP_UNIVEC_FASTQ:
    input:
        READ_DIR+"E{E_num,[0-9]+}_{Pair,[fr]}.fq.gz",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_combine_poormap_del_univec_dup_namelist_ci{min}_cx{max}"
    output:
        temp(BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_{Pair,[fr]}_poormap_del_dup_univec_ci{min}_cx{max}.fq")
    shell:
        "{SEQTK_PATH} subseq {input[0]} {input[1]} > {output}"


# ###################################################################
# # do platanus assembly
# ###################################################################
rule PLATANUS_ASSEMBLY:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_del_dup_univec_ci{min}_cx{max}.fq",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_r_poormap_del_dup_univec_ci{min}_cx{max}.fq"
    output:
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_contig_ci{min}_cx{max}.fa",
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_contigBubble_ci{min}_cx{max}.fa"
    threads:8
    params:
        #prefix of output files
        prefix=ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_ci{min}_cx{max}"
    log:
        LOG_DIR+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_assemble_ci{min}_cx{max}.log"
    shell:
        "platanus assemble -o {params.prefix} -f {input[0]} {input[1]} -t {threads} -k {wildcards.PLATANUS_kmer} -c {wildcards.PLATANUS_c} -tmp tmp_platanus/  >{log} 2>&1"

rule PLATANUS_scaffold:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_del_dup_univec_ci{min}_cx{max}.fq",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_r_poormap_del_dup_univec_ci{min}_cx{max}.fq",
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_contig_ci{min}_cx{max}.fa",
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_contigBubble_ci{min}_cx{max}.fa"
    output:
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_scaffold_ci{min}_cx{max}.fa",
        ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_scaffoldBubble_ci{min}_cx{max}.fa"
    threads:8
    params:
        #prefix of output files
        prefix=ASSEMBLY_RESULT+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_ci{min}_cx{max}"
    log:
        LOG_DIR+"E{E_num,[0-9]+}_platanus_k{PLATANUS_kmer,[0-9]+}_c{PLATANUS_c,[0-9]+}_scaffold_ci{min}_cx{max}.log"
    threads:10
    shell:
        "platanus scaffold -o {params.prefix} -c {input[2]} -b {input[3]} -IP1 {input[0]} {input[1]}  -t 10 -tmp tmp_platanus/  >{log} 2>&1"

# ###################################################################
# # do  spades assembly
# ###################################################################
rule SPADES_ASSEMBLY:
    input:
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_del_dup_univec_ci{min}_cx{max}.fq",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_r_poormap_del_dup_univec_ci{min}_cx{max}.fq"
    output:
        directory(ASSEMBLY_RESULT+"E{E_num,[0-9]+}_spades_ci{min}_cx{max}")
    threads:45
    log:
        LOG_DIR+"E{E_num,[0-9]+}_spades_ci{min}_cx{max}.log"
    shell:
        "spades.py --only-assembler -1 {input[0]} -2 {input[1]} -t {threads} -m 350  -o {output[0]} >{log} 2>&1"

rule MV_scaffold:
    input:
        directory(ASSEMBLY_RESULT+"E{E_num,[0-9]+}_spades_ci{min,[0-9]+}_cx{max,[0-9]+}")
    output:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}.fa"
    params:
        scaff=ASSEMBLY_RESULT+"E{E_num,[0-9]+}_spades_ci{min,[0-9]+}_cx{max,[0-9]+}/scaffolds.fasta"
    shell:
        "sed 's/^>NODE/>E{wildcards.E_num}_node/g' {params.scaff} > {output}"

rule CD_hit_remove_self:
    input:
        Pan_genome+"E{E_num}_ci{min}_cx{max}.fa"
    output:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_similar98.fa"
    log:
        LOG_DIR+"E{E_num,[0-9]+}_ci{min}_cx{max}_cd-hit_self.log"
    threads:15
    shell:
        "cd-hit-est -i {input} -o {output} -c 0.98 -n 11 -d 0 -M 0 -T {threads} >{log} 2>&1"

rule CD_hit_self_over200:
    input:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_similar98.fa"
    output:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_similar98_over200.fa"
    shell:
        "seqkit seq -m 200 {input} -o {output} "

rule get_cd_hit_iterate_reference:
    output:
        temp("get_cd_hit_iterate_reference.done")
    shell:
        "cp {iterate_reference} {Pan_genome}cd_hit_reference.fa > get_cd_hit_iterate_reference.done && "
        "cp {cd_hit_iterate_reference_sh} {Pan_genome} "

rule get_minimap2_iterate_reference:
    output:
        temp("get_minimap2_iterate_reference.done")
    shell:
        "cp {iterate_reference} {Pan_genome}minimap2_reference.fa > get_minimap2_iterate_reference.done && "
        "cp {minimap2_iterate_reference_sh} {Pan_genome} "

rule cd_hit_iterate_reference:
    input:
        expand(Pan_genome+"E{E_num}_ci{min}_cx{max}_similar98_over200.fa",zip,\
        E_num=run_list.Enum,min=run_list.single_copy_min,max=run_list.single_copy_max),
        "get_cd_hit_iterate_reference.done"
    threads:40
    output:
        Pan_genome+"cd_hit_finall.fa"
    shell:
        "bash {Pan_genome}cd_hit_iterate_reference.sh {input}"

rule pre_minimap2_diff_to_E28:
    input:
        temp("get_minimap2_iterate_reference.done")
    output:
        Pan_genome+"ref.mmi"
    shell:
        "minimap2 -d {output} {Pan_genome}minimap2_reference.fa"

rule minimap2_diff_to_E28:
    input:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_similar98_over200.fa"
    output:
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_similar98_over200.diff_E28.fa"
    threads:15
    shell:
        "minimap2 -t 15 -a {Pan_genome}ref.mmi | samtools fasta -f4 -@15 - >  {output}"

rule minimap2_iterate_reference:
    input:
     "get_minimap2_iterate_reference.done",
        expand(Pan_genome+"E{E_num}_ci{min}_cx{max}_similar98_over200.diff_E28.fa",zip,\
        E_num=run_list.Enum,min=run_list.single_copy_min,max=run_list.single_copy_max)
    threads:40
    params:
        diff_path=BASE_DIR+"pan_genome/"
    output:
        Pan_genome+"minimap2_finall.fa"
    shell:
        "bash {Pan_genome}minimap2_iterate_reference.sh {params.diff_path}"

rule report_start:
    output:
        temp(Pan_genome+"title.csv")
    shell:
        "python {write_report_py} start {output} > {output}"

rule report:
    input:
        Pan_genome+"title.csv",
        BASE_DIR+"fastp_report/E{E_num}.json",
        READ_DIR+"E{E_num,[0-9]+}_f_ci{min}_cx{max}.fq",
        READ_DIR+"E{E_num,[0-9]+}_r_ci{min}_cx{max}.fq",
        BASE_DIR+"poor_map_bam/E{E_num,[0-9]+}_f_poormap_del_dup_univec_ci{min}_cx{max}.fq",
        Pan_genome+"E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}.fa",
        Pan_genome+"E{E_num}_ci{min}_cx{max}_similar98_over200.diff_E28.fa"
    output:
        "E{E_num,[0-9]+}_ci{min,[0-9]+}_cx{max,[0-9]+}_report.done"
    shell:
         "python {write_report_py} {wildcards.E_num} {input[0]} {input[1]}  {wildcards.min}  {wildcards.max} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} > {output}"

rule report_done:
    input:
        Pan_genome+"title.csv",
        expand("E{E_num}_ci{min}_cx{max}_report.done",zip,\
        E_num=run_list.Enum,min=run_list.single_copy_min,max=run_list.single_copy_max)
    output:
        Pan_genome+"all_report.csv"
    shell:
        "cat {input} > {output}"
