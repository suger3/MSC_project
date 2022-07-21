import json
import sys
import csv
from subprocess import Popen,PIPE
def get_base_num(file):
    with open(file) as f:
        data = json.load(f)
    return data["summary"]["before_filtering"]["total_bases"],data["summary"]["after_filtering"]["total_bases"]
def get_single_copy(file):  # sourcery skip: use-fstring-for-formatting
    p = Popen("grep -c \"+\" {}".format(file),shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    return int(stdout)

def get_assemble_Data(file):
    p = Popen("seqkit stat -T {} | tail -n 1 | sed 's/,//g' | cut -f4,5 ".format(file),shell=True,stdout=PIPE)
    stdout,stderr=p.communicate()
    return (str(stdout)[2:-3].split("\\t"))

if sys.argv[1]=="start":
    print("Enum,raw_read_Base,clean_read_Base,clean_coverage,ci,cx,single_copy_reads_F+R,Assemble_read_F+R,Assemble_contigs,Assemble_span,Assemble_unmapped_contigs,Assemble_unmapped_span\n")
else:
    #1: Enum string 2:csv file 3:fastp_json 4:ci num 5 cx num 6,7:single_fastq 8 assemble fastq 9:assemble.fa 10:to_assemble.fa
    Enum=sys.argv[1]
    raw_read_Base,clean_read_Base=get_base_num(sys.argv[3])
    clean_coverage=float(int(clean_read_Base)/800000000)
    ci=sys.argv[4]
    cx=sys.argv[5]
    single_copy_reads_F=get_single_copy(sys.argv[6])
    single_copy_reads_R=get_single_copy(sys.argv[7])
    Assemble_read_FR=get_single_copy(sys.argv[8])*2
    Assemble_contigs,Assemble_span=get_assemble_Data(sys.argv[9])
    Assemble_unmapped_contigs,Assemble_unmapped_span=get_assemble_Data(sys.argv[10])
    print("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(Enum,raw_read_Base,clean_read_Base,clean_coverage,ci,cx,single_copy_reads_F+single_copy_reads_R,Assemble_read_FR,Assemble_contigs,Assemble_span,Assemble_unmapped_contigs,Assemble_unmapped_span))
