# -*- coding: utf-8 -*-
from pickle import NONE
import sys
import pysam
import time
import re
Mpattern=re.compile(r"(\d+)M")
Dpattern=re.compile(r"(\d+)D")
Ipattern=re.compile(r"(\d+)I")
Hpattern=re.compile(r"(\d+)H")
Spattern=re.compile(r"(\d+)S")
bwa_file_f=sys.argv[1]
bwa_file_r=sys.argv[2]
output_f=sys.argv[3]
output_r=sys.argv[4]

#print("bwa file: {}\nkmer_dis file: {}\nkmer_len: {}\n".format(bwa_file,kmer_dis,kmer_length))
#this file to get the kmer coverage
#model_file=sys.argv[4]
def get_cigarstring(cigarstring):
    M=Mpattern.findall(cigarstring)
    D=Dpattern.findall(cigarstring)
    I=Ipattern.findall(cigarstring)
    H=Hpattern.findall(cigarstring)
    S=Spattern.findall(cigarstring)
    Mnum=sum(int(i) for i in M )
    Dnum=sum(int(i) for i in D )
    Inum=sum(int(i) for i in I )
    Hnum=sum(int(i) for i in H )
    Snum=sum(int(i) for i in S )
    #print(Mnum,D,I,H,S)
    return Mnum,Dnum,Inum,Hnum,Snum

def get_name_list(bwa_file,output):
    name_list=[]
    samfile=pysam.AlignmentFile(bwa_file,"rb")
    # read the bam file
    for read in samfile.fetch(until_eof=True):
        #print(read)
        #deal with each hit
        if read.is_unmapped:
            name_list.append(read.query_name)
        else:
            #read is mapped , get M number in cigar, if M/sequence length is over <0.9   min 10% difference
            M,D,H,I,S = get_cigarstring(read.cigarstring)
            if M/(M+D+H+I+S)<0.9:
                name_list.append(read.query_name)
    #output identifier of single copy reads    
    outfile = open(output,"w+")
    for i in name_list:
        outfile.write(i+"\n")
    outfile.close()
    return NONE

input_list=[bwa_file_f,bwa_file_r]
output_list=[output_f,output_r]
from multiprocessing import  Process
process_list = []
for i in range(2):  #开启5个子进程执行fun1函数
    p = Process(target=get_name_list,args=(input_list[i],output_list[i],)) #实例化进程对象
    p.start()
    process_list.append(p)

for i in process_list:
    p.join()
