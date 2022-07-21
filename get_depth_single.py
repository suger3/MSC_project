#input is by bash pipe with sam file such as 
#samtools view -@ 10 -h -F 1796 test_sorted.bam | python get_depth_single.py
import pysam
import sys
#bwa_file=sys.argv[1]
#print(bwa_file)
samfile=pysam.AlignmentFile("-", "r")
header = dict(samfile.header)
dict={}
for i in header['SQ']:
    dict[i['SN']]={'LN':int(i['LN']),'base':0,'total_length':0}
for read in samfile.fetch():
    base=0
    total_length=0
    #print(read)
    contigs=read.reference_name
    a=list(read.cigar)
    total_length=int(read.query_length)
    for i in a:
        #print(i)
        tmp=list(i)
        #print(tmp)
        if tmp[0]==0:
            base+=int(tmp[1])
        #if tmp[0]==1:
            #base+=int(tmp[1])
        if tmp[0]==2:
            base+=int(tmp[1])
    # if dict.__contains__(contigs):
    dict[contigs]["base"]=dict[contigs]["base"]+base
    dict[contigs]["total_length"]=dict[contigs]["total_length"]+total_length
    # else:
    #     dict[contigs]={}
    #     dict[contigs]["base"]=base
    #     dict[contigs]["total_length"]=total_length
    
    #if "NODE_3" in contigs:
        #break
print("Contigs\tlength\tmapped_base\tmapped_read_total_length\tdepth")
for key in dict:
    print("{}\t{}\t{}\t{}\t{:.2f}".format(key,dict[key]['LN'],dict[key]['base'],dict[key]['total_length'],dict[key]['base']/int(dict[key]['LN'])))
    #if "NODE_1_" in key:
        #print("{}\t{}\t{}\t{}\t{:.2f}".format(key,dict[key]['LN'],dict[key]['base'],dict[key]['total_length'],dict[key]['base']/int(dict[key]['LN'])))
    #if "NODE_2_" in key:
        #print("{}\t{}\t{}\t{}\t{:.2f}".format(key,dict[key]['LN'],dict[key]['base'],dict[key]['total_length'],dict[key]['base']/int(dict[key]['LN'])))
#print("{}\t{}\t{}".format(name,base,total_length))
#print ("{}\t{}\t{}".format(bwa_file,,base))
# samtools view -H test_sorted.bam >at split/NODE_all_title
# samtools view -H ../test_sorted.bam | cut -f2 > NODE_title
# cat NODE_all_title | cat - $1 > temp && mv temp $1
