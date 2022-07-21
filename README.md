# MSC_project
## k_mer.smk
A snakemake script to calcualte the kmer distribution for paired reads.
dependence: kmc,kmc_tools
input: paired reads
output: KMC database and histogram of kmer distribution. This histogram is used to draw kmer histo by Tetmer or GenomeScopes.
## mapping.smk
A snakemake script to get pan-genome.
dependence:kmc_tools,bwa-mem2,samtools,python3,seqtk,blastn,spades,cd-hit-est,seqkit,minimap2,python scripts.
input:paired reads,KMC database
output: single copy sequence assembled by spades, pan-genome, and report.
## depth.smk
A snakemake script to calcualte the depth for each contigs in pan-genome.
dependence:bwa-mem2,samtools,python script.
input:paired reads
output:tsv file
## get_depth_single.py
A python script to calculate the depth for each contigs.
input: sam file in pipe
output: tsv file
example:
  samtools view -@ 10 -h -F 1796 test_sorted.bam | python get_depth_single.py - > my_result.tsv
or
  bwa-mem2 mem -t 10  E28_assembly_over200.fa E069_f.fq.gz E069_r.fq.gz | samtools view  -@ 5 -h -F 1796 | python get_depth_single.py - > my_result.tsv
 ## get_unmap_read.py
 A python script to filter extract the name of unmapped or poor mapped reads.
 input:bam file for forward reads and reverse reads(seprately)
 output:two name list
 ## link_short_contigs.sh
 A bash script to link the contigs shorter than threshold(500 in example), and add 1000 "N" between them to generate a huge contigs.
 dependence:seqkit,perl
 input:fasta file, threshold
 output:contigs_link_{threshold}.fa
 example:
  bash link_short_contigs.sh pan_genome_35.fa.masked 500
 ##
