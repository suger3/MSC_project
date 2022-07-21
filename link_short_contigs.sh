#$1 is input fasta
#$2 is the length of link contigs
DIR_PATH=$(dirname $1)
name=$(basename $1)
gg=$2
seqkit seq -g -M $2 $1 > $DIR_PATH/$name.need_link_file
wait
seqkit seq -g -m $((gg+1)) $1 > $DIR_PATH/contigs_link_$2.fa
wait
echo ">link_contigs" >> $DIR_PATH/contigs_link_$2.fa
wait
perl -pe '@ones ="N"x1000;s/^>.*/@ones/gm' $DIR_PATH/$name.need_link_file | sed ':a;N;$!ba;s/\n//g' >>$DIR_PATH/contigs_link_$2.fa
wait
rm $DIR_PATH/$name.need_link_file

