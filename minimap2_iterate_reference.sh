path=$(readlink -f "${BASH_SOURCE:-$0}")

DIR_PATH=$(dirname $path)
# minimap2 -d $DIR_PATH/ref.mmi $DIR_PATH/minimap2_reference.fa
# wait
# echo $@
# for arg in "$@";
# do
#     if [ $arg == "get_minimap2_iterate_reference.done" ]
#     then
#     continue
#     fi
#     name=$(basename $arg .fa)
# 	minimap2 -t 40 -a $DIR_PATH/ref.mmi $arg > $DIR_PATH/$name.diff_E28.sam
#     wait
#     samtools fasta -f4 -@40 $DIR_PATH/$name.diff_E28.sam > $DIR_PATH/$name.diff_E28.fa
#     wait
#     rm $DIR_PATH/$name.diff_E28.sam
#     wait
# done


file=$(ls -l $1*.diff_E28.fa | sort -k5nr | rev | cut -d " " -f1|rev | sed "s/.min_diff//g")
echo $file
#for arg in "$@";
for arg in $file;
do
    name=$(basename $arg .diff_E28.fa)
    seqkit stat $DIR_PATH/minimap2_reference.fa >> $DIR_PATH/minmap.record
    echo $name >> $DIR_PATH/minmap.record
    wait
	minimap2 -t 40 -a $DIR_PATH/minimap2_reference.fa $DIR_PATH/$name.fa > $DIR_PATH/$name.diff.sam
    wait
    samtools fasta -f4 -@40 $DIR_PATH/$name.diff.sam > $DIR_PATH/$name.min_diff_to_ref.fa
    wait
    rm $DIR_PATH/$name.diff.sam
    cat $DIR_PATH/$name.min_diff_to_ref.fa >> $DIR_PATH/minimap2_reference.fa
    wait
done
seqkit stat $DIR_PATH/minimap2_reference.fa >> $DIR_PATH/minmap.record
wait
mv $DIR_PATH/minimap2_reference.fa $DIR_PATH/minimap2_finall.fa
wait
grep "reference.fa  FASTA" $DIR_PATH/minmap.record | sed "s/,//g" > $DIR_PATH/minmap.record.hist
wait
Rscript /localdisk/home/s2224743/data/snakemake_pipeline/scripts/draw_histogram.R $DIR_PATH/minmap.record.hist
wait
