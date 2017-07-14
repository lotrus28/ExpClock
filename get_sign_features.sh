#!/bin/bash
# GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct
# ./get_sign_features.sh GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct /mnt/lustre/fgalkin/gtex/tissue_samples/ ./subj_sample_annot.txt

rpkm_file=$1
tissue_folder=$2
metadata=$3
SCRIPT_FOLDER='/mnt/lustre/fgalkin/gtex'

mkdir $tissue_folder

echo "python3 subset_data.py $rpkm_file $tissue_folder" | qsub -N tissues_div -l nodes=node12

declare -A nodes
nodes+=( [0]="node01" [1]="node03" [2]="node07" [3]="node08" [4]="node12" [5]="node15" )
jobID=0
wait='tissues_div'
for f in $tissue_folder/*
do
	node=${nodes[$(($jobID % 6))]}
	echo 'R CMD BATCH ''"--args '"$f $metadata ./"'" '"$SCRIPT_FOLDER/trial.R" | qsub -N TisBor$jobID -l nodes=$node -hold_jid $wait
	jobID=$(($jobID + 1))
	sleep 1s
done

sleep 10s

echo 'R CMD BATCH ''"--args '"$tissue_folder ./subj_sample_annot.txt ./"'" '"$SCRIPT_FOLDER/test_rand_for.R" | qsub -N test_forest -l nodes=node12 -hold_jid TisBor$jobID