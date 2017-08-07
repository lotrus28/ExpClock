#!/bin/bash
# ./get_sign_features.sh GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct /mnt/lustre/fgalkin/gtex/tissue_samples/ ./subj_sample_annot.txt
# ./get_sign_features.sh none /data5/bio/runs-galkin/expclock/tissue_samples /data5/bio/runs-galkin/expclock/subj_sample_annot.txt

rpkm_file=$1
tissue_folder=$2
metadata=$3
# SCRIPT_FOLDER='/mnt/lustre/fgalkin/gtex/Scripts'
# module load python/python-3.5.1

SCRIPT_FOLDER='/data5/bio/runs-galkin/expclock/Scripts'

if [ ! -f $metadata ]; then
	python3 $SCRIPT_FOLDER/join_sample_and_subj_annot.py
fi

if [ ! -d $tissue_folder ]; then
	mkdir $tissue_folder
	echo "python3 $SCRIPT_FOLDER/subset_data.py $rpkm_file $tissue_folder" | qsub -N tissues_div -l nodes=node12
fi

declare -A nodes
# nodes+=( [0]="node01" [1]="node03" [2]="node07" [3]="node08" [4]="node12" [5]="node15" )
nodes+=( [0]="node6.net0.pyxis.ripcm.com" [1]="node8.net0.pyxis.ripcm.com" [2]="node9.net0.pyxis.ripcm.com" [3]="node4.net0.pyxis.ripcm.com" [4]="node11.net0.pyxis.ripcm.com" [5]="node12.net0.pyxis.ripcm.com" )
jobID=0
if [ ! -f $tissue_folder/significant_features.txt ]; then
	for f in $tissue_folder/*
	do
		echo "$f"
		node=${nodes[$(($jobID % 6))]}
		# echo 'R CMD BATCH ''"--args '"$f $metadata $tissue_folder"'" '"$SCRIPT_FOLDER/boruta_clock.R" | qsub -N TisBor$jobID -l nodes=$node
		echo 'R CMD BATCH ''"--args '"$f $metadata $tissue_folder"'" '"$SCRIPT_FOLDER/boruta_clock.R" | qsub -pe make 3 -N TisBor$jobID -cwd -l hostname=$node -r yes
		jobID=$(($jobID + 1))
		sleep 15m
	done
fi

last_job=$(($jobID - 1))

echo "Loop over"

# echo 'R CMD BATCH ''"--args '"$tissue_folder ./subj_sample_annot.txt ./"'" '"$SCRIPT_FOLDER/test_rand_for.R" | qsub -N test_forest -l nodes=node12 -hold_jid TisBor$jobID
echo 'R CMD BATCH ''"--args '"$tissue_folder ./subj_sample_annot.txt ./"'" '"$SCRIPT_FOLDER/test_rand_for.R" | qsub -N test_forest -hold_jid TisBor$last_job -cwd -l hostname=$node -r yes
echo 'Over'
