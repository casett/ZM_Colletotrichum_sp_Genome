#!/usr/bin/bash -l
#SBATCH --ntasks 24 --mem 16G --time 6:00:00 -p intel -N 1 --out logs/phyling.%A.log
#SBATCH -J Colleto_Phyling

module unload miniconda2
module load miniconda3
module load hmmer/3
module unload perl
module load parallel
module load muscle

if [ ! -f config.txt ]; then
	echo "Need config.txt for PHYling"
	exit
fi

source config.txt
if [ ! -z $HMM ]; then
	rm -rf aln/$HMM
fi

./PHYling_unified/PHYling init
./PHYling_unified/PHYling search -q parallel
./PHYling_unified/PHYling aln -c -q parallel


