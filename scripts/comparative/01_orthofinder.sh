#!/bin/bash -l
#SBATCH -p highmem --ntasks 24 --mem 750G --out logs/orthofinder.cont.log
#SBATCH -J sg_colleto_orthofinder
##SBATCH -t 14-00:00:00

module load diamond
module load orthofinder


CPU=24

ulimit -n 262144
ulimit -Hn
ulimit -Sn


#orthofinder -f data -M msa -A mafft -T fasttree -o results -t $CPU
orthofinder -M msa -A mafft -T fasttree -ft results/Results_Sep06 -t $CPU -a $CPU

