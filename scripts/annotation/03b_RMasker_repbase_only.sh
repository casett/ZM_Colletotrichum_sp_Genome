#!/usr/bin/bash
#SBATCH -p stajichlab
#SBATCH -N 1 
#SBATCH -n 64 
#SBATCH --mem 96gb 
#SBATCH --out logs/rmasker_repbase.log
#SBATCH -J Colleto_rmasker_rebase

DIR=genomes
PREFIX=SG_Colletotrichum

module load RepeatMasker/4-1-1
#module load funannotate
#module unload ncbi-rmblast
#module load ncbi-rmblast/2.9.0-p2
#module unload miniconda2
#module load miniconda3


RepeatMasker -pa 64 -s -e ncbi -species Ascomycota ${DIR}/$PREFIX.sorted.fasta > ${DIR}/$PREFIX.RM.out

