#!/bin/bash
##
#SBATCH -o logs/quast_ragtag.log.txt
#SBATCH -e logs/quast_ragtag.log.txt
#SBATCH --nodes=1
#SBATCH --ntasks=12 # Number of cores
#SBATCH --mem=60G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=24:00:00
#SBATCH -J colleto_ragtag_quast


module load busco

INDIR=genomes
SAMPFILE=genomes/samples.csv


#/rhome/cassande/bigdata/software/quast-5.0.2/quast.py $INDIR/Colletotrichum.masked.fasta --threads 12 --eukaryote --space-efficient --conserved-genes-finding
/rhome/cassande/bigdata/software/quast-5.1.0rc1/quast.py  $INDIR/Colletotrichum.masked.fasta --threads 12 --eukaryote --space-efficient --conserved-genes-finding

