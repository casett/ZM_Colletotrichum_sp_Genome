#!/usr/bin/bash
#SBATCH --ntasks 24 --mem 16G --time 2:00:00 -p short -N 1 --out logs/fastANI.%A.log
#SBATCH -J Colleto_ANI


INDIR=data
OUTPUT=ANI_out

/rhome/cassande/bigdata/software/FastANI/fastANI -q $INDIR/Colletotrichum.sorted.fasta -r $INDIR/GCA_001563125.1_CSAL01_genomic.fna -o $OUTPUT
