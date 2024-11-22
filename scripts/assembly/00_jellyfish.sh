#!/usr/bin/bash
#SBATCH -p intel,batch
#SBATCH --nodes=1
#SBATCH --ntasks=16 # Number of cores
#SBATCH --mem=48G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -J Colleto_jellyfish
#SBATCH -o logs/jellyfishlog
#SBATCH -e logs/jellyfish.log



module load jellyfish/2.3.0
CPU=16
SAMPFILE=samples.csv
PREFIX=SG_Colletotrichum
LOCATION=raw_seqs
FREAD=SG_Colletotrichum_R1.fq.gz
RREAD=SG_Colletotrichum_R2.fq.gz

jellyfish count -C -m 21 -s 2G -t 16 <(zcat ${LOCATION}/$FREAD) <(zcat ${LOCATION}/$RREAD) -o $PREFIX.jf

jellyfish histo -t 16 $PREFIX.jf > $PREFIX.jf.histo


