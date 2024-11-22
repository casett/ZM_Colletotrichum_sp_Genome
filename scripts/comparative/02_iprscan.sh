#!/bin/bash
#SBATCH -p intel --ntasks 16 --mem 98G --out logs/interpro.log
#SBATCH -J sg_colleto_interpro


#module load iprscan
module load interproscan/5.55-88.0

mkdir interpro_out

interproscan.sh -i all_proteins.fa --goterms -f TSV -d interpro_out --disable-precalc -t p -cpu 16


#kinfin/scripts/iprs2table.py -i all_proteins.tsv
