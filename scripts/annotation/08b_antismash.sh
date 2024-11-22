#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 16G --out logs/antismash.update.%a.log -J collet_antismash
#SBATCH --array 1

module load antismash/6.1.1
which perl
which antismash
hostname
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
BUSCO=fungi_odb10
INDIR=genomes
OUTDIR=annotate
SAMPFILE=genomes/samples.csv

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
  N=$1
  if [ ! $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
  fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi

INPUTFOLDER=predict_results

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SPECIES STRAIN PHYLUM BIOSAMPLE BIOPROJECT LOCUSTAG
do
  BASE=$(echo -n "$SPECIES" | perl -p -e 's/\s+/_/g')
  name=$BASE
  STRAIN_NOSPACE=$(echo -n "$STRAIN" | perl -p -e 's/\s+/_/g')
  echo "$BASE"
  MASKED=$(realpath $INDIR/$BASE.masked.fasta)

  if [ ! -d $OUTDIR/$name ]; then
    echo "No annotation dir for ${name}"
    exit
  fi
  echo "processing $OUTDIR/$name"
  if [[ ! -d $OUTDIR/$name/antismash_local && ! -s $OUTDIR/$name/antismash_local/index.html ]]; then
    #	antismash --taxon fungi --output-dir $OUTDIR/$name/antismash_local  --genefinding-tool none \
      #    --asf --fullhmmer --cassis --clusterhmmer --asf --cb-general --pfam2go --cb-subclusters --cb-knownclusters -c $CPU \
      #    $OUTDIR/$name/$INPUTFOLDER/*.gbk
    time antismash --taxon fungi --output-dir $OUTDIR/$name/antismash_local \
      --genefinding-tool none --fullhmmer --clusterhmmer --cb-general \
      --pfam2go -c $CPU $OUTDIR/$name/$INPUTFOLDER/*.gbk
  fi
done
