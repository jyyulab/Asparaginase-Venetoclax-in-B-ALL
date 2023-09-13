module load cellranger-arc/2.0.0

i="2828"
bsub -R "rusage[mem=3000]" -P cellRanger -J cellRanger_${i} -oo ${i}.out -eo ${i}.err  -q compbio -n 15 "cellranger-arc count --id=${i} \
                   --reference=./cellranger_reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                   --libraries=./libraries_${i}.csv \
                   --localcores=16 \
                   --localmem=64 \
                   --jobmode=lsf --maxjobs=200"
