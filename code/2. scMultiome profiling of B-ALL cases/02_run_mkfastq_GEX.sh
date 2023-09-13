#!/bin/env bash
#BSUB -P cellranger_arc
#BSUB -J mkfastq
#BSUB -oo cellranger_arc_mkfastq_log.out 
#BSUB -eo cellranger_arc_mkfastq_err.err 
#BSUB -R "rusage[mem=128000]" 
#BSUB -q large_mem
#BSUB -N xin.huang@stjude.org
#BSUB -n 1

indir=./Sequencing_Runs/Illumina/NovaSeq/221228_A00641_0539_AHYCJKDRX2
 module load cellranger-arc/2.0.0
 cellranger-arc mkfastq --id=2828 \
 --run=${indir} \
 --csv=SampleSheet.csv \
 --jobmode=lsf --maxjobs=200
