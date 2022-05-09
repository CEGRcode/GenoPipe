#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l pmem=24gb
#PBS -l walltime=03:00:00
#PBS -A open
#PBS -o logs/tally.log.out
#PBS -e logs/tally.log.err

module load anaconda3
source activate genopipe

WRK=/path/to/GenoPipe/paper/ENCODE-CellLines
cd $WRK

# Compile StrainID results
python scripts/analyze_encode_results.py -i results/ID/ -m 210512_sample_metadata.txt > results/encode_cell_line_results.txt 
