WRK=/path/to/GenoPipe/paper/ENCODEdata-eGFP
cd $WRK

module load anaconda3
source activate genopipe
#module load python/3.6.8
python -V

# Get success rate based on metadata file and results of EpitopeID
METADATA=210227_sample_metadata.txt
CALCULATE=scripts/analyze_encode_results.py
python $CALCULATE -m $METADATA -i results/ID \
	> results/eGFPandFLAG_results_analyzed.txt \
	2> results/eGFPandFLAG_results_analyzed.err

