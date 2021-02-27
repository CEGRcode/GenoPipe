WRK=/path/to/GenoPipe/paper/ENCODEdata-eGFP
cd $WRK

## Get success rate based on metadata file and results of EpitopeID
METADATA=210227_sample_metadata.txt
CALCULATE=scripts/analyze_encode_results.py
python $CALCULATE -m $METADATA -i results/ID \
	> results/eGFP_results_analyzed.txt \
	2> results/eGFP_results_analyzed.err

