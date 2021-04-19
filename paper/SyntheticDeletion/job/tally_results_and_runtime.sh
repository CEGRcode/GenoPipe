# Confirm all simulations are finished
ls -l results/sacCer3_Rap1_500K/ID/*.tab | wc -l 
ls -l results/sacCer3_Rap1_1M/ID/*.tab | wc -l 
ls -l results/sacCer3_Rap1_2M/ID/*.tab | wc -l 
ls -l results/sacCer3_Rap1_3M/ID/*.tab | wc -l 
ls -l results/sacCer3_Rap1_4M/ID/*.tab | wc -l 
ls -l results/sacCer3_Rap1_5M/ID/*.tab | wc -l 

ls -l results/sacCer3_Reb1_500K/ID/*.tab | wc -l 
ls -l results/sacCer3_Reb1_1M/ID/*.tab | wc -l 
ls -l results/sacCer3_Reb1_2M/ID/*.tab | wc -l 
ls -l results/sacCer3_Reb1_3M/ID/*.tab | wc -l 
ls -l results/sacCer3_Reb1_4M/ID/*.tab | wc -l 
ls -l results/sacCer3_Reb1_5M/ID/*.tab | wc -l 

# Get Success Tallies
CHECK=scripts/check_ID_tally.sh
bash $CHECK RAP1 sacCer3_Rap1_500K  > results/RAP1_500K_tally.txt
bash $CHECK RAP1 sacCer3_Rap1_1M    > results/RAP1_1M_tally.txt
bash $CHECK RAP1 sacCer3_Rap1_2M    > results/RAP1_2M_tally.txt
bash $CHECK RAP1 sacCer3_Rap1_3M    > results/RAP1_3M_tally.txt
bash $CHECK RAP1 sacCer3_Rap1_4M    > results/RAP1_4M_tally.txt
bash $CHECK RAP1 sacCer3_Rap1_5M    > results/RAP1_5M_tally.txt

bash $CHECK REB1 sacCer3_Reb1_500K  > results/REB1_500K_tally.txt
bash $CHECK REB1 sacCer3_Reb1_1M  > results/REB1_1M_tally.txt
bash $CHECK REB1 sacCer3_Reb1_2M  > results/REB1_2M_tally.txt
bash $CHECK REB1 sacCer3_Reb1_3M  > results/REB1_3M_tally.txt
bash $CHECK REB1 sacCer3_Reb1_4M  > results/REB1_4M_tally.txt
bash $CHECK REB1 sacCer3_Reb1_5M    > results/REB1_5M_tally.txt


# Get Runtimes from logs
grep 'real' logs/depth.did.Reb1.500K.log.err-* > results/REB1_500K_runtime.txt
grep 'real' logs/depth.did.Reb1.1M.log.err-* > results/REB1_1M_runtime.txt
grep 'real' logs/depth.did.Reb1.2M.log.err-* > results/REB1_2M_runtime.txt
grep 'real' logs/depth.did.Reb1.3M.log.err-* > results/REB1_3M_runtime.txt
grep 'real' logs/depth.did.Reb1.4M.log.err-* > results/REB1_4M_runtime.txt
grep 'real' logs/depth.did.Reb1.5M.log.err-* > results/REB1_5M_runtime.txt

grep 'real' logs/depth.did.Rap1.500K.log.err-* > results/RAP1_500K_runtime.txt
grep 'real' logs/depth.did.Rap1.1M.log.err-* > results/RAP1_1M_runtime.txt
grep 'real' logs/depth.did.Rap1.2M.log.err-* > results/RAP1_2M_runtime.txt
grep 'real' logs/depth.did.Rap1.3M.log.err-* > results/RAP1_3M_runtime.txt
grep 'real' logs/depth.did.Rap1.4M.log.err-* > results/RAP1_4M_runtime.txt
grep 'real' logs/depth.did.Rap1.5M.log.err-* > results/RAP1_5M_runtime.txt

wc -l results/*runtime*
