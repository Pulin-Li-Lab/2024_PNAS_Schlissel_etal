#!/bin/bash
## Configuration values for SLURM job submission.
## One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=friendlyname        # friendly name for job.
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                 # run a single task
#SBATCH --cpus-per-task=12             # number of cores/threads requested.
#SBATCH --mem=24gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
#SBATCH --array=1-9
#SBATCH --output friendlyname-%j.out   # name of output file.  %j is jobid


echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

SENDER_CONC=1000
MAX_TIME=10000
MAP_FOLD=9
FRAC_FAST=0.05

slowRate=0.5
fastRate=20
CELL_SIZE=2

echo
echo "Sender conc is $SENDER_CONC"
echo "Max time is $MAX_TIME"
echo "Frac fast is $FRAC_FAST"
echo "Rate vector is [ $slowRate $fastRate ]"
echo "Cell size is $CELL_SIZE"


for FRAC in $(seq 1 1 $SLURM_ARRAY_TASK_ID)
do
    sleep 10
done
    
srun --nodes=1 -n1 --mem=16G matlab -nodisplay -nosplash -nodesktop -r "sender_conc=$SENDER_CONC;max_time=$MAX_TIME;map_frac=$SLURM_ARRAY_TASK_ID;map_fold=$MAP_FOLD;frac_fast=$FRAC_FAST;slowRate=$slowRate;fastRate=$fastRate;cell_dimension=$CELL_SIZE;run('/PATH/TO/ANALYSIS/DIRECTORY/steady_membrane_toy_grid/steady_state_toy.m');exit;" 


echo finished calling steady_FLUX --  waiting for them to finish
wait
echo script complete
echo



