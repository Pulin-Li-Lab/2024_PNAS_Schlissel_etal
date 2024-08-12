#!/bin/bash
## Configuration values for SLURM job submission.
## One leading hash ahead of the word SBATCH is not a comment, but two are.
#SBATCH --job-name=friendlyname        # friendly name for job.
#SBATCH --nodes=1                      # ensure cores are on one node
#SBATCH --ntasks=1                 # run a single task
#SBATCH --cpus-per-task=2             # number of cores/threads requested.
#SBATCH --mem=8gb                      # memory requested.
#SBATCH --partition=20                 # partition (queue) to use
#SBATCH --array=1-50
#SBATCH --output friendlyname-%j.out   # name of output file.  %j is jobid


echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

NOW=$( date +%s )

MAP_FOLD=50
#CELL_SIZE=3.5355
CELL_AREA=25
CELL_SIDES=4
FRAC_FAST=0.05
N_AGENTS=2000
MAX_TIME=20000

SLOW_RATE=0.5
FAST_RATE=1


for FRAC in $(seq 1 1 $SLURM_ARRAY_TASK_ID)
do
    sleep 5
done

CONSTRAINT=0
srun --nodes=1 -n1 --mem=4G matlab -nodisplay -nosplash -nodesktop -r "pid=$NOW;n_agents=$N_AGENTS;max_time=$MAX_TIME;slowRate=$SLOW_RATE;fastRate=$FAST_RATE;map_frac=$SLURM_ARRAY_TASK_ID;map_fold=$MAP_FOLD;frac_fast=$FRAC_FAST;cell_sides=$CELL_SIDES;cell_area=$CELL_AREA;constraint=$CONSTRAINT;run('/PATH/TO/CODE/DIRECTORY/steady_state_square_bilayer_toy_grid/steady_state_toy.m');exit;" &

sleep 2

CONSTRAINT=1
srun --nodes=1 -n1 --mem=4G matlab -nodisplay -nosplash -nodesktop -r "pid=$NOW;n_agents=$N_AGENTS;max_time=$MAX_TIME;slowRate=$SLOW_RATE;fastRate=$FAST_RATE;map_frac=$SLURM_ARRAY_TASK_ID;map_fold=$MAP_FOLD;frac_fast=$FRAC_FAST;cell_sides=$CELL_SIDES;cell_area=$CELL_AREA;constraint=$CONSTRAINT;run('/PATH/TO/CODE/DIRECTORY/steady_state_square_bilayer_toy_grid/steady_state_toy.m');exit;" &

echo finished calling steady_MSD --  waiting for them to finish
wait
echo script complete
echo
