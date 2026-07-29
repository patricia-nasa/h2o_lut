#!/bin/bash
#SBATCH --job-name=generate_lut_array
#SBATCH --output=logs/lut_job_%A_%a.out  # %A is the job array ID, %a is the task ID
#SBATCH --error=logs/lut_job_%A_%a.err
#SBATCH --array=0-99                # 100 jobs, indexed 0 through 99
#SBATCH --time=12:00:00             # Adjust time limit as needed
#SBATCH --ntasks=1                  # 1 task per job
#SBATCH --cpus-per-task=1           # 1 CPU per task (adjust if python script is multithreaded)
#SBATCH --mem=4G                    # Memory per job (adjust as needed)
#SBATCH --account=s2942
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug

# Load your python environment here if needed (e.g., module load python, conda activate env, etc.)
unalias conda
source /gpfsm/dnb34/pcastell/miniforge3/etc/profile.d/conda.csh
conda activate h2o-lut


# The current array task ID becomes the start channel
START_CHANNEL=$SLURM_ARRAY_TASK_ID

# python script treats "stop" as exclusive (e.g., start 0 stop 1 processes only channel 0), add 1:
STOP_CHANNEL=$((START_CHANNEL + 1))

echo "Starting task $SLURM_ARRAY_TASK_ID"
echo "Running channels: --channel-start $START_CHANNEL --channel-stop $STOP_CHANNEL"

# Execute the python script
python generate_lut.py config.yaml --channel-start $START_CHANNEL --channel-stop $STOP_CHANNEL

echo "Task $SLURM_ARRAY_TASK_ID completed."
