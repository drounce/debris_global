#!/bin/sh
#SBATCH --partition=debug
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24

echo partition: $SLURM_JOB_PARTITION
echo num_nodes: $SLURM_JOB_NUM_NODES nodes: $SLURM_JOB_NODELIST
echo num_tasks: $SLURM_NTASKS tasks_node: $SLURM_NTASKS_PER_NODE

# activate environment
module load lang/Anaconda3/2.5.0
source activate pygem_hpc
  
# run the file on a separate node (& tells the command to move to the next loop for any empty nodes)
srun -N 1 -n 1 python ERA5_preprocess_v2.py


echo -e "\nScript finished"
