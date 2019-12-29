#!/bin/sh
#SBATCH --partition=debug
#SBATCH --ntasks=24
#SBATCH --tasks-per-node=24

echo partition: $SLURM_JOB_PARTITION
echo num_nodes: $SLURM_JOB_NUM_NODES nodes: $SLURM_JOB_NODELIST
echo num_tasks: $SLURM_NTASKS tasks_node: $SLURM_NTASKS_PER_NODE

ORDERED_SWITCH=1

# activate environment
module load lang/Anaconda3/2.5.0
source activate pygem_hpc

ROI="HMA"
# region batch string
latlon_batch_str="${ROI}_latlon_batch"

# split glaciers into batches for different nodes
python spc_split_lists.py -n_batches=$SLURM_JOB_NUM_NODES -option_ordered=$ORDERED_SWITCH

# list  batch filenames
latlon_fns=$(find ${latlon_batch_str}*)
echo latlon filenames:$latlon_fns
# create list
list_latlon_fns=($latlon_fns)
echo first_batch:${list_latlon_fns[0]}


for i in $latlon_fns 
do
  # print the filename
  echo $i
  
  # determine batch number
  BATCHNO="$(cut -d'.' -f1 <<<$(cut -d'_' -f6 <<<"$i"))"
  echo $BATCHNO
  
  # run the file on a separate node (& tells the command to move to the next loop for any empty nodes)
  srun -N 1 -n 1 python meltmodel_global.py -num_simultaneous_processes=$SLURM_NTASKS_PER_NODE -latlon_fn=$i&
done
# wait tells the loop to not move on until all the srun commands are completed
wait

echo -e "\nScript finished"
