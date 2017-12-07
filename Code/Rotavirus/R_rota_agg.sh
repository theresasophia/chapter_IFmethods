#!/bin/bash 

# Set the allocation to be charged for this job
# not required if you have set a default allocation

# The name of the script is myjob
#SBATCH -J R_rota_agg

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 20:00:00

# set the project to be charged for this
# should normally be of the format 2015-1 or 2015-16-1 or similar
#SBATCH -A 2017-1

# Number of nodes
#SBATCH --nodes=1
# Number of MPI processes per node (24 is recommended for most cases)
# 48 is the default to allow the possibility of hyperthreading
#SBATCH --ntasks-per-node=24
# Number of MPI processes.

#SBATCH -e error_file.e
#SBATCH -o output_file.o
#SBATCH --mail-type=ALL

# load the R module
module add R
 
# Run the executable named myexe 
# and write the output into my_output_file
srun -n 1 Rscript rota_agg.R
