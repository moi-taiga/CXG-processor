#!/bin/bash


#SBATCH --nodes=2
#SBATCH --mem=200G 
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --partition=amd
#SBATCH --export=ALL




module load apptainer
module load cuda
apptainer run -H $PWD \
--nv \
--bind "$PWD:/workspace" \
--bind /iridisfs/ddnb/Ahmed/data:/iridisfs/ddnb/Ahmed/data \
/iridisfs/ddnb/images/bionemo-framework_nightly \
python $1


