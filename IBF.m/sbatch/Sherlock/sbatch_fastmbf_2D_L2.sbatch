#!/bin/bash
#SBATCH --job-name=fastmbf_2D_L2
#SBATCH --output=out/out_fastmbf_2D_L2.out
#SBATCH --error=err/err_fastmbf_2D_L2.err
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --mem=1500000
#################

module load matlab
cd ../../test/fio
matlab -nojvm -r "batch_fastmbf_2D_L2;quit"
