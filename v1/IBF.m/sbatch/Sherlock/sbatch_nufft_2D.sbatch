#!/bin/bash
#SBATCH --job-name=nufft_2D
#SBATCH --output=out/out_nufft_2D.out
#SBATCH --error=err/err_nufft_2D.err
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --mem=1500000
#################

module load matlab
cd ../../test/nufft
matlab -nojvm -r "batch_nufft_2D;quit"
