#!/bin/bash
#SBATCH --job-name=nufft_1D
#SBATCH --output=out/out_nufft_1D.out
#SBATCH --error=err/err_nufft_1D.err
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --mem=1500000
#################

module load matlab
cd ../../test/nufft
matlab -nojvm -r "batch_nufft_1D;quit"
