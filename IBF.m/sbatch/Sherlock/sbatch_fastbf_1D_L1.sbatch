#!/bin/bash
#SBATCH --job-name=fastbf_1D_L1
#SBATCH --output=out/out_fastbf_1D_L1.out
#SBATCH --error=err/err_fastbf_1D_L1.err
#SBATCH --time=16:00:00
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --mem=1500000
#################

module load matlab
cd ../../test/fio
matlab -nojvm -r "batch_fastbf_1D_L1;quit"
