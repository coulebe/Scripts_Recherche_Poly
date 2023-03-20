#!/bin/bash
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=1
#SBATCH --time=02-00:00:00 # DD-HH:MM    
#SBATCH --mail-user=<alfred-gaetan.coulebetouba@polymtl.ca>
#SBATCH --mail-type=ALL

cd Python/MoD/


source py38/bin/activate
module load python/3.8.13 scipy-stack

python MoD1_Discovery_CC_Norm.py
