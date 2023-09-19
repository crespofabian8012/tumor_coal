#!/bin/bash

#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -C clk
#SBATCH -c 48     #Máximo de nodos
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --mail-user crespofabian8012@gmail.com
#SBATCH --mail-type FAIL,BEGIN,END
#SBATCH -t 05:59:00
#SBATCH  --mem 144GB

 
module load cesga/2022

Rscript 1run_model_1_population_parallel_asc_16.R   $1









