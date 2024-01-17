#!/bin/bash
#SBATCH --job-name=PGS_project
#SBATCH --output=PGS2.out
#SBATCH --mail-user=h.van.der.pol@umail.leidenuniv.nl
#SBATCH --mail-type=ALL
#SBATCH --partition=cpu-medium,cpu-long,mem
#SBATCH --time=23:45:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --nodes=1

module load R/4.0.5-foss-2020b
Rscript Compare_E_vs_UI_v3.R