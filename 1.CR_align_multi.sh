#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 17
#SBATCH -w gattaca
#SBATCH -t 5800-12
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

#Author Leanne Whitmore run multi 

export PATH="cellranger-9.0.1/bin:$PATH"

cellranger multi --id "O10_RUp16_Wk80-81" --csv ./Multifiles/MultiConfig_RUp16-wks-80-81.csv
wait 

cellranger multi --id "010_RUp16_weeks_18-20-48" --csv ./Multifiles/MultiConfig_RUp16-wks-18-20-48.csv
wait 

# Errored unless feature type was changed to VDJ B, and chemistry was specified
cellranger multi --id "010_RUp16_weeks_26-36-48" --csv ./Multifiles/MultiConfig_RUp16-wks-26-36-48.csv
wait 

cellranger multi --id "010_RUp16_Wk87" --csv ./Multifiles/MultiConfig_RUp16-wks-87.csv
wait 

cellranger multi --id "010-RUp16-Efficacy_Wk83LN_Rsf16_Wk50" --csv ./Multifiles/MultiConfig_RUp16-Efficacy_Wk83LN_Rsf16_Wk50.csv
wait 

cellranger multi --id "010-RUp16-Efficacy_Wk93" --csv ./Multifiles/MultiConfig_RUp16-Efficacy-Wk93.csv
wait 

#Run2
cellranger multi --id "RUp16_PBMC_wk26_36_RSf16_PBMC_wk48" --csv ./Multifiles/MultiConfig_RUp16_PBMC_wk26_36_RSf16_PBMC_wk48.csv
wait 
