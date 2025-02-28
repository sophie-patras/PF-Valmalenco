#!/bin/bash

#cd ~
#spack load parflow
# echo "Using PARFLOW_DIR: ${PARFLOW_DIR}"

export EXEC_DIR=/home/patras/PF-Valmalenco/scripts/run/
cd ${EXEC_DIR}
# mkdir Tmp

export WORK_DIR=/home/patras/PF-Valmalenco/scripts/run/tmp/
# echo "Using WORK_DIR: ${WORK_DIR}"
# input files must be in Tmp

# cp /home/patras/Valmalenco/Data/DataPF/slope* ${WORK_DIR}

run_name=BB.c500.VMRF0IC-0375 # in {PLT,CLM,RR,Box...}
cellsize=500 # in {500; 250; 100}

## RUN SIMULATION
echo "Job started at " `date`

# chmod +rx pfsimulation_yaml.py
python3 pfsimulation_yamlall.py $run_name $cellsize --parflow-directory ${PARFLOW_DIR} --working-directory ${WORK_DIR} --show-line-error --validation-verbose -p 3 -q 3 -r 1

# run_script.py [-h] [--parflow-directory PARFLOW_DIRECTORY] [--parflow-version PARFLOW_VERSION]
#  [--working-directory WORKING_DIRECTORY] [--skip-validation] [--dry-run] [--show-line-error]
#  [--exit-on-error] [--write-yaml] [--validation-verbose] [-p P] [-q Q] [-r R]

echo "Job finished at " `date`

rm tmp/${run_name}.*.dist
#ls tmp/
#cp CLM_V53.yaml tmp/
#rm /home/patras/Valmalenco/Tmp/${run_name}_new/*
mkdir /home/patras/PF-Valmalenco/outputs/${run_name}_new/
cp /home/patras/PF-Valmalenco/scripts/run/tmp/${run_name}.* /home/patras/PF-Valmalenco/outputs/${run_name}_new/ 
cp InputKeys_c500_L11UX.yaml /home/patras/PF-Valmalenco/outputs/${run_name}_new/ 
cp tmp/slope*.c500.v4.pfb /home/patras/PF-Valmalenco/outputs/${run_name}_new/ 
rm /home/patras/PF-Valmalenco/scripts/run/tmp/${run_name}.*
ls /home/patras/PF-Valmalenco/outputs/${run_name}_new/
