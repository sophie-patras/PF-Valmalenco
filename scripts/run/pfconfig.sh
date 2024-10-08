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

run_name=SeepBox_UZ # in {PLT,CLM,...}
cellsize=250 # in {500; 250; 100}

## RUN SIMULATION
echo "Job started at " `date`

# chmod +rx pfsimulation_yaml.py
python3 pfsimulation_yamlall.py $run_name $cellsize --parflow-directory ${PARFLOW_DIR} --working-directory ${WORK_DIR} --show-line-error --validation-verbose -p 4 -q 4 -r 1

# run_script.py [-h] [--parflow-directory PARFLOW_DIRECTORY] [--parflow-version PARFLOW_VERSION]
#  [--working-directory WORKING_DIRECTORY] [--skip-validation] [--dry-run] [--show-line-error]
#  [--exit-on-error] [--write-yaml] [--validation-verbose] [-p P] [-q Q] [-r R]

echo "Job finished at " `date`

rm tmp/*.dist
ls tmp/
cp InputKeys_c250.yaml tmp/${run_name}.yaml
#rm /home/patras/Valmalenco/Tmp/${run_name}_new/*
mkdir /home/patras/PF-Valmalenco/outputs/${run_name}_new/
mv /home/patras/PF-Valmalenco/scripts/run/tmp/${run_name}.* /home/patras/PF-Valmalenco/outputs/${run_name}_new/ 

