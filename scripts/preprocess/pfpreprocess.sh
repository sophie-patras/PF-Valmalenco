#!/bin/bash

cd ~
spack load parflow
# echo "Using PARFLOW_DIR: ${PARFLOW_DIR}"

export EXEC_DIR=/home/patras/Valmalenco/Codes/PreProcess/
cd ${EXEC_DIR}
# mkdir Tmp

# export WORK_DIR=/home/patras/Valmalenco/Codes/PreProcess/
# echo "Using WORK_DIR: ${WORK_DIR}"

# cp /home/patras/Valmalenco/Data/DataPF/slope* ${WORK_DIR}

run_name=PLT
cellsize=500

## RUN SIMULATION

python3 ascdem_to_pfsol.py $run_name -p 2 -q 2 -r 1
# $run_name $cellsize --parflow-directory ${PARFLOW_DIR} --working-directory ${WORK_DIR} --show-line-error --validation-verbose -p 2 -q 2 -r 1

