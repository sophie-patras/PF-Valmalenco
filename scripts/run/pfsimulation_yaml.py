# -*- coding: utf-8 -*-

"""
Python Parflow simulation, using yaml for input keys

Created on Mon Aug 5 11:37:26 2024
@author: Sophie Patras
"""

import os, sys, argparse

from parflow import Run
from parflow.tools.io import write_pfb, read_pfb #read_clm
from parflow.tools.fs import cp, mkdir, chdir, get_absolute_path, rm
from parflow.tools.settings import set_working_directory
from parflow.tools.compare import pf_test_file

# set tcl_precision 17

parser = argparse.ArgumentParser()
parser.add_argument('run_name', type = str, default='test')
parser.add_argument('cellsize', type = float, default=500.0)
# parser.add_argument('-parflow-directory', '--parflow-directory')
parser.add_argument('-working-directory', '--working-directory', default="~/Valmalenco/Codes/RunProcess/Tmp")
parser.add_argument('--show-line-error', action='store_true')
parser.add_argument('--validation-verbose', action='store_true')
parser.add_argument('--write-yaml',action='store_true')
# parser.add_argument('-p', '--p', type=int, default=2)
# parser.add_argument('-q', '--q', type=int, default=2)
# parser.add_argument('-r', '--r', type=int, default=1)
args = parser.parse_args()
# call argument : es Topology.P = args.p

runname = Run(args.run_name, __file__)

#-----------------------------------------------------------------------------
# File input version number
runname.FileVersion = 4
#-----------------------------------------------------------------------------
# Processor topology 
runname.pfset(yaml_file='Process.yaml')
#-----------------------------------------------------------------------------
# Computational Grid
runname.pfset(yaml_file='ComputationalGrid.yaml')
#-----------------------------------------------------------------------------
# Slopes
runname.pfset(yaml_file='TopoSlopes.yaml')
# runname.dist('slopeX.c500.v1.pfb')
# runname.dist('slopeY.c500.v1.pfb')

#-----------------------------------------------------------------------------
# Geometrical definitions
runname.pfset(yaml_file='GeomInput.yaml')
runname.pfset(yaml_file='Domain.yaml')

#-----------------------------------------------------------------------------
# Time
runname.pfset(yaml_file='TimingInfo.yaml')
runname.pfset(yaml_file='TimeStep.yaml')
runname.pfset(yaml_file='Gravity.yaml')
# RF Forcing period
runname.pfset(yaml_file='Cycle.yaml')

runname.pfset(yaml_file='Phase.yaml')
runname.pfset(yaml_file='PhaseSources.yaml')

#-----------------------------------------------------------------------------
# Hydraulics
runname.pfset(yaml_file='Mannings.yaml')
runname.pfset(yaml_file='Perm.yaml')
runname.pfset(yaml_file='SpecificStorage.yaml')

runname.pfset(yaml_file='Contaminants.yaml')
runname.pfset(yaml_file='Wells.yaml')

#-----------------------------------------------------------------------------
# Boundary Conditions (OverlandKinematic)
runname.pfset(yaml_file='BCPressure.yaml')
runname.pfset(yaml_file='Patch.yaml')
runname.pfset(yaml_file='ICPressure.yaml')

#-----------------------------------------------------------------------------
# Recap of physical attributes positioned on geometry
runname.pfset(yaml_file='Geom.yaml')

#-----------------------------------------------------------------------------
# Solver
runname.pfset(yaml_file='Solver.yaml')
#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
runname.KnownSolution = "NoKnownSolution"

#-----------------------------------------------------------------------------
# Distribute inputs, if needs different topology
#-----------------------------------------------------------------------------
#runname.dist("slopesX.c"+str(cellsize)+".v1.pfb")
#runname.dist("slopesY.c"+str(cellsize)+".v1.pfb")
#runname.dist("solid.c500.v1.pfsol")

#-----------------------------------------------------------------------------
# Run and unload pf output files 
#-----------------------------------------------------------------------------

# runname.Run.from_definition(runscript_path)
# print(f"Loaded run with runname: {run.get_name()}")

# runname.validate()

#output_dir = get_absolute_path('./Tmp/')

# output_dir = os.path.join(args.working-directory,"Tmp")
# set_working_directory(f"{output_dir}")

# runname.run(working_directory=output_dir)
runname.run()

# correct_output_dir_name = get_absolute_path('../correct_output')
# passed = True

# sig_digits = 5

# test_files = ["perm_x", "perm_y", "perm_z"]
# for test_file in test_files:
#     filename = f"/{$run_name}.out.{test_file}.pfb"
#     if not pf_test_file(new_output_dir_name + filename, correct_output_dir_name + filename, f"Max difference in {test_file}", sig_digits):
#         passed = False

# for i in range(6):
#     timestep = str(i).rjust(5, '0')
#     filename = f"/{run_name}.out.press.{timestep}.pfb"
#     if not pf_test_file(new_output_dir_name + filename, correct_output_dir_name + filename, f"Max difference in Pressure for timestep {timestep}"):
#         passed = False
#     filename = f"/{run_name}.out.satur.{timestep}.pfb"
#     if not pf_test_file(new_output_dir_name + filename, correct_output_dir_name + filename, f"Max difference in Saturation for timestep {timestep}"):
#         passed = False

# if passed:
#     print(f"{run_name} : PASSED")
# else:
#     print(f"{run_name} : FAILED")
#     sys.exit(1)
