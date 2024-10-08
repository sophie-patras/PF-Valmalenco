# -*- coding: utf-8 -*-

"""
Python Parflow simulation, using yaml for input keys

Created on Mon Aug 5 11:37:26 2024
@author: Sophie Patras
"""

import os, sys, argparse
import numpy as np

from parflow import Run
#from parflow.tools.io import write_pfb, read_pfb, load_patch_matrix_from_asc_file #read_clm
#from parflow.tools.fs import cp, mkdir, chdir, get_absolute_path, rm
#from parflow.tools.settings import set_working_directory
#from parflow.tools.compare import pf_test_file
#from parflow.tools.builders import SolidFileBuilder


#______________________________________________________________________________________
parser = argparse.ArgumentParser()
parser.add_argument('run_name', type = str, default='test')
parser.add_argument('cellsize', type = float, default=500.0)
parser.add_argument('-parflow-directory', '--parflow-directory')
parser.add_argument('-working-directory', '--working-directory', default="~/Valmalenco/Codes/RunProcess/Tmp")
parser.add_argument('--show-line-error', action='store_true')
parser.add_argument('--validation-verbose', action='store_true')
# parser.add_argument('--write-yaml',action='store_true')
parser.add_argument('-p', '--p', type=int, default=2)
parser.add_argument('-q', '--q', type=int, default=2)
parser.add_argument('-r', '--r', type=int, default=1)
args = parser.parse_args()
# call argument : es Topology.P = args.p

runname = Run(args.run_name, __file__)

#_______________________________________________________________________________________
runname.pfset(yaml_file='PLT_35.yaml')

runname.dist("slopeX.c500.v4.pfb")
runname.dist("slopeY.c500.v4.pfb")

runname.run()

#undist (>> rm *.dist) $run_name, slopes.pfb, ip_solid.pfb
