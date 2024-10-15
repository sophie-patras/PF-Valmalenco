# -*- coding: utf-8 -*-

"""
Python Parflow simulation, using yaml for input keys

Created on Mon Aug 5 11:37:26 2024
@author: Sophie Patras
"""

import os, sys, re, argparse
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
parser.add_argument('cellsize', type = float, default=250.0)
parser.add_argument('-parflow-directory', '--parflow-directory')
parser.add_argument('-working-directory', '--working-directory', default="~/Valmalenco/Codes/RunProcess/Tmp")
parser.add_argument('--show-line-error', action='store_true')
parser.add_argument('--validation-verbose', action='store_true')
# parser.add_argument('--write-yaml',action='store_true')
parser.add_argument('-p', '--p', type=int, default=2)
parser.add_argument('-q', '--q', type=int, default=2)
parser.add_argument('-r', '--r', type=int, default=1)
args = parser.parse_args()

runname = Run(args.run_name, __file__)

#_______________________________________________________________________________________
#yamlfn = args.run_name + '_c' + str(int(args.cellsize)) + '.yaml'
#print(yamlfn)
runname.pfset(yaml_file='InputKeys_c250.yaml')

# dist
slopexfn = 'slopeX.c'+str(int(args.cellsize))+'.v1.pfb'
slopeyfn = 'slopeY.c'+str(int(args.cellsize))+'.v1.pfb'
#print(slopexfn)
runname.dist(slopexfn)
runname.dist(slopeyfn)
runname.dist('KS.c250.v2-331.pfb')
runname.dist('MRC_Tr.c250.v1-331.pfb')
runname.dist('MRC_Ts.c250.v1-331.pfb')
runname.dist('MRC_a.c250.v1-331.pfb')
runname.dist('MRC_n.c250.v1-331.pfb')

# clm
#variables = ['DSWR', 'DLWR', 'Press', 'SPFH', 'Temp', 'APCP', 'UGRD', 'VGRD']
#path_clm = '.'
#for variable in variables:
#    p = re.compile('clm.v1-331.'+ variable + '.000[0-9][0-9][0-9]_to_000[0-9][0-9][0-9].pfb')
#    for filename in sorted(os.listdir(path_clm)):
#        if p.match(filename):
#            fn_path = os.path.join(path_clm, filename)
#            runname.dist(fn_path)

# run
runname.run()

#end
