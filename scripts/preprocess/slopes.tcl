set tcl_precision 17

lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
pfset     FileVersion    4

#-----------------------------------------------------------------------------
# Process the slopes
#-----------------------------------------------------------------------------
set slopes_in_x [pfload -sa "/home/patras/Valmalenco/Data/DataElab/slopeX.c250.v1.sa"]
pfsetgrid {98 102 1} {0.0 0.0 0.0} {250.0 250.0 1.0} $slopes_in_x
pfsave $slopes_in_x -pfb "/home/patras/Valmalenco/Data/DataPF/slopeX.c250.v1.pfb"

set slopes_in_y [pfload -sa "/home/patras/Valmalenco/Data/DataElab/slopeY.c250.v1.sa"]
pfsetgrid {98 102 1} {0.0 0.0 0.0} {250.0 250.0 1.0} $slopes_in_y
pfsave $slopes_in_y -pfb "/home/patras/Valmalenco/Data/DataPF/slopeY.c250.v1.pfb"
