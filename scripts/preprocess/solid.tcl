set tcl_precision 17

lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
pfset     FileVersion    4

#-----------------------------------------------------------------------------
# Create solid
#-----------------------------------------------------------------------------
set Mask [pfload -sa "../../Data/DataElab/maskb.c500.v4.sa"]
pfsetgrid {48 49 1} {0.0 0.0 0.0} {500.0 500.0 1.0} $Mask

set DMsk [pfload -sa "../../Data/DataElab/maskd.c500.v2.sa"]
pfsetgrid {48 49 1} {0.0 0.0 0.0} {500.0 500.0 1.0} $DMsk

# set DEM [pfload -sa "../../../Data/SA/burned_dem.sa"]
set DEM [pfload -sa "../../Data/DataElab/hydroDEM.c500.v2.sa"]
pfsetgrid {48 49 1} {0.0 0.0 0.0} {500.0 500.0 1.0} $DEM

set Top [pfcelldiff $DEM $DEM $DMsk]

set Top [pfcellsumconst $Top 1 $DMsk]
set Bottom [pfcellsumconst $Top -1 $DMsk]

pfpatchysolid -top $Top -bot $Bottom -msk $Mask -pfsol "../../Data/DataPF/solid.c500.vtcl3.pfsol" -vtk "../../Data/DataElab/solid.c500.vtcl3.vtk"
