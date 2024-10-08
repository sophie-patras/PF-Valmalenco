
Valmalenco

Repository for hydro-climatic model of Valmalenco basin, using 'parflow-clm (v3.6.0)' 

Reference unit (SI) : [L] = m; [T] = s; [M] = kg
CRS. EPSG:32632

Repository architecture:
	- Data/
	  Data inputs and outputs, distributed and local time series 
		- DataRaw/ for the original version of files downloaded from websites
		- DataElab/ for the elaborated data to fit the computational domain, grid and timing (format: .tif, .asc, ...). EPSG:32632.
		- DataPF/ for the formatted data to parflow binary and solid files (format: .pfb, .pfsol, .c.pfb, .pfsb). Files names are formated in <variable>.c<cellsize>.t<time>.<extension> where the time refers to real date : YYYY-MM-DD(-HH).
		- DataRes/ for the elaborated results and control points
	- Codes/
	  Codes to be compiled in python3, combining 'parflow','gdal', and all other packages.
		- PreProcess/ for GIS data elaboration (DataRaw>DataElab + DataElab>DataPF)
		- RunProcess/ for parflow simulation. main is 'pfsimulation_runname.py', which calls .yaml files for parametrisation
			- Tmp/
			  Temporary dir for input and output of 'parflow' computation (to speed up comunication). To be cleaned after each computation
		- PosProcess/ for results elaboration-post processing (DataPF>DataRes)
