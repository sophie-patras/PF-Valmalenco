
#Repository for hydro-climatic model of Valmalenco basin, using 'parflow-clm (v3.6.0)' 

Reference unit (SI) : [L] = m; [T] = h; [M] = kg
CRS. EPSG:32632

Repository architecture:
	- data/
	  Data inputs and outputs, distributed and local time series 
		- rawdata/ for the original version of files downloaded from websites
		- prepareddata/ for the elaborated data to fit the computational domain, grid and timing (format: .tif, .asc, ...). EPSG:32632.
		- pfdata/ for the formatted data to parflow binary and solid files (format: .pfb, .pfsol, .c.pfb, .pfsb). Files names are formated in <variable>.c<cellsize>.t<time>.<extension> where the time refers to real date : YYYY-MM-DD(-HH).
	- scripts/
	  Codes to be compiled in python3, combining 'parflow','gdal', and all other packages.
		- preprocess/ for GIS data elaboration (DataRaw>DataElab + DataElab>DataPF)
		- runprocess/ for parflow simulation. main is 'pfsimulation_runname.py', which calls .yaml files for parametrisation
			- tmp/
			  Temporary dir for input and output of 'parflow' computation (to speed up comunication). To be cleaned after each computation
		- posprocess/ for results elaboration-post processing (DataPF>DataRes)
	- outputs/
