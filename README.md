# Fmp2Reader
A tool to read in Modflow OWHM2 farm process files and create shapefiles of input data. 

This tool is installed as a python library using pip. Navigate to the trunk directory and use:

```
pip install -e .
```

to install fmp2reader as a python library

fmp2reader can then be imported to a python script using

```python
from fmp2reader import Fmp
from fmp2reader import create_shapefile_from_array
from fmp2reader import create_shapefile_from_transient_array


# Read in an fmp file and set to Fmp as an overridden dictionary
fmp = Fmp('lucerne.fmp', fmp_ws=r"C:/mydirectory", nrow=50, ncol=64,
          crop_lut={1:'alfalfa', 2:'apples'}, cell_area=4000000)

# Build fmp data into static and transient arrays
fmp.to_arrays()

# get dictionarys of farm process inputs
crops = fmp.crop_arrays
areas = fmp.crop_areas

# grab the farm arrays
farms = fmp.transient_arrays['farms']

# create shapefiles... xgrid and ygrid are numpy based model grid 
# vertex arrays in the shape nrow + 1, ncol + 1

# note: all input arrays must be 4d arrays for transient shapefiles
# (nper, nlay, nrow, ncol) nlay = 1 for farm processes.

create_shapefile_from_transient_array(r'C:\mygis\crops.shp', fmp.crop_arrays, 
				  nper=270, nlay=1, xgrid=xgrid, 
				  ygrid=ygrid, no_data=0)
create_shapefile_from_transient_array(r'C:\mygis\farms.shp', fmp.transient_arrays, 
				  nper=270, nlay=1, xgrid=xgrid, 
				  ygrid=ygrid, no_data=0)

```

Further useage can be found in the docstrings of methods and classes.