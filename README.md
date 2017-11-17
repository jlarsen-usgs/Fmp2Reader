# Fmp2Reader
A tool to read in Modflow OWHM2 farm process files and create shapefiles of input data. 

This tool is installed as a python library using pip. Navigate to the trunk directory and use:

```
pip install -e .
```

to install fmp2reader as a python library

fmp2reader can then be imported to a python script using

```
from fmp2reader import Fmp
from fmp2reader import create_shapefile_from_array
from fmp2reader import create_shapefile_from_transient_array
```

Further useage can be found in the docstrings of methods and classes.