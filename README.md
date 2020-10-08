CESM_Reader is a Python processing package to load, process and plot output data from CESM. It was initially intended to process data from the Community Atmosphere Model (CAM).

CESM_Reader is able to perform different types of spatial averaging:
- Zonal average
- Layer average across multiple layers
- Column average
- Global average

It is also able to either compute temporal mean or to load temporal variations of the desired fields.

Below is an example that I extracted from a Jupyter notebook. 

```python
import os
import numpy as np
from datetime import datetime

from CESM_Reader.fileReader import CESM_Reader

rootFolder = '/glade/scratch/fritzt'
folder     = 'CESM-GC_SpeciesConc'
devFolder  = '/path/to/dev/folder'
if '/run' in folder:
    folder    = os.path.join(rootFolder, folder)
    devFolder = os.path.join(rootFolder, devFolder)
else:
    folder    = os.path.join(rootFolder, folder, 'run')
    devFolder = os.path.join(rootFolder, devFolder, 'run')

# Load data from a run
fileReader = CESM_Reader(folder, debug=False);

# Load differences between two CESM runs
#fileReader = CESM_Reader(folder, devFolder=devFolder, debug=False);

# As long as the following species are saved during the CESM run, then they will be loaded
include = ['O3', 'NO*', 'pFe', 'DepFlux_NO*', 'DepVel_O3', 'WDRATE_SO2', 'Jval_O3*']
fileReader.register(include=include, debug=False)
#print(fileReader.include)

## Zonal averaging
fileReader.load(spatialAveraging='Zonal', timeAveraging=True, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

## Layer slicing
#fileReader.load(spatialAveraging='Layer', iLayer=-1, timeAveraging=True, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

## Temporal variations of globally-averaged mixing ratio
#fileReader.load(spatialAveraging='All', timeAveraging=False, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

## Altitude variation - 1D
#fileReader.load(spatialAveraging='V Altitude', timeAveraging=True, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

## Latitude variation - 1D
#fileReader.load(spatialAveraging='V Latitude', timeAveraging=True, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

## Longitude variation - 1D
#fileReader.load(spatialAveraging='V Longitude', timeAveraging=True, minDate=datetime(2005,1,1), maxDate=datetime(2005,1,3), debug=False)

im, cb = fileReader.plot(species=include, plotUnit='ppb', clim=np.array([0,np.NaN]), debug=False)
```
