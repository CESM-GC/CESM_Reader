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
if '/run' in folder:
    folder = os.path.join(rootFolder, folder)
else:
    folder = os.path.join(rootFolder, folder, 'run')
fileReader = CESM_Reader(folder, debug=False);

include = ['O3', 'NO*', 'CFC11']
fileReader.register(include=include, debug=False)
#print(fileReader.include) # This should be equal to ['O3', 'CFC11', 'NO', 'NO2', 'NO3']

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