#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset
from glob import glob
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from datetime import datetime, timedelta
import cartopy.io.shapereader as shpreader
import logging

SUCCESS     = 0
WRONG_UNIT  = 1
WRONG_SHAPE = 2

MIN_DATE = datetime(1,1,1)
MAX_DATE = datetime(9999,12,31)
reader = shpreader.Reader('/glade/u/home/fritzt/.local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_coastline.shp')

class fileHandler:
    def __init__(self, rootFolder, debug=False):

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        self.rootFolder    = rootFolder
        self.fileList      = {}
        self.logList       = {}
        self.wildcard      = '*'
        self.comps         = ['cam', 'clm']
        self.fileExtension = '.nc'
        self.logExtension  = '.log.*'
        _logs  = ['atm', 'cesm', 'cpl', 'ice', 'lnd', 'ocn', 'rof']
        _types = ['r', 'i', 'h']
        _nTape = {}

        CASEROOT = os.path.join(self.rootFolder, 'CASEROOT')
        if os.path.exists(CASEROOT):
            f = open(CASEROOT, 'r')
            self.CASEROOT = f.read()
            self.CASEROOT.replace('\n','')
            f.close()
            print('CASEROOT is {:s}'.format(self.CASEROOT))

        for _comp in self.comps:
            _nTape[_comp]        = 0
            self.fileList[_comp] = {}
            for _type in _types:
                if _type == 'h':
                    isHist      = True
                    _tapeNumber = 0
                    _tape       = str(_tapeNumber)
                    _maxTape    = 10
                else:
                    isHist      = False
                    _tapeNumber = 0
                    _tape       = ''
                    _maxTape    = 1
                while (_tapeNumber < _maxTape):
                    _fileName = self.wildcard + _comp + '.' + _type + _tape + self.wildcard + self.fileExtension
                    keyName = _type + _tape
                    tmp = sorted(glob(os.path.join(self.rootFolder, _fileName)))
                    if len(tmp) > 0:
                        self.fileList[_comp][keyName] = tmp.copy()
                        logging.debug('Found {:3d} files corresponding to component {:4s} and type {:3s}'.format(
                            len(self.fileList[_comp][keyName]), _comp, _type + _tape))
                        _nTape[_comp] += 1
                    _tapeNumber += 1
                    _tape  = str(_tapeNumber)
            logging.debug('Found {:2d} tapes for component {:4s}'.format(_nTape[_comp], _comp))

        _jobIDs   = []
        _runDates = []
        _delim1   = '.log.'
        _delim2   = '.'
        for _log in _logs:
            _fileName = _log + self.logExtension
            self.logList[_log] = sorted(glob(os.path.join(self.rootFolder, _fileName)))
            for _logName in self.logList[_log]:
                _start = _logName.index(_delim1) + len(_delim1)
                try:
                    _jobID = _logName[_start:_logName.index(_delim2,_start)]
                    if _jobID not in _jobIDs:
                        _jobIDs.append(_jobID)
                        if _logName[-3:] == '.gz':
                            _runDates.append(_logName[-16:-3])
                        else:
                            _runDates.append(_logName[-13:])
                except:
                    _jobID = 999999
            if len(self.logList[_log]) > 0:
                logging.debug('Found {:2d} log files for component {:4s}'.format(len(self.logList[_log]), _log))
        logging.debug('.. corresponding to job IDs:')
        for iJob, _jobID in enumerate(_jobIDs):
            logging.debug(' - Job {:8s} on {:13s}'.format(_jobID, _runDates[iJob]))

        # Open most recent ATM log file
        if self.logList['atm']:
            _mRecentFile = max(self.logList['atm'], key=os.path.getctime)
            if _mRecentFile[-3:] != '.gz':
                f = open(_mRecentFile, 'r')
                foundHEMCO = False
                for line in f:
                    if 'HEMCO_Config.rc' in line:
                        foundHEMCO = True
                        line.replace('\n', '')
                        print('Using HEMCO from {:s}'.format(line))
                        break;
                f.close()
                if foundHEMCO == False:
                    print('No mention of HEMCO_Config.rc was found in {:s}'.format(_mRecentFile))

class CESM_Reader:
    def __init__(self, rootFolder, devFolder=None, surfPres=1013.25, debug=False):

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        self.rootFolder = rootFolder
        self.devFolder  = devFolder

        self.diffRun = False
        if devFolder is not None and os.path.exists(devFolder):
            self.diffRun = True

        self.fileInst = fileHandler(rootFolder, debug=debug);
        # Find any CAM output file to load grid characteristics
        _comp = 'cam'
        fId   = None
        for _tape in self.fileInst.fileList[_comp].keys():
            if 'h' in _tape and len(self.fileInst.fileList[_comp][_tape]) > 0:
                fId = Dataset(self.fileInst.fileList[_comp][_tape][0], 'r')
                break;
        if fId == None:
            raise ValueError('Could not figure out grid characteristics for component CAM')
        self.pEdge    = fId['hyai'][:]*fId['P0'][:]/1.0E+02+fId['hybi'][:]*1013.25 #hPa
        self.pMid     = fId['lev'][:] #(self.pEdge[1:]+self.pEdge[:-1])/2.0 # hPa
        self.timeMid  = np.array(())
        self.lat      = fId['lat'][:].data
        self.lon      = fId['lon'][:].data
        self.nLat     = len(self.lat)
        self.nLon     = len(self.lon)
        self.nLev     = len(self.pMid)
        self.nTime    = {}
        self.longName = {}
        self.loaded   = False

        self.isCentered = False
#         if self.lon[0] >= 0:
#             self.isCentered = True
#             self.lon -= 180.0 # Recenter on 0
        self.lonEdge = np.append(self.lon[:], self.lon[-1] + self.lon[-1] - self.lon[-2])
        self.latEdge = np.append(self.lat[:], self.lat[-1] + self.lat[-1] - self.lat[-2])

        print('\nGrid characteristics: {:2d}x{:2d}'.format(len(self.lon), len(self.lat)))
        print('Lowest pressure is  : {:3.2f} hPa'.format(np.min(self.pEdge)))

        if self.diffRun:
            self.fileInst_dev = fileHandler(devFolder, debug=debug)

        self.possUnit = {}
        # Conversion factor to go from current species to mol/mol air
        self.possUnit['mol/mol']  = 1.0E+00
        self.possUnit['-']        = 1.0E+00
        self.possUnit['ppth']     = 1.0E-03
        self.possUnit['ppthv']    = self.possUnit['ppth']
        self.possUnit['ppm']      = 1.0E-06
        self.possUnit['ppmv']     = self.possUnit['ppm']
        self.possUnit['umol/mol'] = self.possUnit['ppm']
        self.possUnit['ppb']      = 1.0E-09
        self.possUnit['ppbv']     = self.possUnit['ppb']
        self.possUnit['nmol/mol'] = self.possUnit['ppb']
        self.possUnit['ppt']      = 1.0E-12
        self.possUnit['pptv']     = self.possUnit['ppt']
        self.possUnit['pmol/mol'] = self.possUnit['ppt']

    def register(self, include=[], exclude=[], excludeTape={}, debug=False):

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        self.include = include
        self.exclude = exclude

        if len(self.include) > 0:
            self.data           = {}
            self.unit           = {}
            self.found          = {}
            self.locat          = {}
            self.compList       = []
            self.tapeList       = {}
            self.includePerTape = {}

            # List tapes to be excluded
            if isinstance(excludeTape, int):
                tmp = excludeTape
                excludeTape = {}
                for _comp in self.fileInst.comps:
                    excludeTape[_comp] = np.array([tmp])
            else:
                for _comp in self.fileInst.comps:
                    if _comp in excludeTape.keys():
                        if isinstance(excludeTape[_comp], int):
                            excludeTape[_comp] = np.array([excludeTape[_comp]])
                    else:
                        excludeTape[_comp] = np.array([])

            # For each species find the component and tape associated to it
            # First, load first file of each component/tape combination
            fId = {}
            for _comp in self.fileInst.comps:
                fId[_comp] = {}
                for _tape in self.fileInst.fileList[_comp]:
                    if (('h' in _tape) and (int(_tape[-1]) not in excludeTape[_comp]) and (len(self.fileInst.fileList[_comp][_tape]) > 0)):
                        fId[_comp][_tape] = Dataset(self.fileInst.fileList[_comp][_tape][0], 'r')

            def IsInExclude(spec):
                found = False
                for exclSpec in self.exclude:
                    if self.fileInst.wildcard in exclSpec:
                        exclSpec2 = exclSpec.replace(self.fileInst.wildcard, '')
                        found = (exclSpec2 in spec)
                    else:
                        found = (exclSpec == spec)
                    if found == True:
                        logging.debug('Species {:s} is excluded because it matches with {:s}'.format(spec, exclSpec))
                        return found
                return found

            # Then find location of each species in these files
            tmp = self.include.copy()
            for spec in tmp:
                self.found[spec] = False
                self.locat[spec] = []
                if self.fileInst.wildcard not in spec:
                    for _comp in self.fileInst.comps:
                        if _comp not in self.includePerTape.keys():
                            self.includePerTape[_comp] = {}
                        for _tape in fId[_comp].keys():
                            if _tape not in self.includePerTape[_comp].keys():
                                self.includePerTape[_comp][_tape] = []
                            if ((not self.found[spec]) and (spec in fId[_comp][_tape].variables.keys()) and (not IsInExclude(spec))):
                                self.found[spec] = True
                                self.locat[spec] = [_comp, _tape]
                                if spec not in self.includePerTape[_comp][_tape]:
                                    self.includePerTape[_comp][_tape] += [spec]
                                logging.debug('{:20s} was found in component {:s} and type {:s}'.format(spec, _comp, _tape))
                elif self.fileInst.wildcard in spec:
                    # Remove wilcard from string
                    specSplit = spec.split(self.fileInst.wildcard)
                    len1 = len(specSplit[0])
                    len2 = len(specSplit[1])
                    # Find all species matching that are at least 2D
                    for _comp in self.fileInst.comps:
                        if _comp not in self.includePerTape.keys():
                            self.includePerTape[_comp] = {}
                        for _tape in fId[_comp].keys():
                            if _tape not in self.includePerTape[_comp].keys():
                                self.includePerTape[_comp][_tape] = []
                            for field in fId[_comp][_tape].variables.keys():
                                if ((len1 > 0) and (len2 > 0)):
                                    willInclude = ((specSplit[0] == field[:len1]) and (specSplit[1] == field[-len2:]))
                                elif ((len1 > 0) and (len2 == 0)):
                                    willInclude = ((specSplit[0] == field[:len1]))
                                elif ((len1 == 0) and (len2 > 0)):
                                    willInclude = ((specSplit[1] == field[-len2:]))
                                elif ((len1 == 0) and (len2 == 0)):
                                    willInclude = True
                                willInclude = willInclude and (len(fId[_comp][_tape][field].dimensions) > 1)
                                if willInclude:
                                    if (field not in self.include) and (not IsInExclude(field)):
                                        self.include.insert(len(self.include), field)
                                        self.found[field] = True
                                        self.locat[field]  = [_comp, _tape]
                                        if field not in self.includePerTape[_comp][_tape]:
                                            self.includePerTape[_comp][_tape] += [field]
                                        logging.debug('{:20s} has been added. It is found in component {:s} and type {:s}'.format(field, _comp, _tape))

            # Remove species that have not been found
            # If a species contains a wildcard, then it will be removed,
            # as the subsequent species should have been added above
            tmp = self.include.copy()
            for spec in tmp:
                if not self.found[spec]:
                    if self.fileInst.wildcard not in spec:
                        logging.warning('{:s} was not found. It will be removed!'.format(spec))
                    self.include.remove(spec)

            for spec in self.include:
                _comp = self.locat[spec][0]
                _tape = self.locat[spec][1]
                if (_comp not in self.compList):
                    self.compList.insert(len(self.compList), _comp)
                    self.tapeList[_comp] = []
                if (_tape not in self.tapeList[_comp]):
                    self.tapeList[_comp].insert(len(self.tapeList[_comp]), _tape)

            # Now close files
            for _comp in fId.keys():
                for _tape in fId[_comp].keys():
                    fId[_comp][_tape].close()
            del fId

        else:
            raise ValueError('No species were selected')

    def load(self, spatialAveraging='', timeAveraging=True,
             maxFile=-1, minDate=MIN_DATE, maxDate=MAX_DATE,
             iLayer=np.array([-1]), doPlot=False, plotUnit=None,
             debug=False):

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        self.spatialAveraging = spatialAveraging
        self.timeAveraging    = timeAveraging

        timeDim = 0
        levDim  = 1
        latDim  = 2
        lonDim  = 3
        doLayer = False

        _avMethods = {}
        _avMethods['zonal']       = 'Data is averaged along the longitudinal axis'
        _avMethods['layer']       = 'Data is extracted at a given layer, with the iLayer argument. (TOA -> iLayer=0, Surface -> iLayer=-1)'
        _avMethods['column']      = 'Data is summed along each column'
        _avMethods['all']         = 'Data is globally averaged (lev, lat, lon)'
        _avMethods['v altitude']  = 'Data is averaged across longitude and latitude'
        _avMethods['v latitude']  = 'Data is averaged across longitude and altitude'
        _avMethods['v longitude'] = 'Data is averaged across latitude and altitude'

        self.zonalSize   = (self.nLev,self.nLat)
        self.layerSize   = (self.nLat,self.nLon)
        self.altSize     = (self.nLev,)
        self.latSize     = (self.nLat,)
        self.lonSize     = (self.nLon,)
        _size_2D         = self.layerSize
        if self.spatialAveraging.lower() == 'zonal':
            logging.debug('Zonal averaging will be performed!')
            _avAxis      = (lonDim,)
            _size        = self.zonalSize
            _scaleFactor = 1.0E+00
        elif self.spatialAveraging.lower() == 'layer':
            doLayer      = True
            _avAxis      = (levDim,)
            _size        = self.layerSize
            _scaleFactor = 1.0E+00
            if isinstance(iLayer, int):
                logging.debug('Extracting single layer!')
                self.iLayer = np.array([iLayer])
            elif isinstance(iLayer, np.ndarray):
                if np.shape(iLayer)[0] > 1:
                    logging.debug('Extracting multiple layers!')
                elif np.shape(iLayer)[0] == 1:
                    logging.debug('Extracting single layer!')
                self.iLayer = iLayer
            else:
                logging.error('Argument iLayer was not understood!')
        elif self.spatialAveraging.lower() == 'column':
            logging.debug('Computing column totals!')
            _avAxis      = (levDim,)
            _size        = self.layerSize
            _scaleFactor = self.nLev
        elif self.spatialAveraging.lower() == 'all':
            logging.debug('Global averaging will be performed!')
            _avAxis      = (levDim,latDim,lonDim,)
            _size        = ()
            _scaleFactor = 1.0E+00
        elif self.spatialAveraging.lower() == 'v altitude':
            logging.debug('Longitude/Latitude averaging!')
            _avAxis      = (latDim, lonDim,)
            _size        = self.altSize
            _scaleFactor = 1.0E+00
        elif self.spatialAveraging.lower() == 'v latitude':
            logging.debug('Longitude/Altitude averaging!')
            _avAxis      = (levDim, lonDim,)
            _size        = self.latSize
            _scaleFactor = 1.0E+00
        elif self.spatialAveraging.lower() == 'v longitude':
            logging.debug('Latitude/Altitude averaging!')
            _avAxis      = (levDim, latDim,)
            _size        = self.lonSize
            _scaleFactor = 1.0E+00
        else:
            # This should correspond to the case where:
            # self.spatialAveraging.lower() not in _avMethods.keys():
            _avAxis      = ()
            _size        = (self.nLev,self.nLat,self.nLon)
            _scaleFactor = 1.0E+00
            logging.warning('No or bad averaging method chosen!')
            logging.warning('Other options include: ')
            for _avMethod in _avMethods.keys():
                logging.warning('- {:20s}: {:s}'.format(_avMethod, _avMethods[_avMethod]))
            logging.warning('For each output field, the data\'s dimension will be in ((time), lev, lat, lon)!')

        # For arrays that are only (time,lat,lon), create a new _avAxis from which levDim is removed
        # _avAxis_2D is designed for (time, lat, lon)
        _avAxis_2D = ()
        _avAxis_3D = ()
        _size_3D   = ()
        _size_2D   = ()
        for _dim in _avAxis:
            _tmp = _dim
            if _dim > levDim:
                _tmp -= 1
                _avAxis_3D += (_tmp,)
                _avAxis_2D += (_tmp,)
            elif _dim < levDim:
                _avAxis_3D += (_tmp,)
        for _tmp in _size:
            if _tmp != self.nLev:
                _size_3D += (_tmp,)
                _size_2D += (_tmp,)

        _avAxis_2D_noT = _avAxis_2D
        _axAxis_3D_noT = _avAxis_3D
        if self.timeAveraging:
            logging.debug('Time averaging will be performed!')
            _avAxis    = (timeDim,) + _avAxis
            _avAxis_2D = (timeDim,) + _avAxis_2D
        else:
            logging.debug('Time variations will be considered!')
            #logging.debug('A time slice will be extracted using the iTime argument')
            # This required to know how many time slices are in each file
            #_size        = (self.nTime,) + _size

        if maxFile >= 0:
            logging.debug('Only loading {:d} files'.format(maxFile))
        if minDate >= maxDate:
            raise ValueError('Minimum date is greater than maximum date!')
        if minDate != MIN_DATE:
            logging.debug('Only loading from  {:s}'.format(minDate.strftime("%Y-%m-%d")))
        if maxDate != MAX_DATE:
            logging.debug('Only loading up to {:s}'.format(maxDate.strftime("%Y-%m-%d")))

        firstFile = {}
        nFiles    = {}
        for spec in self.include:
            firstFile[spec] = True
            nFiles[spec]    = 0

        # Read
        for _comp in self.compList:
            for _tape in self.tapeList[_comp]:
                filesRead = 0
                for iFile, file in enumerate(self.fileInst.fileList[_comp][_tape]):

                    if (maxFile > 0) and (iFile >= maxFile):
                        break

                    fileName = self.fileInst.fileList[_comp][_tape][iFile]
                    if self.diffRun:
                        fileDev  = self.fileInst_dev.fileList[_comp][_tape][iFile]

                    YYYY = int(os.path.basename(fileName)[-19:-15])
                    MM   = int(os.path.basename(fileName)[-14:-12])
                    DD   = int(os.path.basename(fileName)[-11:-9])
                    currDate = datetime(YYYY, MM, DD)

                    # Rather than extracting the date from the file name, we should
                    # extract it from the netCDF variables.
                    _discard = False
                    if currDate >= maxDate:
                        break
                    if currDate < minDate:
                        _discard = True

                    if not _discard:
                        try:
                            fId = Dataset(fileName, 'r')
                            if self.diffRun:
                                fIdDev = Dataset(fileDev, 'r')
                            isAFile = True
                        except:
                            isAFile = False
                            logging.warning('File {:s} was not found'.format(fileName))

                        if isAFile:

                            logging.debug('Dealing with file {:03d} for {:s}'.format(
                                    iFile, currDate.strftime("%Y-%m-%d")))
                            filesRead += 1

                            for iSpec, spec in enumerate(self.includePerTape[_comp][_tape]):
                                if firstFile[spec]:
                                    # Get long name for this species
                                    self.longName[spec] = fId[spec].long_name
                                    # Number of time samples in this file
                                    _localSample = np.shape(fId['time'])[0]
                                    # Number of time samples for all files
                                    self.nTime[spec] = len(self.fileInst.fileList[_comp][_tape]) * _localSample
                                    # Adding an empty tuple allows to copy by memory and not by address!
                                    _tsize_noT = _size + ()
                                    _tsize     = _size + ()
                                    _tsize_2D  = _size_2D + ()
                                    _tsize_3D  = _size_3D + ()
                                    if self.timeAveraging == False:
                                        _tsize    = (_localSample,) + _tsize
                                        _tsize_2D = (_localSample,) + _tsize_2D
                                        _tsize_3D = (_localSample,) + _tsize_3D
                                        if iSpec == 0:
                                            _tmpArray_time = np.array(currDate, dtype=np.datetime64) + np.arange(_localSample)
                                            self.timeMid = _tmpArray_time.copy()
                                            _date    = fId['date'][:].data
                                            _timeSec = fId['datesec'][:].data
                                            self.timeMid[(filesRead-1)*_localSample:filesRead*_localSample] = np.array([datetime.strptime(str(_currDate), '%Y%m%d') +                                                     timedelta(seconds=np.int(_currSec)) for (_currDate, _currSec) in zip(_date, _timeSec)])
                                        _tsize_noT = _size + ()
                                        _tmpArray_noT = np.zeros(_tsize_noT)
                                    _tmpArray       = np.zeros(_tsize)
                                    _tmpArray_2D    = np.zeros(_tsize_2D)
                                    _tmpArray_3D    = np.zeros(_tsize_3D)

                                    self.unit[spec]  = fId[spec].units
                                else:
                                    if ((self.timeAveraging == False) and (iSpec == 0)):
                                        _tmpArray_time = np.array(currDate, dtype=np.datetime64) + np.arange(_localSample)
                                        self.timeMid = np.concatenate((self.timeMid, _tmpArray_time.copy()), axis=0)
                                        _date    = fId['date'][:].data
                                        _timeSec = fId['datesec'][:].data
                                        self.timeMid[(filesRead-1)*_localSample:filesRead*_localSample] = np.array([datetime.strptime(str(_currDate), '%Y%m%d') + \
                                                timedelta(seconds=np.int(_currSec)) for (_currDate, _currSec) in zip(_date, _timeSec)])

                                if (len(np.shape(fId[spec])) == 4):
                                    # This is made for (time, lev, lat, lon) fields
                                    if firstFile[spec]:
                                        self.data[spec] = _tmpArray.copy()
                                    elif not timeAveraging:
                                        self.data[spec] = np.concatenate((self.data[spec], _tmpArray), axis=0)

                                    if timeAveraging:
                                        if doLayer:
                                            self.data[spec] += np.mean(fId[spec][:,self.iLayer,:,:], axis=_avAxis)
                                            if self.diffRun:
                                                self.data[spec] -= np.mean(fIdDev[spec][:,self.iLayer,:,:], axis=_avAxis)
                                        else:
                                            self.data[spec] += np.mean(fId[spec], axis=_avAxis)
                                            if self.diffRun:
                                                self.data[spec] -= np.mean(fIdDev[spec], axis=_avAxis)
                                        nFiles[spec] += 1
                                    else:
                                        currTimeIndex = (filesRead-1) * _localSample
                                        if doLayer:
                                            self.data[spec][currTimeIndex:currTimeIndex+_localSample] = np.mean(fId[spec][:,self.iLayer,:,:], axis=_avAxis)
                                            if self.diffRun:
                                                self.data[spec][currTimeIndex:currTimeIndex+_localSample] -= np.mean(fIdDev[spec][:,self.iLayer,:,:], axis=_avAxis)
                                        else:
                                            self.data[spec][currTimeIndex:currTimeIndex+_localSample] = np.mean(fId[spec], axis=_avAxis)
                                            if self.diffRun:
                                                self.data[spec][currTimeIndex:currTimeIndex+_localSample] -= np.mean(fIdDev[spec], axis=_avAxis)
                                        nFiles[spec] = 1
                                elif (len(np.shape(fId[spec])) == 3):
                                    # This should make distinction between (time, lat, lon) and (lev, lat, lon)
                                    if 'time' in fId.dimensions.keys():
                                        # Then we have a temporally-evolving layer-like field in (time, lat, lon)
                                        if firstFile[spec]:
                                            self.data[spec] = _tmpArray_2D.copy()
                                        elif not timeAveraging:
                                            self.data[spec] = np.concatenate((self.data[spec], _tmpArray_2D), axis=0)
                                        if timeAveraging:
                                            self.data[spec] += np.mean(fId[spec], axis=_avAxis_2D)
                                            if self.diffRun:
                                                self.data[spec] -= np.mean(fIdDev[spec], axis=_avAxis_2D)
                                            nFiles[spec] += 1
                                        else:
                                            currTimeIndex = (filesRead-1) * _localSample
                                            self.data[spec][currTimeIndex:currTimeIndex+_localSample] = np.mean(fId[spec], axis=_avAxis_2D)
                                            if self.diffRun:
                                                self.data[spec][currTimeIndex:currTimeIndex+_localSample] -= np.mean(fIdDev[spec], axis=_avAxis_2D)
                                            nFiles[spec] = 1
                                    else:
                                        # Else, we have a constant (lev, lat, lon) field
                                        if firstFile[spec]:
                                            self.data[spec] = _tmpArray_noT.copy()
                                        elif not timeAveraging:
                                            self.data[spec] = np.concatenate((self.data[spec], _tmpArray_noT), axis=0)

                                        if timeAveraging:
                                            if doLayer:
                                                self.data[spec] += np.mean(fId[spec][self.iLayer,:,:], axis=_avAxis_2D_noT)
                                                if self.diffRun:
                                                    self.data[spec] -= np.mean(fIdDev[spec][self.iLayer,:,:], axis=_avAxis_2D_noT)
                                            else:
                                                self.data[spec] += np.mean(fId[spec], axis=_avAxis_3D_noT)
                                                if self.diffRun:
                                                    self.data[spec] -= np.mean(fIdDev[spec], axis=_avAxis_3D_noT)
                                            nFiles[spec] += 1
                                        else:
                                            currTimeIndex = filesRead-1
                                            if doLayer:
                                                self.data[spec][currTimeIndex] = np.mean(fId[spec][self.iLayer,:,:], axis=_avAxis_2D_noT)
                                                if self.diffRun:
                                                    self.data[spec][currTimeIndex] -= np.mean(fIdDev[spec][self.iLayer,:,:], axis=_avAxis_2D_noT)
                                            else:
                                                self.data[spec][currTimeIndex] = np.mean(fId[spec], axis=_avAxis_3D_noT)
                                                if self.diffRun:
                                                    self.data[spec][currTimeIndex] -= np.mean(fIdDev[spec], axis=_avAxis_3D_noT)
                                            nFiles[spec] = 1
                                else:
                                    logging.warning('{:s} will not be read from file! Found {:d} dimensions'.format(spec, len(np.shape(fId[spec]))))

                                if firstFile[spec]:
                                    firstFile[spec] = False

        for iSpec, spec in enumerate(self.include):
            if nFiles[spec] > 0:
                self.data[spec] *= _scaleFactor / nFiles[spec]
            else:
                raise ValueError('No files were found for {:s}'.format(spec))
            # If we performed a global and temporal average, then just print out
            if len(np.shape(self.data[spec])) == 0:
                if iSpec == 0:
                    print('Globally and temporal averaged mixing ratios')
                    print('----------------+--------------------------')
                if self.unit[spec] in self.possUnit.keys():
                    _val  = self.data[spec]
                    _unit = self.unit[spec]
                    if _val < 1.0E-02:
                        _val *= self.possUnit[_unit] / self.possUnit['ppth']
                        _unit = 'ppth'
                    if _val < 1.0E+00:
                        _val *= self.possUnit[_unit] / self.possUnit['ppmv']
                        _unit = 'ppmv'
                    if _val < 1.0E+00:
                        _val *= self.possUnit[_unit] / self.possUnit['ppbv']
                        _unit = 'ppbv'
                    if _val < 1.0E+00:
                        _val *= self.possUnit[_unit] / self.possUnit['pptv']
                        _unit = 'pptv'
                    print('{:16s}|{:6.2e} {:s}'.format(spec, _val, _unit))

        # Canons ready, captain!
        self.loaded = True

        if doPlot:
            self.plot(species=self.include, plotUnit=plotUnit, debug=debug)

    def save(self, fileName=None,
             debug=False):

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        if not self.loaded:
            logging.error('Data has not been loaded previously')

        if fileName is None:
            logging.error('fileName is None!')

        if not self.timeAveraging:
            logging.error('Function save not yet implemented for timeAveraging = False')

        # Create file
        fId = Dataset(fileName, 'w', format="NETCDF4")

        # Add dimensions
        fId.createDimension('Latitude', self.nLat)
        fId.createDimension('Longitude', self.nLon)
        fId.createDimension('Altitude', self.nLev)

        # Add variables
        vout = fId.createVariable('lat', np.dtype('float32').char, 'Latitude')
        vout[:] = self.lat
        vout = fId.createVariable('lon', np.dtype('float32').char, 'Longitude')
        vout[:] = self.lon
        vout = fId.createVariable('lev', np.dtype('float32').char, 'Altitude')
        vout[:] = self.pMid
        for key in self.data.keys():
            if ( np.shape(self.data[key]) == (self.nLat, self.nLon) ):
                # Layer
                shape = ('Latitude', 'Longitude')
                transpose = False
            elif ( np.shape(self.data[key]) == (self.nLev, self.nLat) ):
                # Zonal
                shape = ('Altitude', 'Latitude')
                transpose = False
            elif ( np.shape(self.data[key]) == (self.nLat) ):
                # Latitude
                shape = ('Latitude')
                transpose = False
            elif ( np.shape(self.data[key]) == (self.nLon) ):
                # Longitude
                shape = ('Longitude')
                transpose = False
            elif ( np.shape(self.data[key]) == (self.nLev) ):
                # Altitude
                shape = ('Altitude')
                transpose = False
            vout = fId.createVariable(key, np.dtype('float32').char, shape)
            if transpose:
               vout[:] = np.transpose(self.data[key])
            else:
               vout[:] = self.data[key]

        # Close file
        fId.close()

        print('Data was saved to {:s}'.format(fileName))

    def plot(self, data=None, species=None, dataUnit=None, plotUnit=None,
             cmap=None, clim=None, ylim=None, xlim=None,
             show_colorbar=True, cbTitle=None, logScale=False,
             labelFtSize=18, labelTickSize=18, isDiff=False,
             debug=False):

        imageDict = {}
        cbDict    = {}
        figDict   = {}

        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        if debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        if not self.loaded:
            logging.error('Data has not been loaded previously')

        if (not isinstance(data, np.ndarray)) and (species == None):
            logging.error('No data was passed to plot')

        spc2Plot = ['None']
        if (not isinstance(data, np.ndarray)):
            if isinstance(species, str):
                spc2Plot = [species]
            elif isinstance(species, list):
                spc2Plot = species
            for spec in spc2Plot:
                if spec not in self.include:
                    logging.error('Species {:s} is not part of loaded species'.format(spec))
                    return imageDict, cbDict, figDict

            if self.spatialAveraging.lower() not in ['zonal', 'layer', 'column', 'all', \
                    'v altitude', 'v latitude', 'v longitude']:
                logging.error('No plotting method exists for spatialAveraging = {:s}'.
                              format(self.spatialAveraging))
                return imageDict, cbDict, figDict

            if ((self.timeAveraging == False) and (self.spatialAveraging.lower() not in ['all', \
                    'v altitude', 'v latitude', 'v longitude'])):
                logging.error('No plotting method exists to plot non-globally averaged temporal variations')
                return imageDict, cbDict, figDict

        custUnit_isStr  = False
        custUnit_isList = False
        custUnit_isDict = False
        targetUnit      = 'ppbv'
        dataUnit_isStr  = False
        dataUnit_isList = False
        dataUnit_isDict = False

        if plotUnit is not None:
            if isinstance(plotUnit, str):
                custUnit_isStr = True
            elif (isinstance(plotUnit, list)) and (len(plotUnit) == len(spc2Plot)):
                custUnit_isList = True
            elif (isinstance(plotUnit, dict)):
                custUnit_isDict = True
        if dataUnit is not None:
            if isinstance(dataUnit, str):
                dataUnit_isStr = True
            elif (isinstance(dataUnit, list)) and (len(dataUnit) == len(spc2Plot)):
                dataUnit_isList = True
            elif (isinstance(dataUnit, dict)):
                dataUnit_isDict = True

        if logScale and clim is not None:
            if clim[0] <= 0:
                logging.error('Specified clim cannot be negative for logScale = True')

        for iSpec, spec in enumerate(spc2Plot):
            im       = None
            fig      = None
            cb       = None
            currUnit = None
            currSpec = None
            if spec != 'None':
                data     = self.data[spec]
                currUnit = self.unit[spec]
                currSpec = spec
            if custUnit_isStr:
                targetUnit = plotUnit
            elif custUnit_isList:
                targetUnit = plotUnit[iSpec]
            elif custUnit_isDict:
                try:
                    targetUnit = plotUnit[spec]
                except:
                    logging.warning('Could not find targetUnit for species {:s}'.format(spec))
                    targetUnit = 'ppbv'
            if currUnit is None:
                if dataUnit_isStr:
                    currUnit = dataUnit
                elif dataUnit_isList:
                    currUnit = dataUnit[iSpec]
                elif dataUnit_isDict:
                    try:
                        currUnit = dataUnit[spec]
                    except:
                        logging.warning('Could not find currUnit for species {:s}'.format(spec))
                        currUnit = 'ppbv'
            nDim = len(np.shape(data))
            RC = -1
            if nDim == 1:
                if np.size(data) == np.size(self.timeMid):
                    im, cb, fig = self.plotTime(data=data, spec=currSpec, unit=currUnit,
                                                targetUnit=targetUnit,
                                                ylim=ylim, xlim=xlim, logScale=logScale,
                                                labelFtSize=labelFtSize,
                                                labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == self.altSize:
                    im, cb, fig = self.plotAltitude(data=data, spec=currSpec, unit=currUnit,
                                                    targetUnit=targetUnit,
                                                    ylim=ylim, xlim=xlim, logScale=logScale,
                                                    labelFtSize=labelFtSize,
                                                    labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == self.latSize:
                    im, cb, fig = self.plotLatitude(data=data, spec=currSpec, unit=currUnit,
                                                    targetUnit=targetUnit,
                                                    ylim=ylim, xlim=xlim, logScale=logScale,
                                                    labelFtSize=labelFtSize,
                                                    labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == self.lonSize:
                    im, cb, fig = self.plotLongitude(data=data, spec=currSpec, unit=currUnit,
                                                     targetUnit=targetUnit,
                                                     ylim=ylim, xlim=xlim, logScale=logScale,
                                                     labelFtSize=labelFtSize,
                                                     labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                else:
                    logging.error('Species {:s} has wrong shape {:s} for 1-D plot'.format(spec, '{}'.format(np.shape(data))))
            elif nDim == 2:
                if np.shape(data) == self.zonalSize:
                    im, cb, fig = self.plotZonal(data=data, spec=currSpec, unit=currUnit,
                                                 targetUnit=targetUnit,
                                                 cmap=cmap, clim=clim, ylim=ylim, xlim=xlim,
                                                 show_colorbar=show_colorbar, logScale=logScale,
                                                 cbTitle=cbTitle, labelFtSize=labelFtSize,
                                                 labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == self.layerSize:
                    im, cb, fig = self.plotLayer(data=data, spec=currSpec, unit=currUnit,
                                                 targetUnit=targetUnit,
                                                 cmap=cmap, clim=clim, ylim=ylim, xlim=xlim,
                                                 show_colorbar=show_colorbar, logScale=logScale,
                                                 cbTitle=cbTitle, labelFtSize=labelFtSize,
                                                 labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == (len(self.timeMid),self.nLev):
                    im, cb, fig = self.plotTimeAlt(data=data, spec=currSpec, unit=currUnit,
                                                   targetUnit=targetUnit,
                                                   cmap=cmap, clim=clim, ylim=ylim, xlim=xlim,
                                                   logScale=logScale, labelFtSize=labelFtSize,
                                                   labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == (len(self.timeMid),self.nLat):
                    im, cb, fig = self.plotTimeLat(data=data, spec=currSpec, unit=currUnit,
                                                   targetUnit=targetUnit,
                                                   cmap=cmap, clim=clim, ylim=ylim, xlim=xlim,
                                                   logScale=logScale, labelFtSize=labelFtSize,
                                                   labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                elif np.shape(data) == (len(self.timeMid),self.nLon):
                    im, cb, fig = self.plotTimeLon(data=data, spec=currSpec, unit=currUnit,
                                                   targetUnit=targetUnit,
                                                   cmap=cmap, clim=clim, ylim=ylim, xlim=xlim,
                                                   logScale=logScale, labelFtSize=labelFtSize,
                                                   labelTickSize=labelTickSize, isDiff=isDiff)
                    RC = SUCCESS
                else:
                    logging.error('Species {:s} has wrong shape {:s} for 2-D plot'.format(spec, '{}'.format(np.shape(data))))
            else:
                logging.error('Data appears to be >2D. Not sure how to plot higher dimensional data')
            if RC == SUCCESS:
                imageDict[spec] = im
                cbDict[spec] = cb
                figDict[spec] = fig

        return imageDict, cbDict, figDict

    def plotZonal(self, data=None, spec=None,
                  cmap=None, clim=None, ylim=None, xlim=None,
                  show_colorbar=True, cbTitle=None, logScale=False,
                  labelFtSize=18, labelTickSize=18, isDiff=False,
                  unit=None, targetUnit=None, latTicks=np.array([-90,-60,-30,0,30,60,90])):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotZonal')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _usr_cmap      = 'viridis'
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        _latTickLabels = self.__getLatTickLabels(latTicks)

        # Set colormap
        if isDiff or (np.min(data) < 0) or (np.max(np.abs(data)) == 0):
            _usr_cmap = 'RdBu_r'
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')
        if cmap is not None:
            _usr_cmap = cmap

        # Plot data
        if not logScale:
            im = ax.pcolormesh(self.latEdge, self.pEdge, data * _convFactor, cmap=_usr_cmap)
        else:
            im = ax.pcolormesh(self.latEdge, self.pEdge, data * _convFactor, cmap=_usr_cmap,
                    norm=colors.LogNorm(vmin=np.min(data * _convFactor),
                                        vmax=np.max(data * _convFactor)))

        # Change color limits
        if _isNeg:
            im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
        if clim is not None:
            _tmpclim = np.array(im.get_clim())
            if np.isfinite(clim[0]):
                _tmpclim[0] = clim[0]
            if np.isfinite(clim[1]):
                _tmpclim[1] = clim[1]
            im.set_clim(_tmpclim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        # Invert pressure axis and log scale
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_xlim([-90, 90])
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)
        ax.set_ylabel('Pressure, hPa', fontsize=labelFtSize)
        ax.set_xlabel('Latitude', fontsize=labelFtSize)
        ax.set_xticks(latTicks)
        ax.set_xticklabels(_latTickLabels)
        if spec is not None:
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('{:s} zonally-averaged, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        if show_colorbar:
            cb = fig.colorbar(im, ax=ax, shrink=1.0, orientation='vertical', pad=0.04)
            if cbTitle is not None:
                cb.set_label(cbTitle, fontsize=labelFtSize)
            else:
                cb.set_label(_displayUnit, fontsize=labelFtSize)
            cb.ax.tick_params(axis='y', which='major', labelsize=labelTickSize)
        else:
            cb = None

        return im, cb, fig

    def plotLayer(self, data=None, spec=None,
                  cmap=None, clim=None, xlim=None, ylim=None,
                  show_colorbar=True, cbTitle=None, logScale=False,
                  labelFtSize=18, labelTickSize=18, isDiff=False,
                  unit = None, targetUnit=None,
                  lonTicks=np.array([-180,-120,-60,0,60,120,180]),
                  latTicks=np.array([-90,-60,-30,0,30,60,90]),
                  show_coastlines=True):

        if (not isinstance(data, np.ndarray))and (spec is None):
            logging.error('No data was passed to plotLayer')
        elif (not isinstance(data, np.ndarray))and (spec is not None):
            data = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8),
                               subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0.0)})

        # Initialize local variables
        _usr_cmap      = 'viridis'
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        #lonShift = 0
        #if (not self.isCentered) and (np.min(lonTicks) < 0):
        #    lonShift = 180
        _latTickLabels = self.__getLatTickLabels(latTicks)
        _lonTickLabels = self.__getLonTickLabels(lonTicks, shift=0, isCentered=self.isCentered)

        # Set colormap
        if isDiff or (np.min(data) < 0) or (np.max(np.abs(data)) == 0):
            _usr_cmap = 'RdBu_r'
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')
        if cmap is not None:
            _usr_cmap = cmap

        # Plot data
        if not logScale:
            im = ax.pcolormesh(self.lonEdge, self.latEdge, data * _convFactor, cmap=_usr_cmap)
        else:
            im = ax.pcolormesh(self.lonEdge, self.latEdge, data * _convFactor, cmap=_usr_cmap,
                    norm=colors.LogNorm(vmin=np.min(data * _convFactor),
                                        vmax=np.max(data * _convFactor)))

        # Change color limits
        if _isNeg:
            im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
        if clim is not None:
            _tmpclim = np.array(im.get_clim())
            if np.isfinite(clim[0]):
                _tmpclim[0] = clim[0]
            if np.isfinite(clim[1]):
                _tmpclim[1] = clim[1]
            im.set_clim(_tmpclim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        ax.set_ylim([-90, 90])
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        ax.set_ylabel('Latitude', fontsize=labelFtSize)
        ax.set_xlabel('Longitude', fontsize=labelFtSize)
        ax.set_yticks(latTicks)
        ax.set_yticklabels(_latTickLabels)
        ax.set_xticks(lonTicks)
        ax.set_xticklabels(_lonTickLabels)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)
        if spec is not None:
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('{:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)

        if show_colorbar:
            cb = fig.colorbar(im, ax=ax, shrink=1.0, orientation='vertical', pad=0.04)
            if cbTitle is not None:
                cb.set_label(cbTitle, fontsize=labelFtSize)
            else:
                cb.set_label(_displayUnit, fontsize=labelFtSize)
            cb.ax.tick_params(axis='y', which='major', labelsize=labelTickSize)
        else:
            cb = None

        if show_coastlines:
            ax.coastlines('50m')

        return im, cb, fig

    def plotTime(self, data=None, spec=None,
                 xlim=None, ylim=None,
                 logScale=False, isDiff=False,
                 labelFtSize=18, labelTickSize=18,
                 unit=None, targetUnit=None):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotTime')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        # Set colormap
        if isDiff or (np.min(data) < 0):
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')

        # Plot data
        im = ax.plot(self.timeMid, data * _convFactor)

        # Change axis limits
        if _isNeg:
            ax.set_ylim(np.array([-1,1])*np.max(np.abs(ax.get_ylim())))
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        if logScale:
            ax.set_yscale('log')
        ax.set_xlabel('Time', fontsize=labelFtSize)
        ax.set_xlim([self.timeMid[0], self.timeMid[-1]])
        if spec is not None:
            ax.set_ylabel('Global {:s}, {:s}'.format(spec, targetUnit), fontsize=labelFtSize)
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('{:s} globally-averaged, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        cb = None

        return im, cb, fig

    def plotTimeAlt(self, data=None, spec=None,
                    cmap=None, clim=None, xlim=None, ylim=None,
                    show_colorbar=True, cbTitle=None, logScale=False,
                    labelFtSize=18, labelTickSize=18, isDiff=False,
                    unit=None, targetUnit=None):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotTimeAlt')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _usr_cmap      = 'viridis'
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        # Set colormap
        if isDiff or (np.min(data) < 0) or (np.max(np.abs(data)) == 0):
            _usr_cmap = 'RdBu_r'
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')
        if cmap is not None:
            _usr_cmap = cmap

        # Plot data
        if not logScale:
           im = ax.pcolormesh(self.timeMid, self.pEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap)
        else:
           im = ax.pcolormesh(self.timeMid, self.pEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap,
                    norm=colors.LogNorm(vmin=np.min(data * _convFactor),
                                        vmax=np.max(data * _convFactor)))


        # Change axis limits
        if _isNeg:
            im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
        if clim is not None:
            _tmpclim = np.array(im.get_clim())
            if np.isfinite(clim[0]):
                _tmpclim[0] = clim[0]
            if np.isfinite(clim[1]):
                _tmpclim[1] = clim[1]
            im.set_clim(_tmpclim)
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        # Invert pressure axis and log scale
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylabel('Pressure, hPa', fontsize=labelFtSize)
        ax.set_xlabel('Time', fontsize=labelFtSize)
        ax.set_xlim([self.timeMid[0], self.timeMid[-1]])
        if spec is not None:
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Temporal altitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        if show_colorbar:
            cb = fig.colorbar(im, ax=ax, shrink=1.0, orientation='vertical', pad=0.04)
            if cbTitle is not None:
                cb.set_label(cbTitle, fontsize=labelFtSize)
            else:
                cb.set_label(_displayUnit, fontsize=labelFtSize)
            cb.ax.tick_params(axis='y', which='major', labelsize=labelTickSize)
        else:
            cb = None

        return im, cb, fig

    def plotTimeLat(self, data=None, spec=None,
                    cmap=None, clim=None, xlim=None, ylim=None,
                    show_colorbar=True, cbTitle=None, logScale=False,
                    labelFtSize=18, labelTickSize=18, isDiff=False,
                    unit=None, targetUnit=None,
                    latTicks=np.array([-90,-60,-30,0,30,60,90])):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotTimeLat')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _usr_cmap      = 'viridis'
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        _latTickLabels = self.__getLatTickLabels(latTicks)

        # Set colormap
        if isDiff or (np.min(data) < 0) or (np.max(np.abs(data)) == 0):
            _usr_cmap = 'RdBu_r'
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')
        if cmap is not None:
            _usr_cmap = cmap

        # Plot data
        if not logScale:
            im = ax.pcolormesh(self.timeMid, self.latEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap)
        else:
            im = ax.pcolormesh(self.timeMid, self.latEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap,
                    norm=colors.LogNorm(vmin=np.min(data * _convFactor),
                                        vmax=np.max(data * _convFactor)))

        # Change axis limits
        if _isNeg:
            im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
        if clim is not None:
            _tmpclim = np.array(im.get_clim())
            if np.isfinite(clim[0]):
                _tmpclim[0] = clim[0]
            if np.isfinite(clim[1]):
                _tmpclim[1] = clim[1]
            im.set_clim(_tmpclim)
        ax.set_ylim([-90, 90])
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)


        ax.set_ylabel('Latitude', fontsize=labelFtSize)
        ax.set_yticks(latTicks)
        ax.set_yticklabels(_latTickLabels)
        ax.set_xlabel('Time', fontsize=labelFtSize)
        ax.set_xlim([self.timeMid[0], self.timeMid[-1]])
        if spec is not None:
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Temporal latitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        if show_colorbar:
            cb = fig.colorbar(im, ax=ax, shrink=1.0, orientation='vertical', pad=0.04)
            if cbTitle is not None:
                cb.set_label(cbTitle, fontsize=labelFtSize)
            else:
                cb.set_label(_displayUnit, fontsize=labelFtSize)
            cb.ax.tick_params(axis='y', which='major', labelsize=labelTickSize)
        else:
            cb = None

        return im, cb, fig

    def plotTimeLon(self, data=None, spec=None,
                    cmap=None, clim=None, xlim=None, ylim=None,
                    show_colorbar=True, cbTitle=None, logScale=False,
                    labelFtSize=18, labelTickSize=18, isDiff=False,
                    unit=None, targetUnit=None,
                    lonTicks=np.array([-180,-120,-60,0,60,120,180])):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotTimeLon')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _usr_cmap      = 'viridis'
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        lonShift = 180
        _lonTickLabels = self.__getLonTickLabels(lonTicks + lonShift, isCentered=self.isCentered)

        # Set colormap
        if isDiff or (np.min(data) < 0) or (np.max(np.abs(data)) == 0):
            _usr_cmap = 'RdBu_r'
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')
        if cmap is not None:
            _usr_cmap = cmap

        # Plot data
        if not logScale:
            im = ax.pcolormesh(self.timeMid, self.lonEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap)
        else:
            im = ax.pcolormesh(self.timeMid, self.lonEdge, np.transpose(data) * _convFactor, cmap=_usr_cmap,
                    norm=colors.LogNorm(vmin=np.min(data * _convFactor),
                                        vmax=np.max(data * _convFactor)))

        # Change axis limits
        if _isNeg:
            im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
        if clim is not None:
            _tmpclim = np.array(im.get_clim())
            if np.isfinite(clim[0]):
                _tmpclim[0] = clim[0]
            if np.isfinite(clim[1]):
                _tmpclim[1] = clim[1]
            im.set_clim(_tmpclim)
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        ax.set_ylabel('Longitude', fontsize=labelFtSize)
        ax.set_yticks(lonTicks + lonShift)
        ax.set_yticklabels(_lonTickLabels)
        ax.set_xlabel('Time', fontsize=labelFtSize)
        ax.set_xlim([self.timeMid[0], self.timeMid[-1]])
        if spec is not None:
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Temporal longitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        if show_colorbar:
            cb = fig.colorbar(im, ax=ax, shrink=1.0, orientation='vertical', pad=0.04)
            if cbTitle is not None:
                cb.set_label(cbTitle, fontsize=labelFtSize)
            else:
                cb.set_label(_displayUnit, fontsize=labelFtSize)
            cb.ax.tick_params(axis='y', which='major', labelsize=labelTickSize)
        else:
            cb = None

        return im, cb, fig

    def plotAltitude(self, data=None, spec=None,
                     xlim=None, ylim=None, isDiff=False,
                     logScale=False, labelFtSize=18, labelTickSize=18,
                     unit=None, targetUnit=None):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotAltitude')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        # Set colormap
        if isDiff or (np.min(data) < 0):
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')

        # Plot data
        im = ax.plot(data * _convFactor, self.pMid)

        # Change axis limits
        if _isNeg:
            ax.set_xlim(np.array([-1,1])*np.max(np.abs(ax.get_xlim())))
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        # Invert pressure axis and log scale
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_ylabel('Pressure, hPa', fontsize=labelFtSize)
        if logScale:
            ax.set_xscale('log')
        if spec is not None:
            ax.set_xlabel('Mean {:s}, {:s}'.format(spec, targetUnit), fontsize=labelFtSize)
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Altitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        cb = None

        return im, cb, fig

    def plotLatitude(self, data=None, spec=None,
                     xlim=None, ylim=None, isDiff=False,
                     logScale=False, labelFtSize=18, labelTickSize=18,
                     unit=None, targetUnit=None,
                     latTicks=np.array([-90,-60,-30,0,30,60,90])):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotLatitude')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        _latTickLabels = self.__getLatTickLabels(latTicks)

        # Set colormap
        if isDiff or (np.min(data) < 0):
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')

        # Plot data
        im = ax.plot(self.lat, data * _convFactor)

        # Change axis limits
        ax.set_xlim([-90, 90])
        if _isNeg:
            ax.set_ylim(np.array([-1,1])*np.max(np.abs(ax.get_ylim())))
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        ax.set_xlabel('Latitude', fontsize=labelFtSize)
        ax.set_xticks(latTicks)
        ax.set_xticklabels(_latTickLabels)
        if logScale:
            ax.set_yscale('log')
        if spec is not None:
            ax.set_ylabel('Mean {:s}, {:s}'.format(spec, targetUnit), fontsize=labelFtSize)
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Latitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        cb = None

        return im, cb, fig

    def plotLongitude(self, data=None, spec=None,
                      xlim=None, ylim=None, isDiff=False,
                      logScale=False, labelFtSize=18, labelTickSize=18,
                      unit=None, targetUnit=None,
                      lonTicks=np.array([-180,-120,-60,0,60,120,180])):

        if (not isinstance(data, np.ndarray)) and (spec is None):
            logging.error('No data was passed to plotLongitude')
        elif (not isinstance(data, np.ndarray)) and (spec is not None):
            data     = self.data[spec]
            unit = self.unit[spec]

        fig, ax = plt.subplots(1, 1, figsize=(15,8))

        # Initialize local variables
        _isNeg         = False
        _isMixingRatio = True
        RC             = SUCCESS
        _convFactor    = 1.0E+00
        _displayUnit   = 'None'

        if unit is not None:
            _convFactor, _displayUnit, RC = self.__checkUnit(unit, targetUnit)

        if RC == WRONG_UNIT:
            targetUnit = unit
            logging.warning('Keeping {:s} as plotting unit'.format(unit))
            _isMixingRatio = False

        #lonShift = 0
        #if (not self.isCentered) and (np.min(lonTicks) < 0):
        lonShift = 180
        _lonTickLabels = self.__getLonTickLabels(lonTicks + lonShift, isCentered=self.isCentered)

        # Set colormap
        if isDiff or (np.min(data) < 0):
            _isNeg = True
            if logScale:
                logging.error('Logarithmic scale is True, but data is negative!')

        # Plot data
        im = ax.plot(self.lon, data * _convFactor)

        # Change axis limits
        ax.set_xlim(np.array([-180, 180]) + lonShift)
        if _isNeg:
            ax.set_ylim(np.array([-1,1])*np.max(np.abs(ax.get_ylim())))
        if ylim is not None:
            _tmplim = np.array(ax.get_ylim())
            if np.isfinite(ylim[0]):
                _tmplim[0] = ylim[0]
            if np.isfinite(ylim[1]):
                _tmplim[1] = ylim[1]
            ax.set_ylim(_tmplim)
        if xlim is not None:
            _tmplim = np.array(ax.get_xlim())
            if np.isfinite(xlim[0]):
                _tmplim[0] = xlim[0]
            if np.isfinite(xlim[1]):
                _tmplim[1] = xlim[1]
            ax.set_xlim(_tmplim)

        if _isMixingRatio:
            Min  = np.min(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMin = 'ppbv'
            if Min < 1.0E+00:
                Min *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMin = 'pptv'
            elif Min > 1.0E+03:
                Min *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMin = 'ppmv'
            Max  = np.max(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMax = 'ppbv'
            if Max < 1.0E+00:
                Max *= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMax = 'pptv'
            elif Max > 1.0E+03:
                Max *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMax = 'ppmv'
            Mean = np.mean(data * _convFactor) * self.possUnit[targetUnit] / self.possUnit['ppbv']
            uMean= 'ppbv'
            if Mean < 1.0E+00:
                Mean*= self.possUnit['ppbv'] / self.possUnit['pptv']
                uMean= 'pptv'
            elif Mean > 1.0E+03:
                Mean *= self.possUnit['ppbv'] / self.possUnit['ppmv']
                uMean = 'ppmv'
        else:
            Min  = np.min(data * _convFactor)
            uMin = targetUnit
            Max  = np.max(data * _convFactor)
            uMax = targetUnit
            Mean = np.mean(data * _convFactor)
            uMean= targetUnit
        uMin  = self.__fancyUnit(uMin)
        uMax  = self.__fancyUnit(uMax)
        uMean = self.__fancyUnit(uMean)

        ax.set_xlabel('Longitude', fontsize=labelFtSize)
        ax.set_xticks(lonTicks + lonShift)
        ax.set_xticklabels(_lonTickLabels)
        if logScale:
            ax.set_yscale('log')
        if spec is not None:
            ax.set_ylabel('Mean {:s}, {:s}'.format(spec, targetUnit), fontsize=labelFtSize)
            _fancySpec = self.__fancySpec(spec)
            ax.set_title('Longitudinal variations of {:s}, {:s}\nMin: {:3.2e} {:s}, Max: {:3.2e} {:s}, Mean: {:3.2} {:s}'.
                     format(_fancySpec, _displayUnit, Min, uMin, Max, uMax, Mean, uMean),
                     fontsize=labelFtSize)
        ax.tick_params(axis='both', which='major', labelsize=labelTickSize)

        cb = None

        return im, cb, fig

    def __checkUnit(self, currUnit, targetUnit):

        _convFactor  = 1.0E+00
        _displayUnit = 'Empty unit'
        RC = SUCCESS

        if currUnit is None:
            RC = WRONG_UNIT
        if (targetUnit is not None) and (currUnit is None):
            print('Required targetUnit = {:s}'.format(targetUnit))
            raise ValueError('Could not figure out current unit...')

        if (targetUnit is not None) and (targetUnit not in self.possUnit):
            logging.warning('Could not parse targetUnit: {:s}'.format(targetUnit))
            RC = WRONG_UNIT
            _displayUnit = currUnit
        elif (currUnit is not None) and (currUnit not in self.possUnit):
            logging.warning('Could not parse currUnit: {:s}'.format(currUnit))
            RC = WRONG_UNIT
            _displayUnit = currUnit
        else:
            if currUnit is not None:
                _displayUnit = currUnit
            if targetUnit is not None:
                _displayUnit = targetUnit
                _convFactor = self.possUnit[currUnit]/self.possUnit[targetUnit]

        _displayUnit = self.__fancyUnit(_displayUnit)

        return _convFactor, _displayUnit, RC

    def __fancyUnit(self, unit):

        for superscript in ['1', '2', '3', '4']:
            unit = unit.replace(superscript, '$^' + superscript + '$')

        return unit

    def __fancySpec(self, spec):

        for subscript in ['1', '2', '3', '4', '5', '6', '7', '8', '9']:
            spec = spec.replace(subscript, '$_' + subscript + '$')

        return spec

    def __getLatTickLabels(self, latTicks):
        # Latitude ticks
        latTickLabels = []
        for tick in latTicks:
            tick_str = '{:d}'.format(np.abs(tick))
            if tick < -0.0001:
                tick_label = tick_str + '$^\circ$' + 'S'
            elif tick > 0.0001:
                tick_label = tick_str + '$^\circ$' + 'N'
            else:
                tick_label = 'Eq'
            latTickLabels.append(tick_label)
        return latTickLabels

    def __getLonTickLabels(self, lonTicks, shift=0, isCentered=True):
        # Longitude ticks
        lonTickLabels = []
        for tick in lonTicks:
            _tick = tick+shift
            if _tick > 180:
                _tick = _tick - 360
            tick_str = '{:d}'.format(np.abs(_tick))
            _threshold = 0
#             if not isCentered:
#                 _threshold = 180
            if _tick < _threshold - 0.0001:
                tick_label = tick_str + '$^\circ$' + 'W'
            elif _tick > _threshold + 0.0001:
                tick_label = tick_str + '$^\circ$' + 'E'
            else:
                tick_label = '0'
            lonTickLabels.append(tick_label)
        return lonTickLabels

