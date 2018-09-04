import glob
import itertools
import xarray as xr

from .atms import Atms
from .viirs import Viirs
from .landsat import Landsat
from .arbitrary import Arbitrary
from ..core import utils

def readRaster(fileName, sensor='infer'):
    if sensor == 'infer':
        sensor = _inferSensor(fileName)

    if sensor == 'viirs':
        out = Viirs.read(fileName)

    elif sensor == 'landsat8':
        out = Landsat.read(fileName,sensor=sensor)

    elif sensor == 'landsat457':
        out = Landsat.read(fileName,sensor=sensor)

    elif sensor == 'atms':
        out = Atms.read(fileName)

    elif sensor == 'unknown':
        try:
            out = Arbitrary.read(fileName)
        except:
            raise NotImplementedError('Specified raster file was not able to be read.\
                                       See documentation for data and sensors available to use.')
    else:
        raise NotImplementedError('Specified sensor was not able to be read.\
                                   See documentation for data and sensors available to use.')

    return out

def readCollection(directory,sensor='infer',spatialMosaic=True):
    possibleFormats = ['.h5','MTL.txt']
    fList = glob.glob(directory+'*.*')
    fList = [f for f in fList if any(format in f for format in possibleFormats)]

    rasters = list(map(readRaster,fList))
    xrRasters = [rast.raster for rast in rasters]

    if sensor == 'infer':
        sensor = _inferSensor(fList[0])

    if ('viirs' in sensor) | ('atms' in sensor):
        inRasters = list(map(lambda x:[x.attrs['date'].strftime("%Y-%m-%d"),x],xrRasters))

    elif 'landsat' in sensor:
        inRasters = list(map(lambda x:[0,x],xrRasters))

    if spatialMosaic:
        result = utils.stich(inRasters)
    else:
        result = xr.concat(xrRasters,dim='time')

    return result


def _inferSensor(fileName):

    if ('LC08' in fileName):
        sensorName = 'landsat8'

    elif ('LT05' in fileName) | ('LT04' in fileName) | ('LE07' in fileName):
        sensorName = 'landsat457'

    elif ('VNP' in fileName):
        sensorName = 'viirs'

    elif ('ATM' in fileName):
        sensorName = 'atms'

    else:
        sensorName = 'unknown'

    return sensorName
