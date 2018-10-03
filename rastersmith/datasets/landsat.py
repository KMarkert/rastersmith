from __future__ import division, print_function

import os
import datetime
import numpy as np
from osgeo import gdal
from scipy import ndimage
import xarray as xr
import math

from ..core import core
from ..core import utils


class Landsat(core.Raster):
    def __init__(self):
        core.Raster.__init__(self)

        return

    @classmethod
    def read(cls,metaFile,sensor='landsat',crs='4326'):
        cls.crs = {'init':'epsg:{}'.format(crs)}
        cls.sensor= sensor

        metadata = cls._parseMetadata(metaFile)

        path = os.path.split(metaFile)[0]+'/'

        north = max(float(metadata['CORNER_UR_PROJECTION_Y_PRODUCT']),
                    float(metadata['CORNER_UL_PROJECTION_Y_PRODUCT']))
        south = min(float(metadata['CORNER_LR_PROJECTION_Y_PRODUCT']),
                    float(metadata['CORNER_LL_PROJECTION_Y_PRODUCT']))
        east  = max(float(metadata['CORNER_UR_PROJECTION_X_PRODUCT']),
                    float(metadata['CORNER_LR_PROJECTION_X_PRODUCT']))
        west  = min(float(metadata['CORNER_UL_PROJECTION_X_PRODUCT']),
                    float(metadata['CORNER_LL_PROJECTION_X_PRODUCT']))

        extent = (west,south,east,north)

        projStr = '+proj=utm +zone={} +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs'\
                    .format(metadata['UTM_ZONE'])

        bandList = ['QUALITY']
        for i in range(2,8):
            bandList.append(str(i))

        'BQA',

        qa = cls._readBand(path,metadata,'QUALITY',calibrate=False)
        mask = cls._extractBits(qa,4,4)

        shape = mask.shape

        if ('8' in sensor) | ('LC08' in metaFile):
            bandNames = ['B1','B2','B3','B4','B5','B6','B7']
            bandValue = range(1,8)
        else:
            bandNames = ['B1','B2','B3','B4','B5','B7']
            bandValue = range(1,8).remove(7)
        bandDim = len(bandNames)+1

        bands = {}

        dataarr = list(map(lambda k: cls._readBand(path,metadata,k),bandValue))
        dataarr.append(mask)
        bandNames.append('mask')

        dataarr = utils.formatDataarr(dataarr)

        lons,lats = cls.geoGrid(extent,shape,projStr,wgsBounds=False)

        date = '{0} {1}{2}'.format(metadata['DATE_ACQUIRED'],metadata['SCENE_CENTER_TIME'][:-3], ' UTC')
        dt = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')

        coords = {'z': range(dataarr.shape[2]),
                  'lat':lats[:,0],
                  'lon':lons[0,:],
                  'band':(bandNames),
                  'time':([np.datetime64(dt)])}

        dims = ('lat','lon','z','band','time')

        attrs = {'projStr': projStr,
                 'bandNames':tuple(bandNames),
                 'extent':(west,south,east,north),
                 'date':dt,
                 'units': 'percent_reflectance',
                 'scale_factor': 10000,
                 'add_offset': 0,
                 'resolution':30
                 }

        ds = xr.DataArray(dataarr,coords=coords,dims=dims,attrs=attrs,name=cls.sensor)

        return ds

    @classmethod
    def _readBand(cls,path,metadata,key,calibrate=True):
        bandKey = 'FILE_NAME_BAND_{}'

        name = metadata[bandKey.format(key)]
        ds = gdal.Open('{0}{1}'.format(path,name))

        tmp = ds.ReadAsArray()

        if calibrate:
            out = cls._calibrate(tmp,key,metadata)
        else:
            out = tmp

        ds = None

        return out

    @classmethod
    def _calibrate(cls,arr,key,metadata,toa=True):
        gainKey = 'RADIANCE_MULT_BAND_{}'#'REFLECTANCE_MULT_BAND_{}'
        biasKey = 'RADIANCE_ADD_BAND_{}'

        sZenith = 90 - float(metadata['SUN_ELEVATION'])
        dSun = float(metadata['EARTH_SUN_DISTANCE'])

        gain = float(metadata[gainKey.format(key)])
        bias = float(metadata[biasKey.format(key)])

        rad = ((gain * arr) + bias)
        rad[np.where(rad<0)] = 0

        if toa:
            esun = cls._esunLookup(key,metadata['SENSOR_ID'])
            out = (math.pi * rad * (dSun**2)) / (esun * math.cos(math.radians(sZenith)))

        else:
            out = rad

        return out * 10000

    @staticmethod
    def _parseMetadata(metadata):
        with open(metadata,'r') as f:
            data = f.read()

        split_metadata = data.split('\n')

        output = {}
        for x in split_metadata:
            if "=" in x:
                line = x.split("=")
                output[line[0].strip()] = line[1].strip()
                clean_output = {key: item.strip('"') for key, item in output.items()}

        return clean_output


    @staticmethod
    def _esunLookup(key,sensor):
        esuns = {"TM":{'B1':1983.0,'B2':1796.0,'B3':1536.0,'B4':1031.0,'B5':220.0,'B7':83.44},
                "ETM":{'B1':1997.0,'B2':1812.0,'B3':1533.0,'B4':1039.0,'B5':230.8,'B7':84.90},
                "OLI_TIRS":{'B1':1895.33,'B2':2004.57,'B3':1820.75,'B4':1549.49,'B5':951.76,'B6':247.55,'B7':85.46}
                }
        return esuns[sensor]['B'+str(key)]
