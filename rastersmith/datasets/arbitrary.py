from __future__ import division, print_function

import datetime
import numpy as np
from osgeo import gdal,osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from pyproj import Proj
from scipy import ndimage
import xarray as xr

from ..core import core
from ..core import utils

class Arbitrary(core.Raster):
    def __init__(self):
        core.Raster.__init__(self)

        return

    @classmethod
    def read(cls,infile,yxzAxes=[0,1,None],bandNames=None,time=None,
             preprocessing=None,preprocessingArgs=None,preprocessingKwargs=None,
             attrs=None,sensor='raster',crs='4326'):
        cls.crs = {'init':'epsg:{}'.format(crs)}
        cls.sensor= sensor

        ds = gdal.Open(infile,GA_ReadOnly)

        if ds.GetDriver().ShortName == 'GTiff':
            nBands = ds.RasterCount
            yDim = ds.RasterYSize
            xDim = ds.RasterXSize
            dims = (yDim,xDim)

            srs = osr.SpatialReference()
            srs.ImportFromWkt(ds.GetProjection())

            projStr = srs.ExportToProj4()
            proj = Proj(projStr)

            mask = np.ones([yDim,xDim])

            if type(bandNames) != list:
                bandNames = []

            if len(bandNames) < nBands:
                bandNames = [bandNames[b] if b < len(bandNames) else 'b{}'.format(b+1) for b in range(nBands)]

            elif len(bandNames) > nBands:
                bandNames = bandNames[:nBands]

            else:
                pass

            args = [[ds,band] for band in range(1,nBands+1)]

            dataarr = list(map(lambda x: cls._readBand(x,preprocessing,preprocessingArgs),args))
            dataarr.append(mask)
            bandNames.append('mask')

            dataarr = utils.formatDataarr(dataarr)

            west, xres, xskew, north, yskew, yres  = ds.GetGeoTransform()
            east = west + (xDim * xres)
            south = north + (yDim * yres)

            extent = (west,south,east,north)

            lons,lats = cls.geoGrid(extent,dims,projStr,wgsBounds=proj.is_latlong())

            if time:
                if type(time) == str:
                    dt= datetime.datetime.strptime(time,'%Y-%m-%d')
                elif type(time) == datetime.datetime:
                    dt = time
                else:
                    raise ValueError('Time either needs to be a datetime object or string in format of YYYY-MM-DD')
            else:
                dt = datetime.datetime(1970,1,1,0,0,0,0)

        # elif ds.GetDriver().ShortName == 'HDF5':
        #     subdata = ds.GetSubDatasets()

        else:
            raise NotImplementedError('Input dataset was not able to be read in')

        coords = {'z': range(dataarr.shape[2]),
                  'lat':lats[:,0],
                  'lon':lons[0,:],
                  'band':(bandNames),
                  'time':([np.datetime64(dt)])}

        dims = ('lat','lon','z','band','time')

        attrs = {'projStr': projStr,
                 'bandNames':tuple(bandNames),
                 'extent':extent,
                 'date':dt,
                 'resolution':(yres,xres)
                 }

        ds = xr.DataArray(dataarr,coords=coords,dims=dims,attrs=attrs,name=cls.sensor)

        return ds

    @staticmethod
    def _readBand(args,func=None,funcArgs=None,funcKwargs=None):

        ds,bandIdx = args

        band = ds.GetRasterBand(bandIdx)
        var = BandReadAsArray(band)

        if func:
            if funcArgs:
                result = func(var,*funcArgs)
            else:
                result = func(var)
        else:
            result = var


        return result

    # def _readSubData(self,):
    #
    #     return
