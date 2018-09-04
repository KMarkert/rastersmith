# VIIRS packge

from __future__ import division, print_function

import datetime
import numpy as np
from osgeo import gdal
from scipy import ndimage
import xarray as xr


from ..core import core

class Viirs(core.Raster):

    def __init__(self):
        core.Raster.__init__(self)

        return

    @classmethod
    def read(cls,infile,sensor='viirs',crs='4326'):
        cls.crs = {'init':'epsg:{}'.format(crs)}
        cls.sensor= sensor

        tree = '//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/'
        field = 'SurfReflect_{0}{1}_1'
        base = 'HDF5:"{0}":{1}{2}'

        m = [i for i in range(12) if i not in [0,6,9]]
        i = [i for i in range(1,4)]
        bands = [m,i]
        flatBands = [item for sublist in bands for item in sublist]

        res = ['1km','500m']
        mode = ['M','I']

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',1)))
        metadata = band.GetMetadata()
        cloudQA = cls._extractBits(band.ReadAsArray(),2,3)
        hiresCloudQA = ndimage.zoom(cloudQA,2,order=0)
        band = None

        band = gdal.Open(base.format(infile,tree.format('1km'),field.format('QF',2)))
        shadowQA = cls._extractBits(band.ReadAsArray(),3,3)
        hiresShadowQA = ndimage.zoom(shadowQA,2,order=0)

        mask = ~(hiresCloudQA>0)&(hiresShadowQA<1)

        shape = mask.shape

        east,west = float(metadata['EastBoundingCoord']), float(metadata['WestBoundingCoord'])
        north,south = float(metadata['NorthBoundingCoord']), float(metadata['SouthBoundingCoord'])

        extent = (west,south,east,north)

        bandNames = ['{0}{1}'.format(mode[i],bands[i][j]) \
                        for i in range(len(res)) \
                        for j in range(len(bands[i]))
                    ]

        subdata = [[infile,res[i],mode[i],bands[i][j]] \
                        for i in range(len(res)) \
                        for j in range(len(bands[i]))
                  ]

        dataarr = list(map(cls._readBand,subdata))
        dataarr.append(mask)
        bandNames.append('mask')

        dataarr = np.moveaxis(np.array(dataarr),0,2)[:,:,:,np.newaxis]
        dataarr = np.moveaxis(dataarr,2,3)[:,:,:,:,np.newaxis]

        nativeCRS = {'init':'epsg:6974'}
        proj = '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'

        lons,lats = cls._geoGrid(extent,shape,proj,wgsBounds=True)

        gt = None

        date = '{0}{1}{2}'.format(metadata['RangeBeginningDate'],metadata['RangeBeginningTime'],' UTC')

        dt = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f %Z')

        coords = {'y': range(dataarr.shape[0]),
                  'x': range(dataarr.shape[1]),
                  'z': range(dataarr.shape[2]),
                  'lat':(['y','x'],lats),
                  'lon':(['y','x'],lons),
                  'band':(bandNames),
                  'time':([np.datetime64(dt)])}

        dims = ('y','x','z','band','time')

        attrs = {'nativeCrs':{'init':'epsg:6974'},
                 'projStr': proj,
                 'bandNames':tuple(bandNames),
                 'extent':(west,south,east,north),
                 'date':dt
                 }

        ds = xr.DataArray(dataarr,coords=coords,dims=dims,attrs=attrs,name=cls.sensor) \
               .chunk({'x': 1000,'y': 1000})

        return ds

    @staticmethod
    def _readBand(subdata):
        tree = '//HDFEOS/GRIDS/VNP_Grid_{}_2D/Data_Fields/'
        field = 'SurfReflect_{0}{1}_1'
        base = 'HDF5:"{0}":{1}{2}'


        infile, r, m, b = subdata
        subdataset = base.format(infile,tree.format(r),field.format(m,b))

        band = gdal.Open(subdataset)
        if m == 'M':
            data = ndimage.zoom(band.ReadAsArray(),2,order=0)

        else:
            data = np.array(band.ReadAsArray())

        return data.astype(np.int16)
