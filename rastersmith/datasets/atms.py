# ATMS package

from __future__ import division, print_function

import os
import datetime
import numpy as np
from osgeo import gdal
import xarray as xr

from ..core import core
from ..core import utils
from ..core import mapping


dir = os.path.dirname(os.path.abspath(__file__))

class Atms(core.Raster):
    def __init__(self):
        core.Raster.__init__(self)

        return

    @classmethod
    def read(cls,infile,ingeo,sensor='atms',crs='4326'):
        cls.crs = {'init':'epsg:{}'.format(crs)}
        cls.sensor= sensor

        trees = ['//All_Data/ATMS-SDR_All/',
                 '//All_Data/ATMS-SDR-GEO_All/'
                ]

        fields = [['GainCalibration','BrightnessTemperature',],
                  ['Latitude','Longitude','SatelliteZenithAngle']
                 ]

        base = 'HDF5:"{0}":{1}{2}'

        gainDs = gdal.Open(base.format(infile,trees[0],fields[0][0]))
        gainCal = gainDs.ReadAsArray()

        gainDs = None

        ds = gdal.Open(base.format(infile,trees[0],fields[0][1]))
        data = ds.ReadAsArray()
        metadata = ds.GetMetadata()

        ds = None

        bandNames = ['C{}'.format(i) for i in range(1,23)]
        bandIdx = range(22)

        nativeCrs = {'init':'epsg:4326'}
        projStr = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

        geos = list(map(cls._readGeo,[[ingeo,trees[1],field] for field in fields[1]]))
        lats,lons,view = geos

        shape = view.shape

        mask = view<50 | np.ones(shape).astype(np.bool)

        dataarr = list(map(lambda x: cls._calibrate(data,gainCal,x),bandIdx))
        dataarr.append(mask)

        bandNames.append('mask')

        dataarr = utils.formatDataarr(dataarr)

        res = utils.meters2dd((lats.mean(),lons.mean()),scale=15000)
        dataarr,gridmask,newX,newY = mapping.swath2grid(dataarr,lons,lats,resolution=res,method='cubic')

        bb = (newX.min(),newY.min(),newX.max(),newY.max())

        data = None

        fileComps = infile.split('_')
        timeComps = fileComps[2][1:] + fileComps[3][1:] + 'UTC'

        dt = datetime.datetime.strptime(timeComps,'%Y%m%d%H%M%S%f%Z')

        coords = {'z': range(dataarr.shape[2]),
                  'lat':newY,
                  'lon':newX,
                  'band':(bandNames),
                  'time':([np.datetime64(dt)])}

        dims = ('lat','lon','z','band','time')

        attrs = {'projStr': projStr,
                 'bandNames':tuple(bandNames),
                 'extent': bb,
                 'date':dt,
                 'units': 'brightness_temperature',
                 'scale_factor': 1,
                 'add_offset': 0,
                 'resolution':res
                 }

        ds = xr.DataArray(dataarr,coords=coords,dims=dims,attrs=attrs,name=cls.sensor)

        ds = ds.raster.updateMask(gridmask)

        return ds


    @staticmethod
    def _readGeo(args):
        base = 'HDF5:"{0}":{1}{2}'

        ingeo, tree, field = args

        band = gdal.Open(base.format(ingeo,tree,field))
        geo = band.ReadAsArray()

        band = None
        data = None

        return geo

    @classmethod
    def _calibrate(self,data,gc, idx):
        calData = data[:,:,idx]

        BTr = np.zeros_like(calData).astype(np.uint16)
        for i in range(BTr.shape[0]):
            BTr[i,:] = (calData[i,:].astype(int) * gc[i,idx]) * 0.0001

        return BTr

    @staticmethod
    def _mlc(dataarr,keys):
        df = np.load(os.path.join(dir,'atms_mlc_coeffs.npy'))
        calFactors = dict(df.item())

        Btr = dataarr

        cal_mlc = np.zeros_like(Btr)

        for c in range(len(keys)):
            chnl = keys[c]
            if "mask" not in chnl:
                tbi = calFactors[chnl]['Tbi']
                preds = calFactors[chnl]['predictors']
                popt = calFactors[chnl]['popt']
                for i in range(cal_mlc.shape[0]):
                    factor = np.zeros(cal_mlc.shape[1])
                    for j in range(popt.size):
                        Tbij = calFactors[keys[j]]['Tbij']
                        factor = factor + (popt[j]*(Btr[i,:,0,preds[j],0]-Tbij))
                    cal_mlc[i,:,0,c,0] = tbi + factor
            else:
                cal_mlc[:,:,0,c,0] = Btr[:,:,0,c,0]

        return cal_mlc
