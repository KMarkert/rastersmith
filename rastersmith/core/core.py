# base functionality package

from __future__ import division, print_function

#import glob
import math
import datetime
from itertools import groupby

import numpy as np
import xarray as xr
import bottleneck as bn

from osgeo import gdal,osr
from pyproj import Proj, transform
from PIL import Image, ImageDraw

import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


from . import utils
from . import countries



class Grid(object):
    def __init__(self,crs='4326',region=(92.3032344909, 9.93295990645, 101.180005324, 28.335945136)
                     ,resolution=500,noData=0,country=None):

        if country:
            if crs != '4326':
                raise ValueError('User defined crs must be EPSG:4326 when using predifined country bounding boxes')
            else:
                self.west,self.south,self.east,self.north = countries.bounding_boxes[country][1]
        else:
            #! Need to add in error handler for correct bounding box format !#
            self.west,self.south,self.east,self.north = region

        if '4326' in crs:
            midPoint = [bn.nanmean([self.north,self.south]),
                        bn.nanmean([self.east,self.west])]
            spacing = utils.meters2dd(midPoint,resolution)

        elif type(resolution) == list:
            spacing = resolution[0]

        else:
            spacing = resolution,resolution

        self.lons = np.arange(self.west,self.east,spacing)
        self.lats = np.arange(self.south,self.north,spacing)

        self.xx,self.yy = np.meshgrid(self.lons,self.lats)

        self.nominalResolution = (spacing)
        self.dims = self.xx.shape

        return


# @geopandas
class Geometry(object):
    def __init__(self,fileName):
        return


@xr.register_dataarray_accessor('raster')
class Raster(object):
    def __init__(self,xarray_obj):
        self._obj = xarray_obj

        return

    @property
    def gt(self):
        rasterobj = self._obj
        ulx = float(rasterobj.coords['lon'].min().values)
        uly = float(rasterobj.coords['lat'].max().values)

        pySize,pxSize = rasterobj.attrs['resolution']

        if pySize > 0:
            pySize = pySize * -1

        gt = (ulx,pxSize,0,uly,0,pySize)

        return gt

    @staticmethod
    def geoGrid(extent,dims,nativeProj,wgsBounds=False):

        west, south, east, north = extent

        gcsProj = Proj(init='epsg:4326')
        native = Proj(nativeProj)

        gcs = native.is_latlong()

        if wgsBounds and ~gcs:
            llx,lly = transform(gcsProj,native,west,south)
            urx,ury = transform(gcsProj,native,east,north)
        else:
            llx,lly = west,south
            urx,ury = east,north

        yCoords = np.linspace(lly,ury,dims[0],endpoint=False)[::-1]
        xCoords = np.linspace(llx,urx,dims[1],endpoint=False)

        xx,yy = np.meshgrid(xCoords,yCoords)

        return xx,yy

    @staticmethod
    def _extractBits(image,start,end):
        """Helper function to convert Quality Assurance band bit information to flag values

        Args:
            image (ndarray): Quality assurance image as a numpy array
            start (int): Bit position to start value conversion
            end (int): Bit position to end value conversion

        Returns:
            out (ndarray): Output quality assurance in values from bit range
        """

        pattern = 0;
        for i in range(start,end+1):
            pattern += math.pow(2, i)

        bits = image.astype(np.uint16) & int(pattern)
        out = bits >> start

        return out


    @classmethod
    def select(cls,bandList,newNames=None):
        if (type(bandList) != list) and (type(bandList) == str):
            bandList = [bandList]

        bandNames = cls._obj.coords['band'].values

        newBands = [new for new in bandList if any(b in new for b in bandNames)]

        out = cls._obj.copy()

        out.raster = out.sel(band=newBands)

        if newNames:
            out.coords['band'] = newNames

        return out


    def clip(self,geom):

        return

    def normalizedDifference(self,band1=None,band2=None,outBandName='nd',appendTo=None):
        rasterobj = self._obj

        if band1 and band2:
            nd = (rasterobj.sel(band=band1) - rasterobj.sel(band=band2)) / \
                 (rasterobj.sel(band=band1) + rasterobj.sel(band=band2))

            nd = nd.expand_dims('band')
            nd.coords['band'] = [outBandName]
            nd = nd.transpose('lat','lon','z','band','time')

        if appendTo:
            out = xr.concat([appendTo,nd],dim='band')
        else:
            out = nd

        return out

    def updateMask(self,maskDa,applyMask=True):
        out = self._obj

        if type(maskDa) == np.ndarray:
            yCoords = out.coords['lat'].values
            xCoords = out.coords['lon'].values
            maskDa = xr.DataArray(maskDa,coords=[yCoords,xCoords],dims=['lat','lon'])

        out.values[:,:,:,-1,:] = self._obj.sel(band='mask').astype(np.bool)\
                                            & maskDa.astype(np.bool)
        if applyMask:
            out = out.raster.applyMask()

        return out

    def applyMask(self):
        out= self._obj.where(self._obj.sel(band='mask')>0)
        return out


    def unmask(self,value=None):
        bNames = [i for i in self.bands.keys() if i != 'mask']
        mask = self.bands['mask']

        out = self.copy()

        if value:
            for i in bNames:
                out.bands[i][np.where(mask==0)] = value

        else:
            for i in bNames:
                out.bands[i] = self.bands[i].data

        return out


    def showMap(self,band=None,showCountries=True,cmap='viridis'):

        ax = plt.axes(projection=ccrs.PlateCarree())

        if band:
            raster = self._obj.sel(band=band)
        else:
            if 'band' in list(self._obj.coords.keys()):
                raster = self._obj.isel(band=0)
            else:
                raster = self._obj

        raster.plot(ax=ax,robust=True,cmap=cmap)
        ax.coastlines()

        ax.add_feature(cartopy.feature.BORDERS)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5,linestyle=':')

        gl.xlabels_top = False
        gl.ylabels_right = False

        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        ax.set_xlabel('')
        ax.set_ylabel('')

        plt.show()

        return


    def writeGeotiff(self,path,prefix):
        rasterobj = self._obj
        drv = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromProj4(rasterobj.attrs['projStr'])

        y = len(rasterobj.coords['lat'])
        x = len(rasterobj.coords['lon'])
        bands = len(rasterobj.coords['band'])

        t = rasterobj.attrs['date'].strftime("%Y%m%d")

        outDs = drv.Create(path + prefix + '_' + t + '.tif',
                           x,y,bands,
                           gdal.GDT_Int16
                           )

        flip = rasterobj.attrs['resolution'][0] < 0

        for b in range(bands):
            band = outDs.GetRasterBand(b+1)
            band.WriteArray(rasterobj[:,:,0,b,0].values)
            band.SetNoDataValue(0)
            band = None

        outDs.SetGeoTransform(rasterobj.raster.gt)

        outDs.SetProjection(srs.ExportToWkt())

        outDs.FlushCache()

        return
