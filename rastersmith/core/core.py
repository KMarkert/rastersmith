# base functionality package

from __future__ import division, print_function

#import glob
import copy
import math
import pickle
import datetime
from itertools import groupby

import numpy as np
import xarray as xr
# import geopandas as gpd
from osgeo import gdal,osr
from pyproj import Proj, transform

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
            midPoint = [np.mean([self.north,self.south]),
                        np.mean([self.east,self.west])]
            ySpacing, xSpacing = utils.meters2dd(midPoint,resolution)

        elif type(resolution) == list:
            xSpacing = resolution[0]
            ySpacing = resolution[1]

        else:
            ySpacing,xSpacing= resolution,resolution

        lons = np.arange(self.west,self.east,xSpacing)
        lats = np.arange(self.south,self.north,ySpacing)

        self.xx,self.yy = np.meshgrid(lons,lats)

        self.nominalResolution = (xSpacing,ySpacing)
        self.dims = self.xx.shape

        return

    def _type(self):
        return 'grid'

    def mapRaster(self,raster,interpMethod='linear'):

        rasterLons = raster.coords['lon']
        rasterLats = raster.coords['lat']


        xSelect = (rasterLons>self.west) & (rasterLons<self.east)
        ySelect = (rasterLats>self.south) & (rasterLats<self.north)
        spatialSelect = ySelect & xSelect

        idx = np.where(spatialSelect == True)

        # Format geolocation coordinates for gridding
        pts = np.zeros((idx[0].size,2))
        pts[:,0] = rasterLons[idx].ravel()
        pts[:,1] = rasterLats[idx].ravel()

        bNames = list(raster.bands.keys())

        qualityBands = np.zeros([self.dims[0],self.dims[1],len(bNames)])

        out = raster.copy()

        for i in range(len(bNames)):
            if 'mask' in bNames[i]:
                iMethod = 'nearest'
            else:
                iMethod = interpMethod

            # Regrid data to common grid
            out.bands[bNames[i]] = interpolate.griddata(pts,
                                     raster.bands[bNames[i]][idx].ravel(),
                                     (self.xx,self.yy), method=iMethod,
                                     )

            qualityBands[:,:,i] = (out.bands[bNames[i]] < 16000) & \
                                  (out.bands[bNames[i]] > -1)
            if i == 0:
                interpMask = np.isnan(out.bands[bNames[i]])

        qualityMask = np.min(qualityBands,axis=2).astype(np.bool)
        out.bands['mask'] = out.bands['mask'] & ~interpMask & qualityMask
        out.updateMask()

        out.coords['Lon'],out.coords['Lat'] = self.xx,self.yy
        out.extent = (self.west,self.south,self.east,self.north)

        out.gt = out._getGt(self.north,self.west,self.nominalResolution)

        return out

# @geopandas
class Geometry(object):
    def __init__(self,fileName):
        return


@xr.register_dataarray_accessor('raster')
class Raster(object):
    def __init__(self,sensor=None,crs='4326'):
        self.sensor = sensor
        self.crs = {'init':'epsg:{}'.format(crs)}

        return

    # @classmethod
    # def copy(cls):
    #     return copy.deepcopy(cls)


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
    def _geoGrid(cls,extent,dims,nativeProj,wgsBounds=None):

        west, south, east, north = extent

        gcsProj = Proj(cls.crs)
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

        if ~gcs:
            lons,lats = transform(native,gcsProj,xx,yy)
        else:
            lons,lats = xx, yy

        return lons,lats

    @classmethod
    def _getGt(self,north,west,gridSize,projStr=None):
        if projStr:
            outProj = Proj(projStr)
            inProj = Proj(init=self.crs)

            ulx,uly = transform(inProj,outProj,west,north)

        else:
            ulx,uly = west,north

        if type(gridSize) == list:
            pxSize = gridSize[0]
            pySize = gridSize[1]
        else:
            pxSize,pySize = gridSize, gridSize

        gt = (ulx,pxSize,0,uly,0,-pySize)

        return gt

    @classmethod
    def select(self,bandList,newNames=None):
        if (type(bandList) != list) and (type(bandList) == str):
            bandList = [bandList]

        bandNames = self.raster.coords['band'].values

        newBands = [new for new in bandList if any(b in new for b in bandNames)]

        out = self.copy()

        out.raster = out.raster.sel(band=newBands)

        if newNames:
            out.raster.coords['band'] = newNames

        return out

    @classmethod
    def normalizedDifference(cls,band1=None,band2=None,appendTo=None,outBandName='nd'):
        if (band1 == None) or (band2 == None):
            raise ValueError('Band1 and Band2 need be defined to calculate normalizedDifference')

        nd = (cls.sel(band=band1) - cls.sel(band=band2)) / \
             (cls.sel(band=band1)+ cls.sel(band=band2))

        nd = nd.expand_dims('band',2)
        nd.coords['band'] = [bandName]

        if appendTo != None:
            result= xr.merge([appendTo, nd],dim='band')

        else:
            result = nd

        return result

    def gridRaster(self,):

        return

    def clip(self,geom):

        return

    @classmethod
    def updateMask(cls,arr):
        if len(arr.shape) < 3:
            if len(arr.shape) == 2:
                arr = arr[:,:,np.newaxis]
            else:
                raise ValueError('Update to mask has to be at least 2-dimensional\
                                  with (y,x) coordinates consitent with DataArray')

        cls.sel(band='mask').values = cls.sel(band='mask') & arr
        out = cls.mask()

        return out

    @classmethod
    def mask(self):
        out=self.copy()
        out.raster = out.raster.where(out.raster.sel(band='mask')==1)
        return out


    @classmethod
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

    @staticmethod
    def writeGeotiff(rasterobj,fileName):
        drv = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(rasterobj.raster.crs['init']))

        t = self.date

        outDs = drv.Create(prefix + '_' + t + '.tif',
                           y,x,bands,
                           gdal.GDT_Int16
                           )

        for b in range(bands):
            band = outDs.GetRasterBand(b+1)
            band.WriteArray(tValues[b,:,:])
            band.SetNoDataValue(0)
            band = None

        outDs.SetGeoTransform(self.data.attrs['gt'])

        outDs.SetProjection(srs.ExportToWkt())

        outDs.FlushCache()

        return


# @xr.register_dataarray_accessor('collection')
class Collection(Raster):
    def __init__(self):

        return

    def _type(self):
        return 'collection'

    def copy(self):
        return copy.deepcopy(self)


    def writeNetCDF(self,fileName):
        self.data.to_netcdf(filname)
        return


    def writeGeotiffs(self,folder,prefix='bump_out'):

        nFiles,y,x = self.data.dims.values()
        dates = self.data.time.values

        drv = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(gr.crs.split(':')[1]))


        for i in range(nFiles):
            t = str(dates[i]).split('.')[0].replace(':','')
            print(t)
            tValues = self.data.isel(time=i).to_array().values

            bands = tValues.shape[0]

            outDs = drv.Create(prefix + '_' + t + '.tif',
                               y,x,bands,
                               gdal.GDT_Int16
                               )

            for b in range(bands):
                band = outDs.GetRasterBand(b+1)
                band.WriteArray(tValues[b,:,:])
                band.SetNoDataValue(0)
                band = None

            outDs.SetGeoTransform(self.data.attrs['gt'])

            outDs.SetProjection(srs.ExportToWkt())

            outDs.FlushCache()

        return


    def filterDate(self,iniTime,endTime):
        if type(iniTime) != datetime.datetime:
            iniTime = datetime.datetime.strptime(iniTime,'%Y-%m-%d')
        if type(endTime) != datetime.datetime:
            endTime = datetime.datetime.strptime(endTime,'%Y-%m-%d')

        out = self.copy()

        out.data = out.data.sel(time=slice(iniTime, endTime))

        return out

    def filterBounds(self, geom):

        return

    def select(self,bands):
        out = self.copy()

        out.data = out.data.sel(time=slice(iniTime, endTime))

        return out

    def appendRaster(self,raster):

        return

# @xr.register_dataset_accessor('mc')
class MultiCollection(Raster):
    def __init__(self,):
        return


def save(obj,filename):
    if filename[-5:] != '.rasm':
        fileName = fileName.split('.')[0]+'.rasm'

    pickle.dump(obj, open( filename, "wb" ) )

    return

def load(filename):
    obj = pickle.load( open( filename, "rb" ) )

    return obj
