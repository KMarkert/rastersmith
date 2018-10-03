from __future__ import print_function, division

import numpy as np
import bottleneck as bn
import xarray as xr
from scipy import interpolate
from pyproj import Proj, transform
from PIL import Image, ImageDraw
import multiprocessing as mp

from ..core import utils
from ..core import core


def _griddingWrapper(args):
    pts,rasterArr,xx,yy,method,i = args
    return np.flipud(interpolate.griddata(pts,rasterArr[:,:,0,i,0].ravel(),(xx,yy),method=method))

# @staticmethod
def swath2grid(rasterArr,xgrid,ygrid,resolution=500,method='nearest'):
    if method not in ['nearest','linear','cubic']:
        raise ValueError('A valid interpolation method was not defined. Options are "nearest", "linear", or "cubic"')

    lats,lons = ygrid,xgrid

    wVerts = [(lons[i,0],lats[i,0]) for i in range(lats.shape[0])][::-1]
    nVerts = [(lons[0,i],lats[0,i]) for i in range(lats.shape[1])]
    eVerts = [(lons[i,-1],lats[i,-1]) for i in range(lats.shape[0])]
    sVerts = [(lons[-1,i],lats[-1,i]) for i in range(lats.shape[1])]
    sVerts.append(wVerts[0])

    verts = wVerts + nVerts + eVerts + sVerts

    xDim = np.arange(bn.nanmin(lons),bn.nanmax(lons),resolution)
    yDim = np.arange(bn.nanmin(lats),bn.nanmax(lats),resolution)

    xx,yy = np.meshgrid(xDim,yDim)

    imgVerts = list(map(lambda x: utils.find_nearest_idx(x,xx,yy),verts))

    # Create a mask from the image polygon vertices
    img = Image.new('L', (xx.shape[1], xx.shape[0]), 0)
    ImageDraw.Draw(img).polygon(imgVerts, fill=1)
    mask = np.flipud(np.array(img))

    # Format geolocation coordinates for gridding
    pts = np.zeros((lons.size,2))
    pts[:,0] = lons.ravel()
    pts[:,1] = lats.ravel()

    args = [[pts,rasterArr,xx,yy,method,i] for i in range(rasterArr.shape[3])]
    args[-1][4] = 'nearest'

    nProcessors = mp.cpu_count()
    pool = mp.Pool(processes=nProcessors)
    dataarr = pool.map(_griddingWrapper,args)

    result = utils.formatDataarr(dataarr)

    return result,mask,xDim,yDim[::-1]


def coregister(raster,to=None,method='nearest'):

    if type(to) == core.Grid:
        xx = to.xx
        yy = to.yy

    elif type(to) == xr.core.dataarray.DataArray:
        xCoords = to.coords['lon'].values
        yCoords = to.coords['lat'].values

        xx, yy = np.meshgrid(xCoords,yCoords)

    else:
        raise NotImplementedError('Could not interpret Variable passed to argument to')

    projX = xr.DataArray(xx,coords={'y':range(xx.shape[0]),'x':range(xx.shape[1])},dims=('y','x'))
    projY = xr.DataArray(yy,coords={'y':range(xx.shape[0]),'x':range(xx.shape[1])},dims=('y','x'))

    rasterSel = raster.sel(dict(lat=slice(yy.max(),yy.min()),lon=slice(xx.min(),xx.max())))

    outDa = rasterSel.interp(lat=projY,lon=projX,method=method)

    outDa = outDa.drop(['lat','lon'])
    outDa['x'],outDa['y'] = xx[0,:],yy[:,0]

    return outDa.rename({'x':'lon','y':'lat'})


def reproject(raster,outEpsg='4326',outResolution=500,method='nearest'):
    xCoords = raster.coords['lon'].values
    yCoords = raster.coords['lat'].values

    xGrid, yGrid = np.meshgrid(xCoords,yCoords)

    if '4326' in str(outEpsg):
        inproj = Proj(raster.attrs['projStr'])
        outproj = Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
    else:
        raise NotImplementedError('Curretly only reprojecting to WGS84 is support, future versions will support all reprojections')

    xx,yy = transform(inproj,outproj,xGrid,yGrid)

    if '4326' in str(outEpsg):
        outRes = utils.meters2dd((xx.mean(),yy.mean()),outResolution)
    else:
        outRes = outResolution

    newX = np.arange(xx.min(),xx.max(),outRes)
    newY = np.arange(yy.min(),yy.max(),outRes)[::-1]

    newXGrid, newYGrid = np.meshgrid(newX,newY)

    projX,projY = transform(outproj,inproj,newXGrid,newYGrid)

    projX = xr.DataArray(projX,coords={'y':range(newXGrid.shape[0]),'x':range(newXGrid.shape[1])},dims=('y','x'))
    projY = xr.DataArray(projY,coords={'y':range(newXGrid.shape[0]),'x':range(newXGrid.shape[1])},dims=('y','x'))

    outDa = raster.interp(lat=projY,lon=projX,method=method)

    outDa = outDa.drop(['lat','lon'])
    outDa['x'],outDa['y'] = newX,newY

    return outDa.rename({'x':'lon','y':'lat'})

# @staticmethod
def reduceNeighborhood(raster,window=3,reducer='mean',func=None,reduceResolution=False):
    method = {'mean':bn.nanmean,'std':bn.nanstd,'variance':np.var,'sum':np.sum,'unique':np.unique}

    if reducer not in list(method.keys()):
        raise NotImplementedError("Selected reducer is not implemented")

    if func:
        resampled = raster.rolling(lat=window) \
                          .rolling(lon=window)

    else:
        resampled = raster.rolling(lat=window).reduce(method[reducer]) \
                          .rolling(lon=window).reduce(method[reducer])

    if reduceResolution:
        start = int(window/2)
        newX = resampled.coords['lon'][start::window]
        newY = resampled.coords['lat'][start::window]

        resampled = resampled.interp(lat=newY,lon=newX,method='nearest')

    resampled.attrs = raster.attrs

    return resampled
