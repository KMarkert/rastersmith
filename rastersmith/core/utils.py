
# import dependency
from __future__ import print_function,division

import math
import warnings
from itertools import groupby
import numpy as np
import xarray as xr
from scipy import interpolate


# warnings.simplefilter("always")

def meters2dd(inPt,scale=30):
    """Function to convert meters to decimal degrees based on the approximation
    given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system

    Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.

    Keywords:
        scale (int): Resolution of the raster value to covert into decimal
        degrees, must be in meters.

    Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees

    """

    lat = inPt[0] # get latitude value

    radLat = math.radians(lat) # convert degree latitude to radians

    a = 6378137 # radius of Earth in meters

    ba = 0.99664719 # constant of b/a

    ss = math.atan(ba*math.tan(radLat)) # calculate the reduced latitude

    # factor to convert meters to decimal degrees for X axis
    xfct = (math.pi/180)*a*math.cos(ss)

    # factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*math.cos(2*radLat)+1.175*math.cos(4*radLat)-
              0.0023*math.cos(6*radLat))

    # get decimal degree resolution
    ydd = scale / yfct
    xdd = scale / xfct

    # return list of converted resolution values
    return ydd

def dd2meters(inPt,scale=0.1):
    """Function to convert decimal degrees to meters based on the approximation
        given by: https://en.wikipedia.org/wiki/Geographic_coordinate_system

        Args:
        inPt (list or array): A Y,X point provided in geographic coordinates
        in that order.

        Keywords:
        scale (int): Resolution of the raster value to covert into meters,
        must be in decimal degrees.

        Returns:
        list: List of Y,X resolution values converted from meters to decimal degrees

        """

    lat = inPt[0] # get latitude value

    radLat = math.radians(lat) # convert degree latitude to radians

    a = 6378137 # radius of Earth in meters

    ba = 0.99664719 # constant of b/a

    ss = math.atan(ba*math.tan(radLat)) # calculate the reduced latitude

    # factor to convert meters to decimal degrees for X axis
    xfct = (math.pi/180)*a*math.cos(ss)

    # factor to convert meters to decimal degrees for Y axis
    yfct = (111132.92-559.82*math.cos(2*radLat)+1.175*math.cos(4*radLat)-
            0.0023*math.cos(6*radLat))

    # get meter resolution
    y_meters = scale * yfct
    x_meters = scale * xfct

    # return list of converted resolution values
    return y_meters

def rasterExpression(expression,lookups,outBandName='bandMath',appendTo=None):

    result = eval(expression,lookups)

    result = result.expand_dims('band')
    result.coords['band'] = [outBandName]
    result = result.transpose('lat','lon','z','band','time')

    if appendTo:
        out = xr.concat([appendTo,result],dim='band')
    else:
        out = result

    return out


def find_nearest_idx(pt,xx,yy):
    xval, yval = pt

    # Find closest image index to the x-y coordinate
    xidx = (np.abs(xx-xval)).argmin()
    yidx = (np.abs(yy-yval)).argmin()

    # Convert the 1-d index to 2-d
    ridx = yidx / xx.shape[1]
    cidx = xidx % xx.shape[1]

    return (cidx,ridx)

def formatDataarr(dataarr):
    result = np.moveaxis(np.array(dataarr),0,2)[:,:,:,np.newaxis]
    result = np.moveaxis(result,2,3)[:,:,:,:,np.newaxis]

    return result
