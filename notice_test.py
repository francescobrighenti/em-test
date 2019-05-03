import os.path
import sys

import urllib.parse

import lxml.etree

import healpy as hp
import numpy as np
from mocpy import MOC, WCS
from math import log
import matplotlib
from matplotlib.path import Path
from matplotlib.patches import PathPatch

import astropy
from  astropy.utils.data import download_file
from astropy import units as u
from astropy.time import Time, TimeDelta,  TimezoneInfo
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle
from astropy.table import Table
from astropy.samp import SAMPIntegratedClient

from datetime import datetime

def moc_confidence_region(prob, id_object, percentage):
    
    """It returns a confidence region that encloses  a given percentage of the localization probability 
       using the Multi Order Coverage map (MOC) method. The sky localization area  (in square degrees)
       is also provided for the confidence region. 
        
    Parameters
    ----------   
    infile : `str`
        the HEALPix fits file name, the local file or the link location
    percentage : `float`
        the decimal percentage of the location probability (from 0 to 1)
    
    Returns
    -------
    A local file in which the MOC is serialized in "fits" format. The file is named :
    "moc_"+'%.1f' % percentage+'_'+skymap
    ----------------------------------------------------
    "moc_" :`str, default`
        default string as file prefix
    '%.1f' % percentage+' : `str, default`
        decimal percentage passed to the function 
    skymap : `str, default`
        string after the last slash and parsing for "." 
        from infile parameter
    ----------------------------------------------------
    area_sq2 : `float, default`
        the area of the moc confidence region in square degrees.
        """
    
    # reading skymap  
    #prob = hp.read_map(infile, verbose = False)
    npix = len(prob)
    nside = hp.npix2nside(npix) #healpix resolution
        
    cumsum = np.sort(prob)[::-1].cumsum() # sort and cumulative sum
    
    # finding the minimum credible region
    how_many_ipixs, cut_percentage = min(enumerate(cumsum),
                                             key = lambda x: abs(x[1] - percentage))
    del(cumsum)
    #print ('number of ipixs',how_many_ipixs,'; cut@', round(cut_percentage,3))
    
    indices = range(0, len(prob))
    prob_indices = np.c_[prob, indices]

    sort = prob_indices[prob_indices[:, 0].argsort()[::-1]]
    ipixs = sort[0:how_many_ipixs+1, [1]].astype(int)
  
    # from index to polar coordinates
    theta, phi = hp.pix2ang(nside, ipixs)

    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
        
    # creating an astropy.table with RA[deg] and DEC[deg]
    contour_ipix = Table([ra, dec], names = ('RA', 'DEC'), meta={'name': 'first table'})

    order = int(log(nside, 2))
    
    # moc from table
    moc = MOC.from_lonlat(contour_ipix['RA'].T* u.deg, 
                          contour_ipix['DEC'].T* u.deg, order)
    
    # getting skymap name 
    skymap = id_object.rsplit('/', 1)[-1].rsplit('.')[0]
    moc_name = "moc_"+'%.1f' % percentage+"_"+skymap
    
    # return moc region in fits format
    moc.write(moc_name, format = "fits",)

    # square degrees in a whole sphere
    from math import pi
    square_degrees_sphere = (360.0**2)/pi
 
    moc = MOC.from_fits(moc_name)

    # printing sky area at the given percentage
    area_sq2 = round( ( moc.sky_fraction * square_degrees_sphere ), 1 )
    print ('The '+str((percentage*100))+'% of '+moc_name+' is ', area_sq2,'sq. deg.')





import gcn
import healpy as hp

# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)
def process_gcn(payload, root):
    # Respond only to 'test' events.
    # VERY IMPORTANT! Replace with the following code
    # to respond to only real 'observation' events.
    # if root.attrib['role'] != 'observation':
    #    return
    if root.attrib['role'] != 'test':
        return

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Respond only to 'CBC' events. Change 'CBC' to "Burst'
    # to respond to only unmodeled burst events.
    if params['Group'] != 'CBC':
        return

    # Print all parameters.
    for key, value in params.items():
        print(key, '=', value)

    if 'skymap_fits' in params:
        # Read the HEALPix sky map and the FITS header.
        prob, header = hp.read_map(params['skymap_fits'],
                                     h=True, verbose=False)
        header = dict(header)

        # Print some values from the FITS header.
        print('Distance =', header['DISTMEAN'], '+/-', header['DISTSTD'])

        id_object = header['OBJECT'] # gracedbid

    # moc creation reading prob from healpix
    moc_confidence_region(prob, id_object, percentage=0.9)

    # query galaxy

    # filtering galaxy by distance

    # save in astropy table

   # francesco


# Listen for GCNs until the program is interrupted
# (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)

# offline
#payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#process_gcn(payload, root)

