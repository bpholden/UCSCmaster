import sys
import optparse
from datetime import datetime
import re
import os
import pytz

import numpy as np
import pickle
import astropy
import astropy.coordinates
import astropy.time
import astropy.units
import astroplan




def makeAPFObs():
    # Generate a astropy.coordinate observer for the APF

    apf_lat = astropy.coordinates.Latitude((37,20,33.1),unit=astropy.units.deg)
    apf_long = astropy.coordinates.Longitude((-121,38,17.7),wrap_angle=astropy.units.deg*180,unit=astropy.units.deg)
    apf_height = 1274 * astropy.units.meter
    apf_loc = astropy.coordinates.EarthLocation.from_geodetic(apf_long,apf_lat,apf_height)

    apf_obs = astroplan.Observer(name='APF Telescope',
               location=apf_loc,
               pressure=870 * astropy.units.hPa,
               relative_humidity=0.3,
               temperature=20 * astropy.units.deg_C,
               timezone=pytz.timezone('US/Pacific'),
               description="APF Telescope on Mount Hamilton, California")
    
    return apf_obs

def sun_times(year,month,day):

    apf_obs = makeAPFObs()

    horizon_deg = -9 * astropy.units.deg
    
    compute_time = astropy.time.Time(datetime(year=year,month=month,day=day,hour=20))
    
    sunset = apf_obs.sun_set_time(compute_time,which='next',horizon=horizon_deg)
    sunrise = apf_obs.sun_rise_time(compute_time,which='next',horizon=horizon_deg)
    return sunset, sunrise, apf_obs
    
def make_obs_sample(fn):
    slow,fwhm = np.loadtxt(fn,unpack=True)
    return slow, fwhm

def gen_seeing(nsize=200,val=-1):
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.52
    
    if val < 0.9:
        mean = np.random.normal(loc=8.,scale=1.0,size=1)
        rms  = np.random.normal(loc=1.9,scale=1.0,size=1)
    else:
        mean = np.random.normal(loc=19.0,scale=1.0,size=1)
        rms  = np.random.normal(loc=4.5,scale=1.0,size=1)
    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean

#    deviates = np.random.normal(loc=mean,scale=rms,size=nsize)
    return deviates
        
def gen_seeing_el(deviate,el):
    zd = 90 - el
    deviate += (0.0903544076597*zd +  -0.00172591889888*zd*zd + 3.3157238117e-05*zd*zd*zd)
    return deviate

def gen_clouds(nsize=200,val=-1):
    if val < 0:
        val = np.random.uniform(size=1)
    alpha = 0.353
    if val < 0.7:
        mean = np.random.normal(loc=2.,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.05,scale=0.01,size=1)
    elif val > 0.9:
        mean = np.random.normal(loc=4.0,scale=0.4,size=1)
        rms  = np.random.normal(loc=0.5,scale=0.1,size=1)
    else:
        mean = np.random.normal(loc=3.4,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.3,scale=0.05,size=1)

    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean
    deviates[deviates < 0.3] = 0.3
    return deviates

def rand_obs_sample(slows,fwhms):
    ls = len(slows) -1
    lf = len(fwhms) -1
    sindx = np.random.randint(0,ls)
    findx = np.random.randint(0,lf)
    return slows[sindx], fwhms[findx]

def compute_el(curtime,star,apf_obs):
    altaz = apf_obs.altaz(curtime,target=star)
    actel = altaz.alt.value
    actaz = altaz.az.value
    return actel,actaz


def checkdate(datestr):
    match = re.match("(\d{4})\-(\d{1,2})\-(\d{1,2})",datestr)
    if not match:
        return False, -1, 0, 0
    if int(match.group(2)) < 1 or int(match.group(2)) > 12:
        return False, -1, 0, 0
    if int(match.group(3)) < 1 or int(match.group(3)) > 31:
        return False, -1, 0, 0
    
    return True, int(match.group(1)), int(match.group(2)), int(match.group(3))


