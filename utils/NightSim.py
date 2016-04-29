import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import ExposureCalculations as ec

import numpy as np
import pickle
import ephem
import optparse
from datetime import datetime
import re
import os

def sun_times(datestr):
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = -9.0*np.pi / 180.0
    apf_obs.date = datestr
    sunset = apf_obs.next_setting(ephem.Sun())
    sunrise = apf_obs.next_rising(ephem.Sun())    
    return sunset,sunrise, apf_obs
    
def make_obs_sample(fn):
    slow,fwhm = np.loadtxt(fn,unpack=True)
    return slow, fwhm

def gen_seeing(nsize=200):
    val = np.random.uniform(size=1)
    alpha = 0.52
    
    if val < 0.9:
        mean = np.random.normal(loc=33.46,scale=1.0,size=1)
        rms  = np.random.normal(loc=1.9,scale=1.0,size=1)
    else:
        mean = np.random.normal(loc=41.0,scale=1.0,size=1)
        rms  = np.random.normal(loc=4.5,scale=1.0,size=1)
    real_rms = np.sqrt((1-alpha**2) * rms**2)
    deviates = np.random.normal(loc=0,scale=real_rms,size=nsize)
    for i in np.arange(1,nsize):
        deviates[i] += alpha*deviates[i-1]
    deviates += mean

#    deviates = np.random.normal(loc=mean,scale=rms,size=nsize)
    return deviates
        
def gen_seeing_el(deviate,el):
    deviate += (-0.823637535775*el + 0.0104081289798*el*el +  -4.45020842236e-05*el*el*el)
    return deviate

def gen_clouds(nsize=200):
    val = np.random.uniform(size=1)
    alpha = 0.353
    if val < 0.7:
        mean = np.random.normal(loc=0.4,scale=0.1,size=1)
        rms  = np.random.normal(loc=0.1,scale=0.05,size=1)
    elif val > 0.9:
        mean = np.random.normal(loc=1.0,scale=0.4,size=1)
        rms  = np.random.normal(loc=0.5,scale=0.1,size=1)
    else:
        mean = np.random.normal(loc=0.4,scale=0.1,size=1)
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
    apf_obs.date = curtime
    star.compute(apf_obs)
    actel = np.degrees(star.alt)
    return actel


def checkdate(datestr):
    match = re.match("(\d{4})\/(\d{1,2})\/(\d{1,2})",datestr)
    if not match:
        return False
    if int(match.group(2)) < 1 or int(match.group(2)) > 12:
        return False
    if int(match.group(3)) < 1 or int(match.group(3)) > 31:
        return False
    
    return True


