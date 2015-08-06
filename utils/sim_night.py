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

def sun_times(datestr):
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = -9.0*np.pi / 180.0
    apf_obs.date = datestr
    sunset = apf_obs.next_setting(ephem.Sun())
    sunset = apf_obs.next_rising(ephem.Sun())    
    return sunset,sunrise, apf_obs
    
def make_obs_sample(fn):
    slow,fwhm = np.loadtxt(fn,unpack=True)
    return slow, fwhm

def rand_obs_sample(slow,fwhm):
    ls = len(slow) -1
    lf = len(fwhm) -1
    sindx = np.random.randint(0,ls)
    findx = np.random.randint(0,lf)
    return slow[sindx], fwhm[findx]

def compute_el(curtime,star,apf_obs):
    apf_obs.date = curtime
    star.compute(apf_obs)
    actel = np.degrees(star.alt)
    return actel

datestr = "2015/08/07"
allnames, star_table, do_flag, stars  = ds.parseGoogledex()
slowdowns, fwhms = make_obs_sample("slowdowns")
lastslow = 5
lastfwhm = 10
otfn = "observed_targets"
ot = open(otfn,"w")
observing = True
curtime, endtime, apf_obs = sun_times(datestr)

while observing:

    result = getNext(curtime, lastfwhm, lastslow, bstar=True, verbose=True)
    idx = allnames.index(result['NAME'])
    for i in range(0,result['COUNTS']):
        actslow, actfwhm = rand_obs_sample(slowdowns,fwhms)
        actel = compute_el(curtime,stars[idx],apf_obs)
        
        meterrate = getEXPMeter_Rate(result['VMAG'],result['BV'],actel,actfwhm)
        meterrate *= 1 + 0.11*np.random.randn(1)
        meterrate *= actslow
        exp_time = getEXPTime(results['I2CNTS'],result['VMAG'],result['BV'],actel,actfwhm)
        exp_time *= 1 + 0.11*np.random.randn(1)
        exp_time *= actslow
        metertime = res['COUNTS'] / meterrate
        if metertime < exp_time:
            curtime += metertime/86400
        else:
            curtime += exp_time/86400
        curtime = ephem.Date(curtime)
        if currtime > endtime:
            observing = False
        
    ot.write("%s\n" % (result["SCRIPTOBS"]))
    ot.close()
