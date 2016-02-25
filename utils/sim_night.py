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


def checkdate(datestr):
    match = re.match("(\d{4})\/(\d{1,2})\/(\d{1,2})",datestr)
    if not match:
        return False
    if int(match.group(2)) < 1 or int(match.group(2)) > 12:
        return False
    if int(match.group(3)) < 1 or int(match.group(3)) > 31:
        return False
    
    return True

#######

parser = optparse.OptionParser()
parser.add_option("-d","--date",dest="date",default="today")
parser.add_option("-f","--fixed",dest="fixed",default="")
parser.add_option("-s","--smartlist",dest="smartlist",default=False,action="store_true")
parser.add_option("-g","--googledex",dest="googledex",default="The Googledex")
parser.add_option("-i","--infile",dest="infile",default="googledex.dat")
parser.add_option("-b","--bstar",dest="bstar",default=True,action="store_false")
(options, args) = parser.parse_args()    

if options.date == "today":
    today = datetime.now()
    datestr = "%d/%02d/%02d" % (today.year,today.month,today.day)
else:
    datestr = options.date

if options.fixed != "":
    if not os.path.isfile(options.fixed):
        print "%s is not a file" % (options.fixed)

if not checkdate(datestr):
    print "%s is not an acceptable date string" % (datestr)
    sys.exit()
    
if options.fixed != "":
    allnames, star_table, lines, stars = ds.parseStarlist(options.fixed)
else:
    allnames, star_table, do_flag, stars  = ds.parseGoogledex(sheetn=options.googledex,outfn=options.infile)

slowdowns, fwhms = make_obs_sample("slowdowns")
lastslow = 5
lastfwhm = 15
otfn = "observed_targets"
ot = open(otfn,"w")
ot.close()
observing = True
curtime, endtime, apf_obs = sun_times(datestr)
bstar = options.bstar
while observing:

    if options.smartlist and options.fixed != "":
        result = ds.smartList(options.fixed, curtime, lastfwhm, lastslow)
    else:
        result = ds.getNext(curtime, lastfwhm, lastslow, bstar=bstar, verbose=True)
    if result:
        if bstar:
            bstar = False
        curtime += 70./86400 # acquisition time
        idx = allnames.index(result['NAME'])
        for i in range(0,int(result['NEXP'])):
            actslow, actfwhm = rand_obs_sample(slowdowns,fwhms)
            actel = compute_el(curtime,stars[idx],apf_obs)
            lastfwhm = actfwhm
            lastslow = actslow
            meterrate = ec.getEXPMeter_Rate(result['VMAG'],result['BV'],actel,actfwhm)
            meterrate *= 1 + 0.11*np.random.randn(1)
            meterrate *= actslow
            specrate = ec.getSpec_Rate(result['VMAG'],result['BV'],actel,actfwhm)
            specrate *= 1 + 0.11*np.random.randn(1)
            specrate *= actslow
            metertime = result['COUNTS'] / meterrate
            exp_time = result['EXP_TIME']
            if metertime < exp_time:
                curtime += (metertime+40.)/86400
                totcounts = metertime * specrate
            else:
                curtime += (exp_time+40.)/86400
                totcounts = exp_time * specrate
            print "%s %s %.1f %.1f %.1f\n" %(result['NAME'] , ephem.Date(curtime), exp_time, metertime, totcounts)
        ot = open(otfn,"a+")
        ot.write("%s\n" % (result["SCRIPTOBS"]))
        ot.close()
    else:
        curtime += 2100./86400 # close for lack of target
        lastslow = 5
        lastfwhm = 10
    if curtime > endtime:
        observing = False
        
    curtime = ephem.Date(curtime)
        
print "sun rose"
