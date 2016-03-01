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

import NightSim as ns

parser = optparse.OptionParser()
parser.add_option("-d","--date",dest="date",default="today")
parser.add_option("-f","--fixed",dest="fixed",default="")
parser.add_option("-s","--smartlist",dest="smartlist",default=False,action="store_true")
parser.add_option("-g","--googledex",dest="googledex",default="The Googledex")
parser.add_option("-i","--infile",dest="infile",default="googledex.dat")
parser.add_option("-o","--outfile",dest="outfile",default=None)
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

if not ns.checkdate(datestr):
    print "%s is not an acceptable date string" % (datestr)
    sys.exit()

if options.outfile == None:
    fdatestr = re.sub("\/","-",datestr)
    outfile = "%s.simout" % (fdatestr )
else:
    outfile = options.outfile
try:
    outfp = open(outfile,"w+")
except Exception as e:
    print "cannot open file %s for output, %s,  exiting" % (outfile,e)
    sys.exit()

hdrstr = "starname date time mjd exptime i2counts fwhm slowdown\n"
outfp.write(hdrstr)
        
if options.fixed != "":
    allnames, star_table, lines, stars = ds.parseStarlist(options.fixed)
else:
    allnames, star_table, do_flag, stars  = ds.parseGoogledex(sheetn=options.googledex,outfn=options.infile)

slowdowns, fwhms = ns.make_obs_sample("slowdowns")
fwhms = ns.gen_seeing()

lastslow = 5
lastfwhm = 15
otfn = "observed_targets"
ot = open(otfn,"w")
ot.close()
observing = True
curtime, endtime, apf_obs = ns.sun_times(datestr)
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
            actel = ns.compute_el(curtime,stars[idx],apf_obs)
            actslow, actfwhm = ns.rand_obs_sample(slowdowns,fwhms)
            actfwhm = ns.gen_seeing_el(actfwhm,actel)
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
            print "%s %s %.1f %.1f %.1f" %(result['NAME'] , ephem.Date(curtime), exp_time, metertime, totcounts)
        ot = open(otfn,"a+")
        ot.write("%s\n" % (result["SCRIPTOBS"]))
        ot.close()
    else:
        curtime += 2100./86400 # close for lack of target
        lastslow = 5
        lastfwhm = 15
    if curtime > endtime:
        observing = False
        
    curtime = ephem.Date(curtime)
        
print "sun rose"
fn = "observed_targets"
if os.path.isfile(fn):
    try:
        os.unlink(fn)
    except:
        print "cannot unlink %s" %(fn)
outfp.close()
