#!/usr/bin/env  /opt/kroot/bin/kpython
import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
from optparse import OptionParser
import datetime

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-s","--slowdown",dest="slowdown",default=0.4,type="float")
    parser.add_option("-f","--fwhm",dest="fwhm",default=13,type="float")
    parser.add_option("-e","--el",dest="el",default=70,type="float")
    (options, args) = parser.parse_args()    

    desiredstars = []
    if args:
        desiredstars = args
#    ws = ds.get_speadsheet(sheetn="The Googledex")
#    vals = ws.get_all_values()
#    texpcol = vals[0].index("APFtexp") 
    
    allnames, star_table, flags, stars  = ds.parseGoogledex()
    if len(desiredstars) == 0:
        desiredstars = allnames
    el = np.zeros_like(star_table[:, ds.DS_BV])
    el += options.el
    fwhm = np.array(options.fwhm)
    i2counts = np.zeros_like(star_table[:, ds.DS_BV])

    precision = star_table[:, ds.DS_ERR]
    i2counts = ds.getI2_K(precision)
    mstars_inds = np.where(star_table[:, ds.DS_BV] > 1.2)
    i2counts[mstars_inds] = ds.getI2_M(precision[mstars_inds])

    exp_times, exp_counts, i2cnts = ds.calculate_ucsc_exposure_time(star_table[:, ds.DS_VMAG],precision,el,fwhm,star_table[:, ds.DS_BV])
    exp_times *= options.slowdown
    etimes, nobs = ds.format_time(exp_times,i2cnts,star_table[:,ds.DS_NSHOTS])
    exp_counts, nobs = ds.format_expmeter(exp_counts,nobs)
    fin_pre = precision
    for star in desiredstars:
        i = allnames.index(star)
        if star_table[i, ds.DS_APFPRI] < 5:
            continue
#        print "%15s %4.1f %3.1f %7.0f %7.0f %.3g %.1f %d" % (allnames[i],star_table[i, ds.DS_APFPRI],precision[i],i2counts[i],exp_times[i],exp_counts[i],etimes[i],nobs[i])
        ret = ds.makeScriptobsLine(allnames[i],star_table[i,:], flags['do'][i], datetime.datetime.utcnow(),decker=flags['decker'][i])
        print ret
