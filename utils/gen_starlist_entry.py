#!/usr/bin/env  /opt/kroot/bin/kpython
import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
from optparse import OptionParser
from datetime import datetime
import time

if __name__ == "__main__":

    THRESHOLD = 60 * 60
    parser = OptionParser()
    parser.add_option("-s","--slowdown",dest="slowdown",default=0.4,type="float")
    parser.add_option("-f","--fwhm",dest="fwhm",default=13,type="float")
    parser.add_option("-e","--el",dest="el",default=70,type="float")
    (options, args) = parser.parse_args()    

    allnames, star_table, flags, stars  = ds.parseGoogledex()

    el = np.zeros_like(star_table[:, ds.DS_BV])
    el += options.el
    fwhm = np.array(options.fwhm)
    i2counts = np.zeros_like(star_table[:, ds.DS_BV])
    totexptimes = np.zeros_like(star_table[:, ds.DS_BV])

    precision = star_table[:, ds.DS_ERR]
    i2counts = ds.getI2_K(precision)
    mstars_inds = np.where(star_table[:, ds.DS_BV] > 1.2)
    i2counts[mstars_inds] = ds.getI2_M(precision[mstars_inds])

    exp_times, exp_counts, i2cnts = ds.calculate_ucsc_exposure_time(star_table[:, ds.DS_VMAG],precision,el,fwhm,star_table[:, ds.DS_BV])
    exp_times *= options.slowdown
    totexptimes += ds.computeMaxTimes(exp_times,star_table[:, ds.DS_MAX])

    mxtime = np.zeros_like(star_table[:,ds.DS_MAX])
    mxtime += ds.MAX_EXPTIME
    shorter = (star_table[:,ds.DS_MAX] < ds.MAX_EXPTIME) & (star_table[:,ds.DS_MAX] >0)
    mxtime[shorter] = star_table[:,ds.DS_MAX][shorter]

    star_table[:, ds.DS_EXPT], exps = ds.format_time(totexptimes,i2counts,star_table[:, ds.DS_NSHOTS],star_table[:, ds.DS_MIN],mxtime)
    star_table[:, ds.DS_COUNTS], star_table[:, ds.DS_NSHOTS] = ds.format_expmeter(exp_counts,exps)

    fin_pre = precision
   
    dt = datetime.utcfromtimestamp(int(time.time()))
    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] < 5:
            continue
        row = []

        row.append( star_table[i, ds.DS_RA])
        row.append( star_table[i, ds.DS_DEC])
        row.append( star_table[i, ds.DS_PMRA])
        row.append( star_table[i, ds.DS_PMDEC])
        row.append( star_table[i, ds.DS_VMAG])
        row.append( star_table[i, ds.DS_EXPT])
        row.append( star_table[i, ds.DS_COUNTS])
        row.append( star_table[i, ds.DS_APFPRI])
        row.append(0)
        row.append( star_table[i, ds.DS_NSHOTS])
        line = ds.makeScriptobsLine(allnames[i],row,flags['do'][i],dt,decker=flags['decker'][i],I2=flags['I2'][i])
        print "%s #  %s" % (line,"pri = %s" % (star_table[i, ds.DS_APFPRI]))

