#!/usr/bin/env  /opt/kroot/bin/kpython
import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
from optparse import OptionParser

if __name__ == "__main__":

    #THRESHOLD = 2700.0
    THRESHOLD = 60 * 60
    #SLOWDOWN = 0.4
    #STAR_EL = 70
    #AVG_FWHM = 11
    parser = OptionParser()
    parser.add_option("-s","--slowdown",dest="slowdown",default=0.4,type="float")
    parser.add_option("-f","--fwhm",dest="fwhm",default=13,type="float")
    parser.add_option("-e","--el",dest="el",default=70,type="float")
    (options, args) = parser.parse_args()    

#    ws = ds.get_speadsheet(sheetn="The Googledex")
#    vals = ws.get_all_values()
#    texpcol = vals[0].index("APFtexp") 
    
    allnames, star_table, do_flag, stars  = ds.parseGoogledex()

    el = np.zeros_like(star_table[:, ds.DS_BV])
    el += options.el
    fwhm = np.array(options.fwhm)
    i2counts = np.zeros_like(star_table[:, ds.DS_BV])

    precision = star_table[:, ds.DS_ERR]
    i2counts = ds.getI2_K(precision)
    mstars_inds = np.where(star_table[:, ds.DS_BV] > 1.2)
    i2counts[mstars_inds] = ds.getI2_M(precision[mstars_inds])
#        exp_counts = ds.getEXPMeter(i2counts, star_table[:, ds.DS_BV])
#    exp_time = ds.getEXPTime(i2counts, star_table[:, ds.DS_VMAG], star_table[:, ds.DS_BV], el, fwhm)

    exp_times, exp_counts, i2cnts = ds.calculate_ucsc_exposure_time(star_table[:, ds.DS_VMAG],precision,el,fwhm,star_table[:, ds.DS_BV])
    exp_times *= options.slowdown
    etimes, nobs = ds.format_time(exp_times,i2cnts)
    exp_counts, nobs = ds.format_expmeter(exp_counts,nobs)
    fin_pre = precision
#    exp_counts /= nobs
    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] < 5:
            continue

        print "%15s %4.1f %3.1f %7.0f %7.0f %.3g %.1f %d" % (allnames[i],star_table[i, ds.DS_APFPRI],precision[i],i2counts[i],exp_times[i],exp_counts[i],etimes[i],nobs[i])

    # plt.scatter(times, err, c=pri, cmap=plt.get_cmap("jet"), edgecolor='none')
    # plt.xlabel("Exposure Time")
    # plt.ylabel("Desired Precision")
    # plt.show()
