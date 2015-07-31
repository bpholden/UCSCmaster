#!/usr/bin/env  /opt/kroot/bin/kpython
import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
import pickle


#THRESHOLD = 2700.0
THRESHOLD = 2* 60 * 60
SLOWDOWN = 0.4
STAR_EL = 70
AVG_FWHM = 11

        
if __name__ == "__main__":

    allnames, star_table, do_flag, stars  = ds.parseGoogledex()

    el = np.zeros_like(star_table[:, ds.DS_BV])
    el += STAR_EL
    fwhm = np.array(AVG_FWHM)
    i2counts = np.zeros_like(star_table[:, ds.DS_BV])

    precision = star_table[:, ds.DS_ERR]
    i2counts = ds.getI2_K(precision)
    mstars_inds = np.where(star_table[:, ds.DS_BV] > 1.2)
    i2counts[mstars_inds] = ds.getI2_M(precision[mstars_inds])
#        exp_counts = ds.getEXPMeter(i2counts, star_table[:, ds.DS_BV])
#    exp_time = ds.getEXPTime(i2counts, star_table[:, ds.DS_VMAG], star_table[:, ds.DS_BV], el, fwhm)

    exp_times, exp_counts, i2cnts = ds.calculate_ucsc_exposure_time(star_table[:, ds.DS_VMAG],precision,el,fwhm,star_table[:, ds.DS_BV])

    etimes, nobs = ds.format_time(exp_times,i2cnts)
    fin_pre = precision
    exp_counts /= nobs
    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] < 5:
            continue
        if star_table[i, ds.DS_APFPRI] > 9.9:
            etimes[i] = 1200.0

        print "%15s %4.1f %3.1f %7.0f %7.0f %.3g %.1f %d" % (allnames[i],star_table[i, ds.DS_APFPRI],precision[i],i2counts[i],exp_times[i],exp_counts[i],etimes[i],nobs[i])

    # plt.scatter(times, err, c=pri, cmap=plt.get_cmap("jet"), edgecolor='none')
    # plt.xlabel("Exposure Time")
    # plt.ylabel("Desired Precision")
    # plt.show()
