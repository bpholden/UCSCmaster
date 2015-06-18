import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
import pickle


#THRESHOLD = 2700.0
THRESHOLD = 2* 60 * 60
SLOWDOWN = 0.3
STAR_EL = 70
AVG_FWHM = 8

def calc_fin_pre(i2counts,exp_time,bmv):
    ratio = exp_time / THRESHOLD
    ni2counts = i2counts / ratio
    if bmv > 1.2:
        A = 4.14
        B = -1.73
    else:
        A = 4.47
        B = -1.58

    log_npre = (np.log10(ni2counts) - A )/ B
    return 10**log_npre,ni2counts

        
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

    exp_times, exp_counts = ds.calculate_ucsc_exposure_time(star_table[:, ds.DS_VMAG],precision,el,fwhm,star_table[:, ds.DS_BV])

    etimes, nobs = ds.format_time(exp_times)
    fin_pre = precision
    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] < 5:
            continue
        
#        if exp_time[i] > THRESHOLD:
#            fin_pre,ni2counts = calc_fin_pre(i2counts,exp_time,star_table[i, ds.DS_BV])
#            nexp_counts = ds.getEXPMeter(ni2counts,star_table[i, ds.DS_BV])
#        else:
#            
#            nexp_counts = exp_counts
            
        print "%15s %4.1f %3.1f %7.0f %7.0f %.1g %3.1f %.1g %.1f %d" % (allnames[i],star_table[i, ds.DS_APFPRI],precision[i],i2counts[i],exp_times[i],nobs[i],fin_pre[i],exp_counts[i],etimes[i],nobs[i])

    # plt.scatter(times, err, c=pri, cmap=plt.get_cmap("jet"), edgecolor='none')
    # plt.xlabel("Exposure Time")
    # plt.ylabel("Desired Precision")
    # plt.show()
