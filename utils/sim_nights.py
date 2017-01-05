import sys
#from ExposureCalc import *

import numpy as np
import pickle
import ephem
import optparse
from datetime import datetime, timedelta
import random
import re
import os
import shutil

import NightSim as ns
import UCSCScheduler_V3 as ds
import ExposureCalculations as ec
import Generate_Errors as ge

def compute_simulation(curtime,result,star,apf_obs,slowdowns,fwhms,star_tab):
    actel,actaz = ns.compute_el(curtime,star,apf_obs)
    actslow, actfwhm = ns.rand_obs_sample(slowdowns,fwhms)
    if actslow < 0.3:
        actslow = 0.3
    actfwhm = ns.gen_seeing_el(actfwhm,actel)
    lastfwhm = actfwhm
    lastslow = actslow
    metersig = np.random.randn(1)
    specsig = np.random.randn(1)
    if abs(specsig) > 3:
        specsig = 3.
    if abs(metersig) > 3:
        metersig = 3.
    
    meterrate = ec.getEXPMeter_Rate(result['VMAG'],result['BV'],actel,actfwhm)
    meterrate *= 1 + 0.11*metersig
    meterrate /= actslow
    specrate = ec.getSpec_Rate(result['VMAG'],result['BV'],actel,actfwhm)
    specrate *= 1 + 0.11*specsig
    specrate /= actslow
    metertime = result['COUNTS'] / meterrate
    exp_time = result['EXP_TIME']
    barycentertime = curtime
    if metertime < exp_time:
        fexptime = metertime
    else:
        fexptime = exp_time
        
    curtime += (fexptime+40.)/86400
    barycentertime += fexptime/(2.*86400)
    totcounts = fexptime * specrate

#    precision, deviation, true_error = ge.compute_real_uncertainty(totcounts,result['BV'],jitter)
    precision, deviation, true_error = ge.compute_real_uncertainty_sinnoise(totcounts,ephem.julian_date(curtime),star_tab)
    
    outstr = "%s %s %.5f %.1f %.1f %.2f %.2f %.2f %.2f %.2f" %(result['NAME'] , ephem.Date(curtime), ephem.julian_date(ephem.Date(barycentertime)), fexptime, totcounts, precision, deviation, true_error, actfwhm, actslow)
#    print outstr
#    for outfp in outfps:
#        outfp.write(outstr + "\n")
    return curtime, lastfwhm, lastslow, outstr


def read_datefile(datefn):

    try:
        datefile = open(options.datefile)
    except Exception as e:
        print "cannot open %s, and we died trying: %s" % (args[0], e)
        sys.exit()

    datelist = []
    for line in datefile:
        datestr, = line.split()
        if not ns.checkdate(datestr):
            print "%s is not an acceptable date string" % (datestr)
            sys.exit()
        datelist.append(datestr)
        
    return datelist

def gen_datelist(startstr,endstr,double=False):
    datelist = []
    start = datetime.strptime(startstr,"%Y/%m/%d")
    end  = datetime.strptime(endstr,"%Y/%m/%d")

    rng = [2,3,4]

    cur = start
    while cur < end:
        breakbeg = datetime(cur.year,12,24)
        breakend = datetime(cur.year+1,1,1)
        td = timedelta(random.choice(rng))
        cur = cur + td
        if cur < end and (cur < breakbeg or cur > breakend):
            datestr = "%d/%02d/%02d" % (cur.year,cur.month,cur.day)
            datelist.append(datestr)
        if double:
            cur = cur + timedelta(1)
            if cur < end and (cur < breakbeg or cur > breakend):
                datestr = "%d/%02d/%02d" % (cur.year,cur.month,cur.day)
                datelist.append(datestr)
                        
    return datelist


def write_sim_file_results(star_strs,outdir="../SimFiles"):
    outpath = os.path.join(os.curdir,outdir)
    for starname in star_strs.keys():
        ofn = "%s.sim" % (starname)
        ofn = os.path.join(outpath,ofn)
        if os.path.exists(ofn):
            ofp = open(ofn,"a+")
        else:
            ofp = open(ofn,"w")
        for outstr in star_strs[starname]:
            (name,date,time,mjd,expt,i2cnts,prec,dev,tru_err,fwhm,slow) = outstr.split()
            outstr = "%s %s %s %s\n" % (mjd,i2cnts,prec,dev)
            ofp.write(outstr)
        ofp.close()
    

def sim_results(outstr,star_strs,star_dates):
    (name,date,time,mjd,expt,i2cnts,prec,dev,tru_err,fwhm,slow) = outstr.split()
    try:
        mjd = float(mjd)
        if name in star_strs.keys():
            star_strs[name].append(outstr)
            star_dates[name].append( mjd )
        else:
            star_strs[name] = [outstr]
            star_dates[name] = [ mjd ]
    except:
        pass
    return


def prep_master(outdir,mastername):
    star_strs = dict()
    star_dates = dict()
    mastername = os.path.join(outdir,mastername)
    newmaster = False
    if not os.path.exists(mastername):
        newmaster = True
    else:
        try:
            masterfp = open(mastername)
        except Exception as e:
            print "Cannot open file %s for output, %s,  exiting" % (mastername, e)
            sys.exit()
        
        for ln in masterfp:
            sim_results(ln,star_strs,star_dates)
        masterfp.close()
    
    try:
        masterfp = open(mastername,"a+")
    except Exception as e:
        print "Cannot open file %s for output, %s,  exiting" % (mastername, e)
        sys.exit()

    return masterfp,star_strs, star_dates


def parse_args():
    parser = optparse.OptionParser()
    parser.add_option("-g","--googledex",dest="googledex",default="The Googledex")
    parser.add_option("-i","--infile",dest="infile",default="newgoogledex.csv")
    parser.add_option("-f","--file",dest="datefile",default="")
    parser.add_option("-s","--seed",dest="seed",default=None)
    parser.add_option("-o","--outdir",dest="outdir",default=".")        
    parser.add_option("-d","--double",dest="double",default=False,action="store_true")
    parser.add_option("-m","--masterfile",dest="master",default="sim_master.simout")
    parser.add_option("-p","--priority",dest="method",default="inquad",choices=["inquad","outquad","uniform","random"])
    (options, args) = parser.parse_args()    

    if len(args) < 2 and options.datefile == "":
        print "needs either a date pair or an input file which lists the dates"
        sys.exit()    

    if options.datefile != "":
        df = os.path.join(options.outdir,options.datefile)
        datelist = read_datefile(df)
    else:
        datelist = gen_datelist(args[0],args[1],double=options.double)

    if options.seed:
        random.seed(int(options.seed))
        np.random.seed(int(options.seed))

    if not os.path.isdir(options.outdir):
        try:
            os.mkdir(options.outdir)
        except Exception as e:
            print "cannot make output directory: %s - %s" % (options.outdir,e)
            sys.exit()

    gd = os.path.join(options.outdir,options.infile)
    if not os.path.exists(gd):
        try:    
            shutil.copyfile(options.infile, gd)
        except Exception as e:
            print "cannot copy %s to %s: %s" % (options.infile,options.outdir,e)
            sys.exit()
            
    return options, datelist

###

if __name__ == "__main__":


    options,datelist = parse_args()
    
    masterfp,star_strs, star_dates = prep_master(options.outdir,options.master)
    
    for datestr in datelist:

        allnames, star_table, do_flag, stars  = ds.parseGoogledex(sheetn=options.googledex,outfn=os.path.join(options.outdir,options.infile))
    
        fwhms = ns.gen_seeing()
        slowdowns = ns.gen_clouds()

        lastslow = 5
        lastfwhm = 15
        otfn = os.path.join(options.outdir,"observed_targets")
        ot = open(otfn,"w")
        ot.close()
        observing = True
        curtime, endtime, apf_obs = ns.sun_times(datestr)
        bstar = True
        standardstar=False
        while observing:

            result = ds.getNext(curtime, lastfwhm, lastslow, star_dates, bstar=bstar, standardstar=standardstar, verbose=True,googledex_file=options.infile,method=options.method)
            if result:
                if standardstar:
                    standardstar = False
                if bstar:
                    bstar = False
                    standardstar = True
                
                curtime += 70./86400 # acquisition time
                idx, = np.where(star_table['starname'] == result['NAME'])
                for i in range(0,int(result['NEXP'])):
                    (curtime,lastfwhm,lastslow,outstr) = compute_simulation(curtime,result,stars[idx],apf_obs,slowdowns,fwhms,star_table[idx])
                    sim_results(outstr,star_strs,star_dates)
                    
                ot = open(otfn,"a+")
                ot.write("%s\n" % (result["SCRIPTOBS"]))
                ot.close()
                print outstr
                masterfp.write("%s\n" % (outstr))
            else:
                curtime += 2100./86400 # close for lack of target
                lastslow = 5
                lastfwhm = 15
            if curtime > endtime:
                observing = False
        
            curtime = ephem.Date(curtime)
        
        print "sun rose"

        if os.path.isfile(otfn):
            try:
                os.unlink(otfn)
            except:
                print "cannot unlink %s" %(otfn)
    if masterfp:
        masterfp.close()
    simoutdir = os.path.join("..","SimFiles",options.outdir)
    if not os.path.isdir(simoutdir):
        try:
            os.mkdir(simoutdir)
        except Exception as e:
            print "cannot make output directory: %s - %s" % (simoutdir,e)
            sys.exit()
        
    write_sim_file_results(star_strs,outdir=simoutdir)
# done
