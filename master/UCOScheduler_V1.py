# UCOScheduler_V1.py
from __future__ import print_function
import os
import sys
import time
import re
from datetime import datetime, timedelta

import numpy as np
import ephem
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter_Rate, getEXPTime
import ParseUCOSched
import ObservedLog
import Coords
from SchedulerConsts import * # I know

try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import Visible


last_objs_attempted = []


def computePriorities(star_table,available,cur_dt,flags):
    # make this a function, have it return the current priorities, than change references to the star_table below into references to the current priority list
    if any(star_table[available, DS_DUR] > 0):
        new_pri = np.zeros_like(star_table[:, DS_APFPRI])
        new_pri[available] += star_table[available,DS_APFPRI]
        delta_pri = np.zeros_like(new_pri[available])
        timedependent, = np.where(star_table[available, DS_DUR] > 0)
        for tdinx in timedependent:
            sdt = datetime(cur_dt.year,cur_dt.month,cur_dt.day,int(star_table[:, DS_UTH][available][tdinx]),int(star_table[:, DS_UTM][available][tdinx]),0)
            durdelt = timedelta(0,star_table[:, DS_DUR][available][tdinx],0)
            if (cur_dt - sdt < durdelt) and (cur_dt - sdt > timedelta(0,0,0) ):
                delta_pri[tdinx] += PRI_DELTA
        new_pri[available] += delta_pri
    else:
        new_pri = star_table[:, DS_APFPRI]
    return new_pri



def makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y",owner='Vogt',focval=0):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y")
    name - name of star, first column in line
    row - star_table row for star that begins with name, cotains all of the data needed for the line except
    do_flag - a string for whether or not scriptob needs to do a pointing check before slewing to the target
    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    """

    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = name + ' '
    # Add the RA as three elements, HR, MIN, SEC
    rastr = Coords.getCoordStr(np.degrees(row[DS_RA]), isRA=True)
    ret += rastr + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = Coords.getCoordStr(np.degrees(row[DS_DEC]))
    ret += decstr + ' '
    # Epoch
    ret += '2000 '
    # Proper motion RA and DEC
    ret += 'pmra=' + str(row[DS_PMRA]) + ' '
    ret += 'pmdec=' + str(row[DS_PMDEC]) + ' '
    # V Mag
    ret += 'vmag=' + str(row[DS_VMAG]) + ' '
    # T Exp
    ret += 'texp=' + str(int(row[DS_EXPT])) + ' '
    # I2
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=' + str(t.hour) + ' '
    ret += 'utm=' + str(t.minute) + ' '
    # Exp Count
    if row[DS_COUNTS] > 3e9:
        ret += 'expcount=%.3g' % (3e9) + ' '
    else:
        ret += 'expcount=%.3g' % (row[DS_COUNTS]) + ' '
    # Decker
    ret += 'decker=%s ' % (decker)
    # do flag
    if do_flag:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    ret += 'count=' + str(int(row[DS_NSHOTS]))

    ret += ' foc=' + str(int(focval))

    if owner != '':
        ret += ' owner=' + str(owner)

    return ret

def calc_elevations(stars, observer):
    els = []
    for s in stars:
        observer.date = ephem.Date(observer.date)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        els.append(cur_el)
    return np.array(els)

def compute_datetime(ctime):
    if type(ctime) == float:
        dt = datetime.utcfromtimestamp(int(ctime))
    elif type(ctime) == datetime:
        dt = ctime
    elif type(ctime) == ephem.Date:
        dt = ctime.datetime()
    else:
        #punt and use current UT
        dt = datetime.utcfromtimestamp(int(time.time()))
    return dt


def makeAPFObs(dt,horizon=str(TARGET_ELEVATION_MIN)):
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = horizon
    apf_obs.date = dt

    return apf_obs

def compute_sunset_n_rise(dt,horizon='0'):
    # computes time in seconds before sunset
    apf_obs = makeAPFObs(dt,horizon=horizon)
    sunset = apf_obs.next_setting(ephem.Sun())
    sunset -= ephem.Date(dt)
    sunset *= 86400.0 # convert to seconds
    
    sunrise = apf_obs.next_rising(ephem.Sun())
    sunrise -= ephem.Date(dt)
    sunrise *= 86400.0 # convert to seconds
    return sunset, sunrise

def compute_sunset(dt,horizon='0'):

    sunset, sunrise = compute_sunset_n_rise(dt,horizon=horizon)
    return sunset
    
def compute_sunrise(dt,horizon='0'):
    sunset, sunrise = compute_sunset_n_rise(dt,horizon=horizon)
    return sunrise


def templateConditions(moon, seeing, slowdown):

    if seeing < 15 and slowdown < 0.5:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon.alt),echo=True,level='debug')
        if moon.phase < 50 and float(moon.alt) < 0:
            return True
        elif moon.phase < 25 and float(moon.alt) < 0.7:
            return True
        else:
            return False
    else:
        return False

def findClosest(ras,decs,ra,dec):


    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()

    return min_ind

def makeTempRow(star_table,ind,bstar=False):

    row = []

    row.append(star_table[ind, DS_RA])
    row.append( star_table[ind, DS_DEC])
    row.append(star_table[ind, DS_PMRA])
    row.append(star_table[ind, DS_PMDEC])
    row.append(star_table[ind, DS_VMAG])
    row.append(1200)
    row.append(1e9)
    row.append( star_table[ind, DS_APFPRI])
    row.append(0)
    if bstar:
        row.append(2)
    else:
        if star_table[ind, DS_VMAG] > 10:
            row.append(9)
        elif star_table[ind, DS_VMAG] < 8:
            row.append(5)
        else:
            row.append(7)
    return row

def enoughTime(star_table,stars,idx,row,apf_obs,dt):
    tot_time = row[DS_NSHOTS]*row[DS_EXPT]
    tot_time += 210 + (2*40 + 40*(row[DS_NSHOTS]-1)) # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, star_elevations, fin_els = Visible.is_visible(apf_obs,[stars[idx]],[tot_time])
    time_left_before_sunrise = compute_sunrise(dt,horizon='-9')

    try:
        apflog( "enoughTime(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time,  time_left_before_sunrise),echo=True)
    except:
        apflog("enoughTime(): cannot log times!?!",echo=True)
        
    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False
        
    
def findBstars(snames,star_table,idx, bstars):

    near_idx = findClosest(star_table[:,DS_RA][bstars],star_table[:,DS_DEC][bstars],star_table[idx,DS_RA],star_table[idx,DS_DEC])
    row = makeTempRow(star_table[bstars],near_idx,bstar=True)
    
    end_idx = findClosest(star_table[:,DS_RA][bstars],star_table[:,DS_DEC][bstars],(star_table[idx,DS_RA]+45*np.pi/180.),star_table[idx,DS_DEC])
    finrow = makeTempRow(star_table[bstars],end_idx,bstar=True)
    
    return snames[near_idx],row,snames[end_idx],finrow


def makeResult(stars,star_table,flags,totexptimes,sn,dt,idx,focval=0):
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV']     = star_table[idx, DS_BV]
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NEXP'] = star_table[idx, DS_NSHOTS]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['NAME']   = sn[idx]
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['DECKER'] = flags['decker'][idx]
    res['isTemp'] =    False
    res['owner'] =    flags['owner'][idx]
    res['SCRIPTOBS'] = []
    scriptobs_line = makeScriptobsLine(sn[idx], star_table[idx,:], flags['do'][idx], dt, decker=flags['decker'][idx], I2=flags['I2'][idx], owner=flags['owner'][idx],focval=focval) + " # end"
    res['SCRIPTOBS'].append(scriptobs_line)
    return res


def getNext(ctime, seeing, slowdown, bstar=False,template=False,sheetns=["Bstars"],owner='public',outfn="googledex.dat",toofn="too.dat",outdir=None,focval=0,inst=''):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    if not outdir:
        outdir = os.getcwd()

    dt = compute_datetime(ctime)

    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['mode']=''
    config['obsblock']=''
    config['Bstar']='N'
    config['owner']=owner
    config['inst']='levy'
    

    apflog( "getNext(): Finding target for time %s" % (dt),echo=True)

    if slowdown > SLOWDOWN_MAX:
        apflog( "getNext(): Slowndown value of %f exceeds maximum of %f at time %s" % (slowdown,SLOWDOWN_MAX,dt),echo=True)
        return None


    try:
        apfguide = ktl.Service('apfguide')
        stamp = apfguide['midptfin'].read(binary=True)
        ptime = datetime.utcfromtimestamp(stamp)
    except:
        if type(dt) == datetime:
            ptime = dt
        else:
            ptime = datetime.utcfromtimestamp(int(time.time()))

    observed = ParseUCOSched.update_local_googledex(ptime,googledex_file=os.path.join(outdir,"googledex.dat"), observed_file=os.path.join(outdir,"observed_targets"))

    # List of targets already observed

    global last_objs_attempted
    try:
        lastline = ktl.read("apftask","SCRIPTOBS_LINE")
        if not bstar:             # otherwise from previous night
            lastobj = lastline.split()[0]
        else:
            lastobj = None

    except:
        lastobj = None

    if lastobj:
        if lastobj not in observed and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)
            
            apflog( "getNext(): Last objects attempted %s" % (last_objs_attempted),echo=True)

            if len(last_objs_attempted) > 5:
                apflog( "getNext(): 5 failed acquisition attempts",echo=True)
                last_objs_attempted = []
                return None
        else:
            last_objs_attempted = []
            # we had a succes so we are zeroing this out

            
    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###

    apf_obs = makeAPFObs(dt)
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180.

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    do_templates = template and templateConditions(moon, seeing, slowdown)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    apflog("getNext(): Parsing the Googledex...",echo=True)
    sn, star_table, flags, stars = ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)
    sn = np.array(sn)
    deckers = np.array(flags['decker'])
    bstar_array = np.array(flags['Bstar'])
    mode = np.array(flags['mode'])
    obsblock = np.array(flags['obsblock'])
    targNum = len(sn)
    
    apflog("getNext(): Parsed the Googledex...",echo=True)


    apflog("getNext(): Finding B stars",echo=True)
    # Note which of these are B-Stars for later.
    bstars = (bstar_array == 'Y')|(bstar_array == 'y')

    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    available = np.ones(targNum, dtype=bool)
    totexptimes = np.zeros(targNum, dtype=float)
    cur_elevations = np.zeros(targNum, dtype=float)
    scaled_elevations = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?

    apflog("getNext(): Culling stars behind the moon",echo=True)
    moon_check = moonDist > minMoonDist
    available = available & moon_check

    #    totobs_check = (star_table[:,DS_NOB] < star_table[:,DS_TOT]) | (star_table[:,DS_TOT] <= 0)
    #    available = available & totobs_check

    # We just need a B star, so restrict our math to those
    if bstar:
        
        apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars

        f = available
        
        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations = Visible.is_visible(apf_obs, fstars, [400]*len(bstars[f]))

        available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]

        star_table[available, DS_COUNTS] = 1e9
        star_table[available, DS_EXPT] = 900
        star_table[available, DS_NSHOTS] = 2
        totexptimes[available] = 400

    # Just need a normal star for observing
    else:
        # Available and not a BStar

        apflog("getNext(): Culling B stars",echo=True)
        available = available & np.logical_not(bstars)

        # has the star been observed - commented out as redundant with cadence
        if len(last_objs_attempted)>0:
            for n in last_objs_attempted:
                attempted = (sn == n)
                available = available & np.logical_not(attempted) # Available and not observed

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function


        apflog("getNext(): Computing exposure times",echo=True)
        exp_times = star_table[available,DS_NSHOTS] * star_table[available,DS_EXPT]
        exp_counts = star_table[available,DS_COUNTS]

        totexptimes[available] += exp_times
        
        # Is the exposure time too long?
        apflog("getNext(): Removing really long exposures",echo=True)
        time_check = exp_times < TARGET_EXPOSURE_TIME_MAX

        available = available & time_check

        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,available) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.is_visible_se(apf_obs, fstars, exp_times)
        
        available = available & vis
        

    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None


    final_priorities = computePriorities(star_table,available,dt,flags)

    cadence_check = (ephem.julian_date(dt) - star_table[:, DS_LAST]) / star_table[:, DS_CAD]
    good_cadence = np.where(cadence_check >  1.0, True, False)
    good_cadence_available = available & good_cadence

    if any(good_cadence_available):
        try:
            pri = max(final_priorities[good_cadence_available])
            sort_i = (final_priorities == pri) & good_cadence_available
        except:
            pri = max(final_priorities[available])
            sort_i = (final_priorities == pri) & available
    elif any(available):
        apflog( "getNext(): No new stars available, going back to the previously observed list.",level="warn",echo=True)
        pri = max(final_priorities[available])
        sort_i = (final_priorities == pri) & available
    else:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

    starstr = "getNext(): star table available: %s" % (sn[sort_i])
    apflog(starstr,echo=True)

    starstr = "getNext(): star table available priorities: %s" % (final_priorities[sort_i])
    apflog(starstr,echo=True)

    if bstar:
        sort_j = cur_elevations[sort_i].argsort()[::-1]
        focval=2
    else:
        sort_j = scaled_elevations[sort_i].argsort()[::-1]
        cstr= "getNext(): cadence check: %s" % (cadence_check[sort_i][sort_j][0])
        apflog(cstr,echo=True)

    t_n = sn[sort_i][sort_j][0]

    elstr= "getNext(): star elevations %s" % (cur_elevations[sort_i][sort_j])
    apflog(elstr,echo=True)

    t_n = sn[sort_i][sort_j][0]

    apflog("getNext(): selected target %s" % (t_n) )

    idx, = np.where(sn == t_n)
    idx = idx[0]

    stars[idx].compute(apf_obs)
    cstr= "getNext(): cadence check: %f (%f %f %f)" % (((ephem.julian_date(dt) - star_table[idx, DS_LAST]) / star_table[idx, DS_CAD]), ephem.julian_date(dt), star_table[idx, DS_LAST], star_table[idx, DS_CAD])
    apflog(cstr,echo=True)

    res =  makeResult(stars,star_table,flags,totexptimes,sn,dt,idx,focval=focval)
    if do_templates and flags['template'][idx] == 'N' and flags['I2'][idx] == 'Y':
        bname,brow,bnamefin,browfin = findBstars(sn,star_table,idx,bstars)
        row = makeTempRow(star_table,idx)
        if enoughTime(star_table,stars,idx,row,apf_obs,dt):
            bline = makeScriptobsLine(bname,brow,'',dt,decker="N",I2="Y", owner='public',focval=2)
            line  = makeScriptobsLine(sn[idx],row,flags['do'][idx],dt,decker="N",I2="N", owner=flags['owner'][idx])
            bfinline = makeScriptobsLine(bnamefin,browfin,'',dt,decker="N",I2="Y", owner='public',focval=2)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            apflog("Attempting template observation of %s" % (sn[idx]),echo=True)

    return res

if __name__ == '__main__':

#    sheetn=["2018B"]
    sheetn="Bstars_test,A004_PRobertson_test,A014_SVogt_2019B_test"

    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 7.99, 0.4, bstar=True,sheetns=sheetn.split(","))
    ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()
    starttime += 400
    for i in range(5):

        result = getNext(starttime, 7.99, 0.4, bstar=False,sheetns=sheetn,template=True)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print("Get None target")
        else:
            for k in result:
                print(k, result[k])
        while len(result["SCRIPTOBS"]) > 0:
            ot = open(otfn,"a")
            ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
            ot.close()
            starttime += result["EXP_TIME"]

    print("Done")
