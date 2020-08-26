# UCOScheduler_V1.py
from __future__ import print_function
import os
import sys
import time
import re
from datetime import datetime, timedelta

import numpy as np
import astropy
import astropy.io
import astropy.coordinates
import astropy.time
import astroplan

from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter_Rate, getEXPTime
import ParseUCOSched
import Coords
from SchedulerConsts import * # I know

try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import Visible

# some globals

last_objs_attempted = []

# some constants
ACQUIRE = 'A'
BLANK = 'B'
FIRST = '1'
LAST = 'L'

def computePriorities(star_table,available,cur_dt,observed=None,hour_table=None,rank_table=None):
    # make this a function, have it return the current priorities, than change references to the star_table below into references to the current priority list
    new_pri = np.zeros_like(star_table['APFpri'])
    new_pri += star_table['APFpri']

    if hour_table is not None:
        too_much = hour_table['cur']  > hour_table['tot']
        done_sheets = hour_table['sheetn'][too_much]
    else:
        done_sheets = []

    if rank_table is not None:
        for sheetn in rank_table['sheetn']:
            if sheetn not in done_sheets:
                cur = star_table['sheetn'] == sheetn
                new_pri[cur] += rank_table['rank'][rank_table['sheetn'] == sheetn]
    
        
    return new_pri

def updateHourTable(hour_table,observed,ctime,outfn='hour_table'):
    '''
    updateHourTableobserved_logs,outfn='hour_table')

    Updates hour_table with history of observations.
    
    '''


    hours = dict()

    dt = ctime.datetime
    
    # observed objects have lists as attributes
    # in reverse time order, so most recent target observed is first.
    nobj = len(observed.names)
    for i in range(0,nobj):
            own = observed.owners[i]
            if own not in hours.keys():
                    hours[own] = 0.0
                    
    cur = dt
    for i in range(0,nobj):
            hr, mn = observed.times[i]
            prev = datetime(dt.year,dt.month,dt.day,hr,mn)
            diff = cur - prev
            hourdiff = (diff.days * 24 + diff.seconds / 3600.)
            if hourdiff > 0:
                hours[observed.owners[i]] += hourdiff
            cur = prev

    for ky in hours.keys():
        hour_table['cur'][hour_table['sheetn'] == ky] = hours[ky]

    try:
        hour_table.write(outfn,format='ascii',overwrite=True)
    except Exception as e:
        apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)

    return hour_table
    

def makeHourTable(sheet_table_name,dt,outfn='hour_table',outdir=None,frac_fn='frac_table',hour_constraints=None):

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)

    if os.path.exists(outfn):
        hour_table =  astropy.table.Table.read(outfn,format='ascii')
        return hour_table

    frac_fn = os.path.join(outdir,frac_fn)
    if os.path.exists(frac_fn):
        frac_table = astropy.table.Table.read(frac_fn,format='ascii')
    else:
        sheetns, fracs = ParseUCOSched.parseFracTable(sheet_table_name=sheet_table_name,outfn=frac_fn)
        frac_table = []
        for i in range(0,len(fracs)):
            frow = []
            frow.append(sheetns[i])
            frow.append(fracs[i])
            frac_table.append(frow)
        frac_table = astropy.table.Table(rows=frac_table,names=['sheetn','frac'])
        try:
            frac_table.write(frac_fn,format='ascii')
        except Exception as e:
            apflog("Cannot write table %s: %s" % (frac_fn,e),level='error',echo=True)


    hour_table= astropy.table.Table(frac_table,names=['sheetn','frac'])

    sunset,sunrise = computeSunsetRise(dt,horizon='-9')
    if sunrise < sunset:
        sunrise += 86400
    tot = sunrise - sunset
    tot /= 3600.

    hour_table['tot'] =tot*hour_table['frac']
    hour_table['cur'] =0.0*hour_table['frac']

    if hour_constraints is not None:
        if 'runname' in hour_constraints.keys() and 'left' in hour_constraints.keys():
            for runname in hour_constraints['runname']:
                if hour_constraints['left'][hour_constraints['runname']==runname] < hour_table['tot'][hour_table['sheetn']==runname]:
                     hour_table['tot'][hour_table['sheetn']==runname] = hour_constraints['left'][hour_constraints['runname']==runname]
                     

    try:
        hour_table.write(outfn,format='ascii')
    except Exception as e:
        apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)
    return hour_table

def makeRankTable(sheet_table_name,outfn='rank_table',outdir=None):

    if not outdir :
        outdir = os.getcwd()

    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn):
        rank_table = astropy.table.Table.read(outfn,format='ascii')
    else:
        sheetns, ranks = ParseUCOSched.parseRankTable(sheet_table_name=sheet_table_name)

        rank_table= astropy.table.Table([sheetns,ranks],names=['sheetn','rank'])
        try:
            rank_table.write(outfn,format='ascii')
        except Exception as e:
            apflog("Cannot write table %s: %s" % (outfn,e),level='error',echo=True)

    return rank_table


def makeScriptobsLine(star_table_row, t, decker="W", I2="Y", owner='public', focval=0, coverid='',temp=False):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(star_table_row, t, decker="W",I2="Y")
    star_table_row -contains all of the data needed for the line except

    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    temp - a boolean for whether or not this is a template observation
    """

    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = str(star_table_row['name']) + ' '
    # Add the RA as three elements, HR, MIN, SEC
    rastr = "%s %s %s" % (star_table_row['RA hr'],star_table_row['RA min'],star_table_row['RA sec'])
    ret += rastr + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = "%s %s %s" % (star_table_row['Dec deg'],star_table_row['Dec min'],star_table_row['Dec sec'])

    ret += decstr + ' '
    # Epoch
    ret += '2000 '
    # Proper motion RA and DEC
    ret += 'pmra=' + str(star_table_row['pmRA']) + ' '
    ret += 'pmdec=' + str(star_table_row['pmDEC']) + ' '
    # V Mag
    ret += 'vmag=' + str(star_table_row['Vmag']) + ' '

    # T Exp
    if temp:
        ret += 'texp=' + str(1200) + ' '
    else:
        ret += 'texp=' + str(int(star_table_row['texp'])) + ' '

    # I2
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=' + str(t.datetime.hour) + ' '
    ret += 'utm=' + str(t.datetime.minute) + ' '
    # Exp Count
    if star_table_row['expcount'] > EXP_LIM:
        ret += 'expcount=%.3g' % (EXP_LIM) + ' '
    elif temp:
        ret += 'expcount=%.3g' % (1e9) + ' '
    else:
        ret += 'expcount=%.3g' % (star_table_row['expcount']) + ' '
    # Decker
    ret += 'decker=%s ' % (decker)
    # do flag
    if star_table_row['do']:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    if temp:
        if star_table_row['Vmag'] > 10:
            count = 9
        elif star_table_row['Vmag'] < 8:
            count = 5
        else:
            count = 7
    else:
        count = int(star_table_row['APFnshots'])

    ret += 'count=' + str(count)

    ret += ' foc=' + str(int(focval))

    if owner != '':
        if owner == 'RECUR_A100':
            owner = 'public'
        ret += ' owner=' + str(owner)

    if coverid != '':
        ret += ' coverid=' + str(coverid)

    if star_table_row['mode'] != None:
        if star_table_row['mode'] == BLANK:
            ret += ' blank=Y'
        elif star_table_row['mode'] == ACQUIRE:
            ret += ' guide=Y'
    else:
        ret += ''

    raoff  = star_table_row['raoff']
    decoff = star_table_row['decoff']
    if raoff == 'None':
        raoff = ''
    if decoff == 'None':
        decoff = ''
#    if raoff is not '' and decoff is not '':
#        ret += ' raoff=' + str(raoff) + ' decoff=' + str(decoff)

    return str(ret)

def computeDatetime(ctime):
    if type(ctime) == float:
        dt = astropy.time.Time(ctime,format='unix')
    elif type(ctime) == datetime:
        dt = astropy.time.Time(ctime)
    elif type(ctime) == astropy.time.Time:
        dt = ctime
    else:
        #punt and use current UT
        dt = astropy.time.Time.now()
    return dt


def computeSunsetRise(compute_time,horizon='0'):
    # computes time in seconds before sunset

    horizon_deg = float(horizon) * astropy.units.degree
    
    apf_obs = Visible.makeAPFObs()
    sunset = apf_obs.sun_set_time(compute_time,which='next',horizon=horizon_deg)
    sunset_sec = sunset.jd - compute_time.jd
    sunset_sec *= 86400.0 # convert to seconds

    sunrise = apf_obs.sun_rise_time(compute_time,which='next',horizon=horizon_deg)
    sunrise_sec = sunrise.jd - compute_time.jd
    sunrise_sec *= 86400.0 # convert to seconds
    return sunset_sec, sunrise_sec

def computeSunset(dt,horizon='0'):

    sunset, sunrise = computeSunsetRise(dt,horizon=horizon)
    return sunset

def computeSunrise(dt,horizon='0'):
    sunset, sunrise = computeSunsetRise(dt,horizon=horizon)
    return sunrise


def conditionCuts(moon_phase,seeing,slowdown,star_table):
    """ available = conditionCuts(moon, seeing, slowdown, star_table)

    Checks if columns are in the star_table, then cuts on those, returns a boolean numpy array

    available - Boolean numpy array of available targets

    moon_phase - phase value from astroplan, ranges from 0 (full) to pi (new) as an angle in radians
    seeing - size in pixels
    transparency - magnitudes of extinction

    """

    if 'seeing' in star_table.colnames:
        available = star_table['seeing']/0.109 > seeing
    else:
        available = np.ones(len(star_table['ra'], dtype=bool))

    if 'moon' in star_table.colnames:
        available = (star_table['moon'] > moon_phase) & available

    if 'transparency' in star_table.colnames:
        ext = 2.5 * np.log10(slowdown)
        available = (star_table['transparency'] > ext) & available


    return available


def templateConditions(apf_obs, dt, seeing, slowdown):
    """ istrue = templateCondition(obs, date, seeing, slowdown)

    Checks to see if moon, seeing and slowdown factor are within template conditions

    istrue - a simple Boolean

    seeing - size in pixels
    slowdown - relative to clear

    """

    moon_po<s   = apf_obs.moon_altaz(dt)
    moon_phase = apf_obs.moon_phase(dt)

    
    if seeing < 15 and slowdown < 1.25:
        apflog("moon.phase=%.2f moon.alt=%.2f" % (moon.phase,moon_pos.alt),echo=True,level='debug')
        if moon_phase.value > np.pi/2 and moon_pos.alt.value < 0:
            return True
        elif moon_phase.value > 3*np.pi/2 and moon_pos.alt.value < 45:
            return True
        else:
            return False
    else:
        return False

def findClosest(ras,decs,ra,dec):


    distances = np.sqrt((ra - ras)**2 + (dec - decs)**2)

    min_ind = distances.argmin()

    return min_ind


def enoughTime(star_table,stars,idx,apf_obs,dt):
    tot_time = star_table['APFnshots'][idx]*star_table['texp'][idx]
    tot_time += 210 + (2*40 + 40*(star_table['APFnshots'][idx]-1)) + 2400 # two B star exposures + three 70 second acquisitions and the actual observation readout times
    vis, star_elevations, fin_els = Visible.visible(apf_obs,[stars[idx]],[tot_time])
    time_left_before_sunrise = computeSunrise(dt,horizon='-9')

    try:
        apflog( "enoughTime(): time for obs= %.1f  time until sunrise= %.1f " % (tot_time,  time_left_before_sunrise),echo=True)
    except:
        apflog("enoughTime(): cannot log times!?!",echo=True)

    if tot_time < time_left_before_sunrise  and vis and time_left_before_sunrise < 14*3600.:
        return True
    else:
        return False


def findBstars(star_table,idx, bstars):

    near_idx = findClosest(star_table['ra'][bstars],star_table['dec'][bstars],star_table['ra'][idx],star_table['dec'][idx])

    end_idx = findClosest(star_table['ra'][bstars],star_table['dec'][bstars],(star_table['ra'][idx]+15*np.pi/180.),star_table['dec'][idx])


    return near_idx,end_idx


def makeObsBlock(star_table, idx, dt, focval):

    rv = []

    cur_obsblock = star_table['obsblock'][idx]

    allinblock = (star_table['obsblock'] == cur_obsblock)
    allinblock = allinblock & (star_table['sheetn'] == star_table['sheetn'][idx])

    if np.any(star_table['mode'][allinblock] == FIRST):
        first = (star_table['mode'][allinblock] == FIRST)
    elif np.any(star_table['mode'][allinblock] == ACQUIRE):
        first = (star_table['mode'][allinblock] == ACQUIRE)
    else:
        first = None

    if np.any(star_table['mode'][allinblock] == LAST):
        last = (star_table['mode'][allinblock] == LAST)
    else:
        last = None

    rest = (star_table['mode'][allinblock] != FIRST)
    rest = rest & (star_table['mode'][allinblock] != ACQUIRE)
    rest = rest & (star_table['mode'][allinblock] != LAST)
    rest_idxs, = np.where(rest)


    if np.any(first):
        first_idxs, = np.where(first)
        for idx in first_idxs:
            scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)

    for idx in rest_idxs:
        scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                               owner=star_table['sheetn'][allinblock][idx], \
                                               I2=star_table['I2'][allinblock][idx], focval=focval)
        rv.append(scriptobs_line)

    if np.any(last):
        last_idxs, = np.where(last)
        for idx in last_idxs:
            scriptobs_line = makeScriptobsLine(star_table[allinblock][idx], dt, decker=star_table['decker'][allinblock][idx], \
                                                owner=star_table['sheetn'][allinblock][idx], \
                                                I2=star_table['I2'][allinblock][idx], focval=focval)
            rv.append(scriptobs_line)


    rv.reverse()
    rv[0] += ' # obsblock=%s end' % (cur_obsblock)
    return(rv)

def makeResult(stars,star_table,totexptimes,dt,idx,focval=0,bstar=False,mode=''):
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table['pmRA'][idx]
    res['PM_DEC'] = star_table['pmDEC'][idx]
    res['VMAG']   = star_table['Vmag'][idx]
    res['BV']     = star_table['B-V'][idx]
    res['COUNTS'] = star_table['expcount'][idx]
    res['EXP_TIME'] = star_table['texp'][idx]
    res['NEXP'] = star_table['APFnshots'][idx]
    res['TOTEXP_TIME'] = totexptimes[idx]
    res['NAME']   = star_table['name'][idx]
    res['PRI']    = star_table['APFpri'][idx]
    res['DECKER'] = star_table['decker'][idx]
    res['I2']     = star_table['I2'][idx]
    res['isTemp'] =    False
    res['isBstar'] =    bstar
    res['mode']   =   ''
    res['owner'] =    star_table['sheetn'][idx]

    res['obsblock'] = star_table['obsblock'][idx]

    res['SCRIPTOBS'] = []
    scriptobs_line = makeScriptobsLine(star_table[idx], dt, decker=res['DECKER'], owner=res['owner'], I2=star_table['I2'][idx], focval=focval)
    scriptobs_line = scriptobs_line + " # end"
    res['SCRIPTOBS'].append(scriptobs_line)

    return res

def lastAttempted(observed):
    global last_objs_attempted
    try:
        lastresult = ktl.read("apftask","SCRIPTOBS_LINE_RESULT",binary=True)
    except:
        return []

    if lastresult == 2:
        try:
            lastline = ktl.read("apftask","SCRIPTOBS_LINE")
            lastobj = lastline.split()[0]
        except:
            lastobj = None
    else:
        lastobj = None

    if lastobj:
        if lastobj not in observed.names and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)

            apflog( "getNext(): Last objects attempted %s" % (last_objs_attempted),echo=True)


        else:
            last_objs_attempted = []
            # we had a succes so we are zeroing this out

    return last_objs_attempted

def behindMoon(apf_obs,dt,ras,decs):
    moon_pos   = apf_obs.moon_altaz(dt)
    moon_phase_ang = apf_obs.moon_phase(dt)

    moon_phase = 1.0 - (moon_phase_ang.value / np.pi) 
    
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ( moon_phase  * md) + TARGET_MOON_DIST_MIN
    moonDist = np.sqrt((moon_pos.icrs.ra.value - ras)**2 + (moon_pos.icrs.dec.value - decs)**2)

    apflog("getNext(): Culling stars behind the moon",echo=True)
    moon_check = moonDist > minMoonDist

    return moon_check

def getNext(ctime, seeing, slowdown, bstar=False,template=False,sheetns=["RECUR_A100"],owner='public',outfn="googledex.dat",toofn="too.dat",outdir=None,focval=0,inst='',rank_sheetn='rank_table',frac_sheet=None):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    if not outdir:
        outdir = os.getcwd()

    dt = computeDatetime(ctime) # this is an astropy.time.Time object 

    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['mode']=''
    config['obsblock']=''
    config['Bstar']='N'
    config['owner']=owner
    config['inst']='levy'
    config['raoff'] = ''
    config['decoff'] = ''


    apflog( "getNext(): Finding target for time %s" % (dt.isot),echo=True)

    if slowdown > SLOWDOWN_MAX:
        apflog( "getNext(): Slowndown value of %f exceeds maximum of %f at time %s" % (slowdown,SLOWDOWN_MAX,dt.isot),echo=True)
        return None


    try:
        apfguide = ktl.Service('apfguide')
        stamp = apfguide['midptfin'].read(binary=True)
        ptime = datetime.utcfromtimestamp(stamp)
    except:
        if type(dt) == astropy.time.Time:
            ptime = dt
        else:
            ptime =  astropy.time.Time.now()

    apflog("getNext(): Updating star list with previous observations",echo=True)
    observed, star_table = ParseUCOSched.updateLocalStarlist(ptime,outfn=outfn,toofn=toofn,observed_file="observed_targets")

    hour_table = None
    if frac_sheet is not None:
        hour_table = makeHourTable(frac_sheet,dt)
        hour_table = updateHourTable(hour_table,observed,dt)
    
        
    
    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians

    if star_table is None:
        apflog("getNext(): Parsing the star list",echo=True)
        star_table, stars = ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)
    else:
        stars = ParseUCOSched.genStars(star_table)
    targNum = len(stars)

    # List of targets already observed

    last_objs_attempted = lastAttempted(observed)
    if len(last_objs_attempted) > 5:
        apflog( "getNext(): 5 failed acquisition attempts",echo=True)
        last_objs_attempted = []

    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###

    apf_obs = Visible.makeAPFObs()
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180.

    # Calculate the moon's location

    do_templates = template and templateConditions(apf_obs, dt, seeing, slowdown)

    apflog("getNext(): Parsed the Googledex...",echo=True)

    apflog("getNext(): Finding B stars",echo=True)
    # Note which of these are B-Stars for later.
    bstars = (star_table['Bstar'] == 'Y')|(star_table['Bstar'] == 'y')

    totexptimes = np.zeros(targNum, dtype=float)
    totexptimes = star_table['texp'] * star_table['APFnshots'] + 40 * (star_table['APFnshots']-1)

    available = np.ones(targNum, dtype=bool)
    cur_elevations = np.zeros(targNum, dtype=float)
    scaled_elevations = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?
    moon_check = behindMoon(apf_obs,dt,star_table['ra'],star_table['dec'])
    available = available & moon_check
    if len(last_objs_attempted)>0:
        for n in last_objs_attempted:
            attempted = (star_table['name'] == n)
            available = available & np.logical_not(attempted) # Available and not observed
            
    cadence_check = ( - star_table['lastobs']) 
    good_cadence = cadence_check >  star_table['APFcad']
    available = available & good_cadence

    # We just need a B star, so restrict our math to those
    if bstar:

        apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars

        f = available

        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations = Visible.visible(apf_obs, dt, fstars, [400]*len(bstars[f]))

        available[f] = available[f] & vis
        cur_elevations[f] += star_elevations[vis]


    # Just need a normal star for observing
    else:
        # Available and not a BStar

        apflog("getNext(): Culling B stars",echo=True)
        available = available & np.logical_not(bstars)

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function

        apflog("getNext(): Computing exposure times",echo=True)
        exp_counts = star_table['expcount']

        # Is the exposure time too long?
        apflog("getNext(): Removing really long exposures",echo=True)

        time_left_before_sunrise = computeSunrise(dt,horizon='-9')
        maxexptime = TARGET_EXPOSURE_TIME_MAX
        if maxexptime > time_left_before_sunrise:
            maxexptime = time_left_before_sunrise
        if maxexptime < TARGET_EXPOSURE_TIME_MIN:
            maxexptime = TARGET_EXPOSURE_TIME_MIN # this will try a target in case we get lucky
            
        time_check = totexptimes <= maxexptime

        available = available & time_check

        apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,available) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.visibleSE(apf_obs, dt, fstars, totexptimes[available],shiftwest=True)
        currently_available = available
        currently_available[available] = currently_available[available] & vis

        cur_elevations[available] += star_elevations[vis]
        scaled_elevations[available] += scaled_els[vis]

        
        if slowdown > SLOWDOWN_THRESH or seeing > SEEING_THRESH:
            bright_enough = star_table['Vmag'] < SLOWDOWN_VMAG_LIM
            available = available & bright_enough

    # Now just sort by priority, then cadence. Return top target
    if len(star_table['name'][available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None


    final_priorities = computePriorities(star_table,available,dt,rank_table=makeRankTable(rank_sheetn),hour_table=hour_table,observed=observed)

    try:
        pri = max(final_priorities[available])
        sort_i = (final_priorities == pri) & available
    except:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

    if bstar:
        sort_j = cur_elevations[sort_i].argsort()[::-1]
        focval=2
    else:
        sort_j = scaled_elevations[sort_i].argsort()[::-1]

    allidx, = np.where(sort_i)
    idx = allidx[sort_j][0]
        
    t_n = star_table['name'][idx]
    o_n = star_table['sheetn'][idx]
    p_n = final_priorities[idx]

    apflog("getNext(): selected target %s for program %s at priority %.0f" % (t_n,o_n,p_n) )
    nmstr= "getNext(): star names %s" % (np.asarray(star_table['name'][sort_i][sort_j]))
    shstr= "getNext(): star sheet names %s" % (np.asarray(star_table['sheetn'][sort_i][sort_j]))
    if bstar:
        elstr= "getNext(): Bstar current elevations %s" % (cur_elevations[sort_i][sort_j])
    else:
        elstr= "getNext(): star scaled elevations %s" % (scaled_elevations[sort_i][sort_j])
    apflog(nmstr,echo=True)
    apflog(shstr,echo=True)
    apflog(elstr,echo=True)

    stars[idx].compute(apf_obs)

    res =  makeResult(stars,star_table,totexptimes,dt,idx,focval=focval,bstar=bstar,mode=config['mode'])
    if do_templates and star_table['Template'][idx] == 'N' and star_table['I2'][idx] == 'Y':
        bidx,bfinidx = findBstars(star_table,idx,bstars)

        if enoughTime(star_table,stars,idx,apf_obs,dt):
            bline = makeScriptobsLine(star_table[bidx],dt,decker="N",I2="Y", owner=res['owner'],focval=2)
            line  = makeScriptobsLine(star_table[idx],dt,decker="N",I2="N", owner=res['owner'],temp=True)
            bfinline = makeScriptobsLine(star_table[bfinidx],dt,decker="N",I2="Y",owner=res['owner'],focval=0)
            res['SCRIPTOBS'] = []
            res['SCRIPTOBS'].append(bfinline + " # temp=Y end")
            res['SCRIPTOBS'].append(line + " # temp=Y")
            res['SCRIPTOBS'].append(bline + " # temp=Y")
            res['isTemp'] = True
            apflog("Attempting template observation of %s" % (star_table['name'][idx]),echo=True)

    return res

if __name__ == '__main__':

    dt = datetime.now()
    
    frac_tablen='2020B_frac'
    hour_table = makeHourTable(frac_tablen,dt)
    
    rank_tablen='2020B_ranks'
    rank_table = makeRankTable(rank_tablen)

#    sheetn=["2018B"]
    sheetn="RECUR_A100,2020B_A000,2020B_A008,2020B_A009,2020B_A010"

    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 7.99, 0.4, bstar=True,sheetns=sheetn.split(","),rank_sheetn=rank_tablen,frac_sheet=frac_tablen)
    ot.write("%s\n" % (result["SCRIPTOBS"].pop()))
    ot.close()
    starttime += 400
    for i in range(5):

        result = getNext(starttime, 7.99, 0.4, bstar=False,sheetns=sheetn,template=True,rank_sheetn=rank_tablen,frac_sheet=frac_tablen)
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
