#! /usr/bin/python
# UCSCScheduler.py

import calendar
from datetime import datetime, timedelta
import ephem
import gspread
import numpy as np
import os
import pickle
import sys
import time
from x_gaussslit import *

from apflog import *

# A few constants to make accessing the star table more readable
DS_RA     = 0
DS_DEC    = 1
DS_PMRA   = 2
DS_PMDEC  = 3
DS_VMAG   = 4
DS_EXPT   = 5
DS_COUNTS = 6
DS_APFPRI = 7
DS_CAD    = 8
DS_NSHOTS = 9
DS_LAST   = 10
DS_BV     = 11
DS_ERR    = 12

MAX_EXPTIME = 1200.

slit_size = {'M': (9.52, 74.07), 'W': (9.26, 27.78), 'N': (4.63, 74.07), 'B': (18.52, 74.07),
             'S': (6.94, 74.07), 'P': (1, 1) }

def parseStarsAPF(filename):
    """Parse out the starsAPF file. Returns an (n x m) numpy array."""
    star_table = []
    names = []
    do_flag = []
    stars = []
    with open(filename,'r') as f:
        for line in f:
            ls = line.split()
            row = []
            names.append(ls[0])
            row.append(getRARad(ls[1], ls[2], ls[3]))
            row.append(getDECRad(ls[5], ls[6], ls[7], ls[4]))
            for v in ls[8:-1]: row.append(float(v))
            if ls[-1] == '_':
                do_flag.append(False)
            else:
                do_flag.append(True)
            star_table.append(row)
            # Try creating a pyephem object for each star
            star = ephem.FixedBody()
            star._ra = ephem.hours(":".join(ls[1:4]))
            star._dec = ephem.degrees(":".join([ls[4]+ls[5], ls[6], ls[7]]))
            stars.append(star)

    return (names, np.array(star_table), do_flag, stars)

def parseStarlist(starlist):
    """ Parse a scriptobs-compatible starlist for the scheduler."""
    names = []
    lines = []
    stars = []
    star_table = []
    try:
        f = open(starlist,'r')
    except IOError:
        apflog( "Warning: Could not open %s. No target can be selected." % starlist)
        return None
    else:
        for line in f:
            ls = line.split()
            names.append(ls[0])
            row = []
            # RA value in radians
            row.append(getRARad(ls[1], ls[2], ls[3]))
            # Dec value in radians
            row.append(getDECRad(ls[4], ls[5], ls[6]))
            # PM RA
            row.append(float(ls[8].split('=')[-1]))
            # PM Dec
            row.append(float(ls[9].split('=')[-1]))
            # V mag
            row.append(float(ls[10].split('=')[-1]))
            # Exposure time
            row.append(float(ls[11].split('=')[-1]))
            # Desired Counts
            row.append(float(ls[16].split('=')[-1]))
            # Filler not used here
            row.append(0.)
            row.append(0.)
            # Number of exposures
            row.append(int(ls[19].split('=')[-1]))

            star_table.append(row)

            # Save the scriptobs line for later
            lines.append(line)

            # Generate a pyEphem object for this target
            star = ephem.FixedBody()
            star._ra = ephem.hours(":".join([ls[1], ls[2], ls[3]]))
            star._dec = ephem.degrees(":".join([ls[4], ls[5], ls[6]]))
            stars.append(star)
            
    return names, np.array(star_table), lines, stars
            

def parseGoogledex(usr, pwd):

    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    try:
        f = open("/u/rjhanson/master/googledex.dat",'r')
    except IOError:
        apflog( "Starting Googledex parse")
        gs = gspread.login(usr,pwd)
        apflog( "Successfully logged in.")
        spreadsheet = gs.open("The Googledex")
        apflog( "Loaded Main Googledex")
        worksheet = spreadsheet.sheet1
        apflog( "Got spreadsheet")
        start = datetime.now()
        full_codex = worksheet.get_all_values()
        #time = (datetime.now() - start).total_seconds()
        #print "Loaded Values. Took {0:f} seconds.".format(time)
        f = open("/u/rjhanson/master/googledex.dat",'w')
        pickle.dump(full_codex, f)
        f.close()
    else:
        full_codex = pickle.load(f)
        f.close()
        
    col_names = full_codex[0]
    codex = full_codex[1:]

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFtexp", "APFpri", "APFcad", "APFnshots", "lastobs", \
                "B-V", "APF Desired Precision"
                ]

    try:
        idx = [col_names.index(v) for v in req_cols]
    except ValueError:
        errstr = v + " Not found in list"    
        apflog( errstr,level='Error' )
    names = []
    star_table = []
    do_flag = []
    stars = []
    # Build the star table to return to 
    for ls in codex:
        if float(ls[idx[11]]) < 0.5: continue
        row = []
        # Get the star name
        names.append(ls[idx[0]])
        # Get the RA
        row.append(getRARad(ls[idx[1]], ls[idx[2]], ls[idx[3]]))
        # Get the DEC
        row.append(getDECRad(ls[idx[4]], ls[idx[5]], ls[idx[6]]))
        for i in idx[7:11]: row.append(float(ls[i]))
        # For now use the old 1e9 count value
        row.append(1.e9)
        for i in idx[11:-1]: row.append(float(ls[i]))
        try:
            row.append(float(ls[idx[-1]]))
        except ValueError:
            row.append(1.5)
        
        star_table.append(row)
        do_flag.append("")
        star = ephem.FixedBody()
        star._ra = ephem.hours(":".join([ls[idx[1]], ls[idx[2]], ls[idx[3]]]))
        star._dec = ephem.degrees(":".join([ls[idx[4]], ls[idx[5]], ls[idx[6]]]))
        stars.append(star)

    return (names, np.array(star_table), do_flag, stars)

def getI2_K(unc):
    A = 4.47
    B = -1.58
    return 10 ** ( A + B*np.log10(unc))
    
def getI2_M(unc):
    A = 4.14
    B = -1.73
    return 10 ** ( A + B*np.log10(unc))
    
def getEXPMeter(i2, bv):
    delta = 4.52
    epsilon = -0.196
    squiggle = 0.262
    x = delta + epsilon*bv + squiggle*bv**2
    x += 0.05
    return i2 * 10**x

def getEXPMeter_Rate(v, bv, el, seeing, decker="W"):
    alpha = -0.908
    beta = 0.0852
    Const = -22.55
    if seeing == 0:
        apflog( "Warning: AVG_FWHM seems to be 0. Using 15 instead.",level="Warn")
        
        seeing = 15
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272
    
    VC = v - 2.5*np.log10(light)
    x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    return (10 ** x)


def getEXPTime(cnts, v, bv, el, seeing, decker="W"):
    alpha = -0.0311
    beta = 0.158
    Const = -11.95
    if seeing == 0:
        apflog( "Warning: AVG_FWHM seems to be 0. Using 15 instead.",level='Warn')
        seeing = 15
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272
    
    if el < 15.0:
        el = 15.0
        # bogus exposure time but the APF does not work this low anyway
    VC = v - 2.5*np.log10(light)
    x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    cnt_rate = 10**x
    time = 0
    time = cnts/cnt_rate
    return time


def getLST(date, longitude):
    """Take a datetime and longitude and calculate the Local Sidereal Time."""
    # Assumes date is a datetime object, and that the longitude is formatted as in PyEphem 

    ll = [float(v) for v in longitude.split(':')]
    if ll[0] > 0:
        sign = 1
    else:
        sign = -1
    ut = date.hour + date.minute/60. + date.second/3600.
    lng = ll[0] + sign*ll[1]/60. + sign*ll[2]/3600.
    d  = ephem.julian_date() - 2451545.0
    lst = 100.46 + 0.985647 * d + lng + 15*ut
    return lst % 360.


def getRADeg(hr, mn, sec):
    ra_hr = float(hr) + float(mn)/60. + float(sec)/3600.
    return ra_hr * 15

def getRARad(hr, mn, sec):
    ra_hours = float(hr) + float(mn)/60. + float(sec)/3600.
    return ra_hours * 15 * np.pi/180.0

def getDECDeg(deg, mn, sec, sign=None):
    if float(deg) < 0:
        sign = '-'
    x = abs(float(deg)) + float(mn)/60. + float(sec)/3600.
    if sign == '-':
        return x*-1
    else:
        return x

def getDECRad(deg, mn, sec, sign=None):
    if float(deg) < 0:
        sign = '-'
    x = abs(float(deg)) + float(mn)/60. + float(sec)/3600.
    x = x * np.pi/180.
    if sign == '-':
        return x*-1
    else:
        return x


def getElAz(ra, dec, lat, lng, time):
    """Given RA, DEC, Latitude, and a time, returns the corresponding elevation and azimuth angles
       Works with single values, or numpy arrays
       """
    lst = getLST(time, lng)
    ha = ((lst- np.degrees(ra)) % 360.) * np.pi/180.
    el = np.arcsin(np.sin(dec) * np.sin(lat) + \
                   np.cos(dec) * np.cos(lat) * np.cos(ha))
    az = np.arccos( (np.sin(dec) - np.sin(el)*np.sin(lat)) / \
                         (np.cos(el) * np.cos(lat)))
    return (np.degrees(el), np.degrees(az))

def makeScriptobsLine(name, row, do_flag, t):
    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = name + ' '
    # Add the RA as three elements, HR, MIN, SEC
    ra = np.degrees(row[DS_RA]) / 15.
    ret += str(int(ra)) + ' '
    ra_min = (ra % 1) * 60
    ret += str(int(ra_min)) + ' '
    ret += str(round((ra_min % 1) * 60., 1)) + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    dec = np.degrees(row[DS_DEC])
    ret += str(int(dec)) + ' '
    dec_min = abs((dec - int(dec)) * 60)
    ret += str(int(dec_min)) + ' '
    ret += str(round( (dec_min % 1) * 60, 1)) + ' '
    # Epoch
    ret += '2000 '
    # Proper motion RA and DEC
    ret += 'pmra=' + str(row[DS_PMRA]) + ' '
    ret += 'pmdec=' + str(row[DS_PMDEC]) + ' '
    # V Mag
    ret += 'vmag=' + str(row[DS_VMAG]) + ' '
    # T Exp
    if row[DS_EXPT] > MAX_EXPTIME:
        ret += "texp=%d " % (int(MAX_EXPTIME))
    else:
        ret += 'texp=' + str(int(row[DS_EXPT])) + ' '
    # I2
    ret += 'I2=Y '
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=' + str(t.hour) + ' '
    ret += 'utm=' + str(t.minute) + ' '
    # Exp Count
    ret += 'expcount=' + str(round(row[DS_COUNTS],1)) + ' '
    # Decker
    ret += 'decker=W '
    # do flag
    if do_flag:
        ret += 'do=Y '
    else:
        ret += 'do= '
    # Count
    ret += 'count=' + str(int(row[DS_NSHOTS])) 
    
    return ret

def getObserved(filename):
    obs = []
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename, level='Error')
        return obs
    else:
        for line in f:
            if line.strip()[0] == '#': continue
            obs.append(line.split()[0])
    return obs


def smartList(starlist, time, seeing, slowdown, az, el):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    dt = datetime.utcfromtimestamp(int(time))    

    observed = getObserved("/u/rjhanson/master/observed_targets")

    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = '30.'
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180.

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Parse the starlist
    try:
        sn, star_table, lines, stars = parseStarlist(starlist)
    except ValueError:
        # This will be raised when the starlist could not be parsed successfully.
        apflog( "No target could be selected because starlist could not be parsed.",level='Error')
        return None
    targNum = len(sn)

    # A place to hold star scores. Everything starts with a 0
    score = np.zeros(targNum) + 100

    # Minimum Brightness based on conditions
    #VMIN = 12 - 2.5 * np.log10(slowdown)
    VMIN = 18

    # Distance to stay away from the moon [Between 15 and 25 degrees]
    minMoonDist = ((moon.phase / 100.) * 10.) + 15  

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    elevation = np.zeros(targNum)

    done = True
    # Loop over each star 
    for i in range(targNum):
        # Have we already observed this target?
        if sn[i] in observed:
            score[i] = 5.
        else:
            done = False

        # If seeing is bad, only observe bright targets ( Large VMAG is dim star )       
        if star_table[i, DS_VMAG] > VMIN:
            score[i] = 0.

        # Is the star visible?
        # Version 2.0
        obs_length = star_table[i,DS_EXPT] * star_table[i,DS_NSHOTS] + 45 * star_table[i,DS_NSHOTS]

        # Finishing position
        apf_obs.date = dt + timedelta(seconds = obs_length)
        stars[i].compute(apf_obs)
        el = np.degrees(stars[i].alt)
        if el < 30. or el > 80.:
            score[i] = 0.
            continue
        apf_obs.date = dt
        stars[i].compute(apf_obs)
        el = np.degrees(stars[i].alt)
        if el < 30. or el > 80.:
            score[i] = 0.
            continue

        
        elevation[i] = el
        score[i] += 5
        

        # Is the star behind the moon?
        if moonDist[i] < minMoonDist:
            score[i] = 0.
            #print sn[i], "behind moon"

    if done:
        apflog( "All targets have been observed")
        return None

    # Prioritize higher elevation targets
    if max(elevation) > 0:
        score += elevation/max(elevation) * 100
    else:
        apflog( "No observable targets found.",level='Warn')
        return None

    idx = score.argsort()[-1]

    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV'] = 0.8
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NAME']   = sn[idx]
    res['SCORE']  = score[idx]
    res['SCRIPTOBS'] = lines[idx]
    return res
    
def format_time(total,hitthemall=False):
    if total < 300:
        return 300, min([np.ceil(300/total), 4])
    elif total > MAX_EXPTIME:
        if hitthemall:
            return MAX_EXPTIME, 1
        else:
            return MAX_EXPTIME, np.ceil(total/MAX_EXPTIME)

    else:
        return np.ceil(total), 1


def getNext(time, seeing, slowdown, wind=None, bstar=False, verbose=True):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time (Unix timestamp), seeing, slowdown, tele az, tele el
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
        """
    # Convert the unix timestamp into a python datetime
    dt = datetime.utcfromtimestamp(int(time))
    v = verbose

    # List of targets already observed
    observed = getObserved('/u/rjhanson/master/observed_targets')
    if observed == []:
        apflog( "getObserved is empty, setting bstar to true")
        bstar = True

    #
    ##
    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###
    ##
    #
         
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = '30.'
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180. 

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians
    with open("/u/rjhanson/master/apf_google_login",'r') as f:
        usr, pwd = pickle.load(f)
    sn, star_table, do_flag, stars = parseGoogledex(usr, pwd)
    targNum = len(sn)

    
    # Note which of these are B-Stars for later.
    bstars = []
    for n in sn:
        if 'HR' in n:
            bstars.append(True)
        else:
            bstars.append(False)


    # A place to hold star scores. Everything starts with a 0
    score = np.zeros(targNum)

    # Minimum Brightness based on conditions
    #VMIN = 12 - 2.5 * np.log10(slowdown)
    #print "VMIN = ", VMIN

    # Distance to stay away from the moon [Between 15 and 25 degrees]
    minMoonDist = ((moon.phase / 100.) * 10.) + 15  

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    
    # Loop over each star 
    ai2counts = []
    for i in range(targNum):

        #
        ##
        ###
        # Calculate the expected exposure time
        # Use this to get a better finishing time
        # Number of exposures needs to be NShots = Math.Ceiling(600 / expTime )
        ###
        ##
        #
        stars[i].compute(apf_obs)
        el = np.degrees(stars[i].alt)

        # Turn desired precision into counts/ exp time
        if bstars[i] == False:
            precision = star_table[i, DS_ERR]
            if precision < 1.5 and star_table[i, DS_APFPRI] < 10:
                precision = 1.5
            if star_table[i, DS_BV] > 1.2:
                i2counts = getI2_M(precision)
            else:
                i2counts = getI2_K(precision)
#            i2counts = getI2(precision)

            exp_time   = getEXPTime(i2counts, star_table[i, DS_VMAG], star_table[i, DS_BV], el, seeing)
            # the above is the desired value
            # now, for real-time control, we need to calculate the expected
            # number of counts on the exposure meter itself
            exp_time *= slowdown
            exp_counts = getEXPMeter(i2counts, star_table[i, DS_BV])

#            star_table[i, DS_COUNTS] = exp_counts
            #star_table[i, DS_EXPT] = 900
            star_table[i, DS_EXPT], star_table[i, DS_NSHOTS] = format_time(exp_time)
            star_table[i, DS_COUNTS] = exp_counts / star_table[i, DS_NSHOTS]

            #print sn[i], star_table[i, DS_EXPT], star_table[i, DS_NSHOTS]
        else:
            star_table[i, DS_EXPT] = 200
            star_table[i, DS_NSHOTS] = 2  
            i2counts = 0
        
        ai2counts.append(i2counts)        
        # --- Reasons to prioritize a target ---

        # APF Priority
        apf_pri = star_table[i, DS_APFPRI]
        if apf_pri >= 10:
            score[i] += apf_pri * 10
        elif apf_pri == 9.9:
            score[i] += 80
        elif apf_pri == 9.8:
            score[i] += 60
        else:
            score[i] += 40


        # --- Things to disqualify a star ---

        # Have we already observed this target?
        if sn[i] in observed:
            score[i] = 0.

        # If seeing is bad, only observe bright targets ( Large VMAG is dim star )       
        #if star_table[i, DS_VMAG] > VMIN:
        #    score[i] = 0.
        
        # If it would take too much time to observe this target
        #  (# exp) * (exp time) is longer than an hour
        if star_table[i, DS_EXPT] * star_table[i, DS_NSHOTS] > 1e8:
            outstr = "Insane values in getNext(): i2counts = %.1f, slowdown = %.2f, seeing=%.2f, el=%.1f, vmag=%.2f b-v=%.2f" % (i2counts,slowdown,seeing,el,star_table[i, DS_VMAG], star_table[i, DS_BV])
            apflog(outstr)
        if star_table[i, DS_EXPT] * star_table[i, DS_NSHOTS] > 2.0 * 60 * 60:
            if apf_pri > 9.9:
                if star_table[i, DS_EXPT] * star_table[i, DS_NSHOTS] > 3.0 * 60 * 60: # 3 hours for a priority 10 target
                    score[i] = 0
                    if v: apflog("getNext(): exposure %.1f too long for %s" % ( exp_time, sn[i] ))
                    continue
            else:
                score[i] = 0
                if v: apflog("getNext(): exposure %.1f too long for %s" % ( exp_time, sn[i]))
                continue

        # Is the star visible?
        # Version 2.0
        obs_length = star_table[i,DS_EXPT] * star_table[i,DS_NSHOTS] + 45 * star_table[i,DS_NSHOTS]
        # print obs_length, star_table[i,DS_NSHOTS], star_table[i,DS_EXPT], star_table[i,DS_COUNTS], sn[i]
        try:
            apf_obs.date = dt + timedelta(seconds = obs_length)
        except:
            apf_obs.date = dt
        if el < 30. or el > 85.:
            if v: apflog("getNext(): %s not visible with elevation %f" % ( sn[i], el))
            score[i] = 0.
            continue
        apf_obs.date = dt
        stars[i].compute(apf_obs)
        el = np.degrees(stars[i].alt)
        if el < 30. or el > 85.:
            if v: apflog("getNext(): %s not visible with elevation %f" % ( sn[i], el))
            score[i] = 0.
            continue
        
        score[i] += 5
        

        # Is the star behind the moon?
        if moonDist[i] < minMoonDist:
            if v: apflog("getNext(): %s is behind the moon - distance = %f" % ( sn[i], moonDist[i] ) )
            score[i] = 0.
            #print sn[i], "behind moon"

        # Do we need a B Star?
        if bstar:
            if not bstars[i]:
                score[i] = 0.
                continue
        if not bstar:
            if v: apflog("getNext(): %s is observable with exposure time %.1f" % ( sn[i], exp_time ) )


    # Index of our target
    idx = score.argsort()[-1]
#    print sn[idx], score[idx]
    if score[idx] < 5:
        return None
    
    if bstars[idx] != bstar:
        apflog("bstars[idx] = "+ str(bstars[idx]) + " bstar = " + str(bstar),level='Warn')
        apflog( "getNext() returning None",level='Warn')
        return None
    
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV']     = star_table[idx, DS_BV]
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NAME']   = sn[idx]
    res['SCORE']  = score[idx]
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['I2COUNTS'] = ai2counts[idx]
    res['SCRIPTOBS'] = makeScriptobsLine(sn[idx], star_table[idx,:], do_flag[idx], dt)
    return res




if __name__ == '__main__':
    # For some test input what would the best target be?
    
    result = getNext(time.time(), 13.99, 4., None, 75, wind=(8.2, 273.9), bstar=False)
    #result = smartList("observed_targets.3", time.time(), 13.5, 2.4, 120, 60)
    
    if result is None:
        print "Get None target"
    else:
        for k in result:
            print k, result[k]

    print "Done"








