#! /usr/bin/python
# UCSCScheduler.py

import calendar
from datetime import datetime, timedelta
import ephem
import gspread
import json
from oauth2client.client import SignedJwtAssertionCredentials
from ExposureCalculations import getI2_M, getI2_K, getEXPMeter, getEXPMeter, getEXPTime

import numpy as np
import os
import pickle
import sys
import time
from apflog import *

# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX = 2 * 60 * 60 # 2 hours
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

# Maximum single exposure time in seconds
MAX_EXPTIME = 1200.


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

slit_size = {'M': (9.52, 74.07), 'W': (9.26, 27.78), 'N': (4.63, 74.07), 'B': (18.52, 74.07),
             'S': (6.94, 74.07), 'P': (1, 1) }

def parseStarlist(starlist):
    """ Parse a scriptobs-compatible starlist for the scheduler."""
    names = []
    lines = []
    stars = []
    star_table = []
    try:
        f = open(starlist,'r')
    except IOError:
        apflog("Warning: Could not open %s. No target can be selected." % starlist,echo=True)
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

def get_spreadsheet(filen='UCSC Dynamic Scheduler-5b98d1283a95.json'):

    # this downloads the googledex from the Google Drive
    # the certificate must be available
    # these certificates are generated through the Google Developer Interface
    # the developer must select the correct API for access

    # the certificate has an email associated with it, that email must
    # have the document shared with it to allow access 

    json_key = json.load(open(filen))
    scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gs = gspread.authorize(credentials)

    apflog("Successfully logged in.", echo=True)
    spreadsheet = gs.open("The Googledex")
    apflog("Loaded Main Googledex",echo=True)
    worksheet = spreadsheet.sheet1
    apflog("Got spreadsheet", echo=True)
    start = datetime.now()

    return worksheet

def parseGoogledex():

    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    try:
        f = open("./googledex.dat",'r')
    except IOError:
        apflog( "Starting Googledex parse")
        worksheet = get_spreadsheet()
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
        apflog("%s Not found in list" % (v) , level="warn")
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
    
def update_googledex_lastobs(filename, time=None):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
    """
    names, times = getObserved(filename)
    
    if time is None:
        time = datetime.utcnow()
    

    ws = get_spreadsheet()
    vals = ws.get_all_values()

    col = vals[0].index("lastobs") + 1
    
    for i, v in enumerate(vals):
        # Did we observe this target tonight?
        if v[0] in names:
            # We observed this target, so update the cell in the worksheet
            # update_cell(row, col, val) - col and row are 1 indexed
            hr, min = times[names.index(v[0])]
            
            t = datetime(time.year, time.month, time.day, hr, min)
            jd = float(ephem.julian_date(t))
            ws.update_cell(i+1, col, round(jd, 1) )
    apflog( "Updated Googledex")

def update_local_googledex(googledex_file="/u/rjhanson/master/googledex.dat", observed_file="/u/rjhanson/master/observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time. 
    """
    names, times = getObserved(observed_file)

    time = datetime.utcnow()

    try:
        g = open(googledex_file, 'r')
    except IOError:
        apflog("googledex file did not exist, so can't be updated")
        return

    full_codex = pickle.load(g)
    g.close()

    codex_cols = full_codex[0]

    starNameIdx = codex_cols.index("Star Name")
    lastObsIdx = codex_cols.index("lastobs")
    
    for i in range(1, len(full_codex)):
        row = full_codex[i]
        if row[starNameIdx] in names:
            # We have observed this star, so lets update the last obs field
            hr, min = times[names.index(row[starNameIdx])]
            
            t = datetime(time.year, time.month, time.day, hr, min)
            # This keeps the JD precision to one decimal point. There is no real reason for this other than
            # the googledex currently only stores the lastObs field to one decimal precision. Legacy styles FTW.
            jd = round(float(ephem.julian_date(t)), 1) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)))
            row[lastObsIdx] = str(jd)
            full_codex[i] = row

    with open(googledex_file, 'w') as f:
        pickle.dump(full_codex, f)
            



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
        ret += 'texp=%d ' % (int(MAX_EXPTIME))
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
    times = []
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename,level="warn")
        return obs, times
    else: 
        for line in f:
            if line.strip()[0] == '#' or line.strip() == "": continue
            ls = line.split()
            obs.append(ls[0])
            times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
            
    return obs, times
	
def calculate_ucsc_exposure_time( vmag, precision, elevation, seeing, bv=None, priority=None, decker="W"):
	"""
		Args:
			vmag: List of V magnitudes for targets
			precision: List of desired RV precisions for targets
			elevation: The last target elevation ( Used for the last airmass )
			seeing: Current FWHM
		Optional Args:
			priority: List of priorities for targets
			bv: List of B-V magnitudes for targets
			decker: The current decker being used, defaults to the standard UCSC 'W' decker
		Returns: tuple( exp_time, exp_counts )
			exp_time: List of exposure times in seconds for the given targets
			exp_counts: List of expected exposure meter counts for given target
	"""
	vmag = np.array(vmag)
	precision = np.array(precision)
	
	if bv is None:
	    bv = np.array( [1.0] * len(vmag) )
	else:
	    bv = np.array(bv)
		
	# Now lets calculate the exposure times
	
	# Desired I2 counts for precision
	i2counts = getI2_K(precision)
    mstars = np.where(bv > 1.2)
    if len(mstars) > 0:
        i2counts[mstars] = getI2_M(precision[mstars])
	
	# Exposure Meter counts to reach desired I2 counts
	exp_counts = getEXPMeter(i2counts, bv)
	
	# Exposure time to reach desired I2 counts
	exp_time = getEXPTime(i2counts, vmag, bv, elevation, seeing, decker=decker)
	
	return exp_time, exp_counts

def is_visible(stars, observer, obs_len, min_el, max_el):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). This is the time frame for which visibility is checked
            min_el: The minimum body elevation to be visible ( degrees )
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    prev_horizon = observer.horizon

    ret = []

    observer.horizon = min_el
    # Now loop over each body to check visibility
    for s, dt in zip(stars, obs_len):
        s.compute(observer)

        # Is the target visible now?
        cur_el = np.degrees(s.alt)
        if cur_el < min_el or cur_el > max_el:
            ret.append(False)
            continue
		
        # Does the target remain visible through the observation?
        # The next setting/rising functions throw an exception if the body never sets or rises
        # ex. circumpolar 
        try:
            next_set = observer.next_setting(s)
        except:
            # If it never sets, no need to worry.
            pass
        else:
            # Making the assumption that next_set is a datetime object. Might not be the case
            if next_set < dt:
                # The object will set before the observation finishes
                ret.append(False)
                continue

        observer.horizon = max_el
        s.compute(observer)

        try:
            next_rise = observer.next_rising(s)
        except:
            # If the body never rises above the max limit no problem
            pass
        else:
            if next_rise < dt:
                # The object rises above the max el before the observation finishes
                ret.append(False)
                continue
		
        # Everything seems to be fine, so the target is visible!
        ret.append(True)
	
    observer.horizon = prev_horizon
    return ret


def smartList(starlist, time, seeing, slowdown, az, el):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    dt = datetime.utcfromtimestamp(int(time))    

    observed, _ = getObserved("/u/rjhanson/master/observed_targets")

    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(TARGET_ELEVATION_MIN)
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
        apflog( "No target could be selected because starlist could not be parsed.", level="warn")
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
        apflog( "All targets have been observed",level="warn")
        return None

    # Prioritize higher elevation targets
    if max(elevation) > 0:
        score += elevation/max(elevation) * 100
    else:
        print "No observable targets found."
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
    
def format_time(total, hitthemall=False):
    total = np.array(total)
    idx = np.where(total < 300, True, False)
    jdx = np.where(total > MAX_EXPTIME, True, False)
    #kdx = not idx and not jdx
    kdx = np.logical_not(np.logical_or(idx, jdx))

    times = np.zeros(len(total))
    exps  = np.zeros(len(total))

    times[idx] = 300
    times[jdx] = MAX_EXPTIME
    times[kdx] = np.ceil(total[kdx])

    exps[idx] = [ min([np.ceil(300/t), 4]) for t in total[idx] ]
    
    if hitthemall:
        exps[jdx] = 1
    else:
        exps[jdx] = np.ceil(total[jdx]/MAX_EXPTIME)
    exps[kdx] = 1

    return times, exps


def getNext(time, seeing, slowdown, bstar=False, verbose=False):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time (Unix timestamp), seeing, slowdown, tele az, tele el
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """
    # Convert the unix timestamp into a python datetime
    dt = datetime.utcfromtimestamp(int(time))
    if verbose:
        print "getNext(): Finding target for time ", dt

    update_local_googledex()

    # List of targets already observed
    observed, _ = getObserved('/u/rjhanson/master/observed_targets')
    if observed == []:
        if verbose:
            print "getNext(): getObserved is empty, setting bstar to true"
        bstar = True

    ###
    # Need to update the googledex with the lastObserved date for observed targets
    # Scriptobs line uth utm can be used for this
    # Need to convert a uth and utm to a JD quickly.
    # timedelta = now - uth,utm : minus current JD?
    ###
         
    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(TARGET_ELEVATION_MIN)
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = (37 + 20/60. + 33.1/3600.) * np.pi/180. 

    # Calculate the moon's location
    moon = ephem.Moon()
    moon.compute(apf_obs)

    # Parse the Googledex
    # Note -- RA and Dec are returned in Radians
    if verbose:
        print("getNext(): Parsing the Googledex...")
    sn, star_table, do_flag, stars = parseGoogledex()
    sn = np.array(sn)
    targNum = len(sn)

    # Note which of these are B-Stars for later.
    bstars = np.array([ True if 'HR' in n else False for n in sn ], dtype=bool)


    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MAX
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  


    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    available = np.ones(targNum, dtype=bool)

    # Is the target behind the moon?
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check
    
    
    star_elevations = []
    if len(stars) != len(available): print "Error, arrays didn't match"
    for idx in range(len(stars)):
        stars[idx].compute(apf_obs)
        
        star_el = np.degrees(stars[idx].alt)
        
        if star_el < TARGET_ELEVATION_MIN or star_el > TARGET_ELEVATION_MAX:
            # This star is outside our visible range
            available[idx] = False

        star_elevations.append(star_el)
    star_elevations = np.array(star_elevations)

    print "Pre loop elevations"
    print star_elevations[available]
    print len(star_elevations[available])

    # We just need a B star, so restrict our math to those
    if bstar:
        available = available & bstars
        f = available
        vis = is_visible([s for s,_ in zip(stars,f) if _ ], apf_obs, [400]*len(bstars[f]), TARGET_ELEVATION_MIN, TARGET_ELEVATION_MAX)
		
        available[f] = available[f] & vis
		
        star_table[available, DS_COUNTS] = 1e9
        star_table[available, DS_EXPT] = 200
        star_table[available, DS_NSHOTS] = 2
		
	# Just need a normal star for observing
    else:
        # Available and not a BStar
        available = np.logical_and(available, np.logical_not(bstars))
        # has the star been observed 
        done = [ True if n in observed else False for n in sn ]
        available = available & np.logical_not(done) # Available and not observed

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function
        f = available
        
        print len(f), len(star_elevations)
        
        exp_times, exp_counts = calculate_ucsc_exposure_time( star_table[f,DS_VMAG], \
                                        star_table[f,DS_ERR], star_elevations[f], seeing, \
                                        star_table[f,DS_BV], star_table[f,DS_APFPRI])
        
        exp_times = exp_times * slowdown

        star_table[f, DS_COUNTS] = exp_counts
        star_table[f, DS_EXPT], star_table[f, DS_NSHOTS] = format_time(exp_times)

        # Is the exposure time too long?
        time_check = np.where( exp_times < TARGET_EXPOSURE_TIME_MAX, True, False)
        print time_check
        available[f] = available[f] & time_check
        f = available

        # Is the star currently visible?
        vis = is_visible([s for s,_ in zip(stars,f) if _ ], apf_obs, exp_times, TARGET_ELEVATION_MIN, TARGET_ELEVATION_MAX)
        if vis != []:
            available[f] = available[f] & vis

    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        print "Couldn't find any suitable targets!"
        return None

    cadence_check = (ephem.julian_date(dt) - star_table[:, DS_LAST]) / star_table[:, DS_CAD]

    pri = max(star_table[available, DS_APFPRI])
    sort_i = np.where(star_table[available, DS_APFPRI] == pri, True, False)

    print sn[available][sort_i]
    print star_table[available, DS_APFPRI][sort_i]
 
    sort_j = cadence_check[available][sort_i].argsort()[::-1]
    print cadence_check[available][sort_i][sort_j][0]

    t_n = sn[available][sort_i][sort_j][0]

    print "Star elevations"
    print star_elevations[available][sort_i][sort_j]


    print t_n

    idx, = np.where(sn == t_n)
    print idx[0]
    idx = idx[0]

    stars[idx].compute(apf_obs)
    
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
    res['SCORE']  = 0.0
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['SCRIPTOBS'] = makeScriptobsLine(sn[idx], star_table[idx,:], do_flag[idx], dt)
    return res




if __name__ == '__main__':
    # For some test input what would the best target be?
    
    result = getNext(time.time(), 13.99, 1.8, bstar=False, verbose=True)
    #result = smartList("observed_targets.3", time.time(), 13.5, 2.4, 120, 60)

    print "testing googledex updater"
    update_googledex_lastobs('observed_targets')

    if result is None:
        print "Get None target"
    else:
        for k in result:
            print k, result[k]

    print "Done"
