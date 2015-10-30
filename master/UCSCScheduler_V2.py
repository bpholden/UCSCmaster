# UCSCScheduler_V2.py

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
try:
    from apflog import *
    import ktl
except:
    from fake_apflog import *
import re

# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX =  60 * 60 # 1 hour
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

# Maximum single exposure time in seconds
MAX_EXPTIME = 1200.
MIN_EXPTIME = 300.
MIN_SHUTTERTIME = 60
MAX_I2 = 40000
MAX_EXPMETER = 1e9

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

SLOWDOWN_MIN = 0.4
SLOWDOWN_MAX = 10.0

###
# arGGGGGG!!
###

last_objs_attempted = []

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
            if not re.search("\A\#",line):
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

def get_spreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):

    # this downloads the googledex from the Google Drive
    # the certificate must be available
    # these certificates are generated through the Google Developer Interface
    # the developer must select the correct API for access

    # the certificate has an email associated with it, that email must
    # have the document shared with it to allow access 

    certificate_path = os.path.dirname(__file__)
    
    json_key = json.load(open(os.path.join(certificate_path, certificate)))
    scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gs = gspread.authorize(credentials)

    apflog("Successfully logged in.", echo=True)
    spreadsheet = gs.open(sheetn)
    apflog("Loaded Main %s" % (sheetn),echo=True)
    worksheet = spreadsheet.sheet1
    apflog("Got spreadsheet", echo=True)

    return worksheet

def findColumns(col_names,req_cols):
    
    idx = []
    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) , level="Alert",echo=True)

    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in didx.keys():
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Error",echo=True)

    return didx

def parseGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat"):

    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    try:
        f = open(os.path.join(os.getcwd(),outfn),'r')
    except IOError:
        apflog( "Starting Googledex parse",echo=True)
        worksheet = get_spreadsheet(sheetn=sheetn,certificate=certificate)
        full_codex = worksheet.get_all_values()
        #time = (datetime.now() - start).total_seconds()
        #print "Loaded Values. Took {0:f} seconds.".format(time)
        f = open(os.path.join(os.getcwd(),"googledex.dat"),'w')
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
                "B-V", "APF Desired Precision", "Close Companion"
                ]
    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    do_flag = []
    stars = []
    # Build the star table to return to 
    for ls in codex:
        if ls[0] == '':
            continue
        if float(ls[didx["APFpri"]]) < 0.5: continue
        row = []
        # Get the star name
        names.append(ls[didx["Star Name"]])
        # Get the RA
        row.append(getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]))
        # Get the DEC
        row.append(getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]))
        for i in ("pmRA", "pmDEC", "Vmag","APFtexp"):
            try:
                row.append(float(ls[didx[i]]))
            except ValueError:
                row.append(0.0)
        # For now use the old 1e9 count value
        row.append(1.e9)
        for i in ("APFpri", "APFcad", "APFnshots", "lastobs", "B-V", "APF Desired Precision" ):
            try:
                row.append(float(ls[didx[i]]))
            except ValueError:
                row.append(0.0)

        match = re.search("\A(y|Y)",ls[didx["Close Companion"]])
        if match:
            do_flag.append("y")
        else:
            do_flag.append("")
        
        star_table.append(row)
        star = ephem.FixedBody()
        star._ra = ephem.hours(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]]))
        star._dec = ephem.degrees(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]]))
        stars.append(star)

    return (names, np.array(star_table), do_flag, stars)
    
def update_googledex_lastobs(filename, sheetn="The Googledex",time=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
    """
    names, times = getObserved(filename)
    if len(names) == 0:
        return
    if time is None:
        time = datetime.utcnow()
    

    ws = get_spreadsheet(sheetn=sheetn,certificate=certificate)
    vals = ws.get_all_values()

    col = vals[0].index("lastobs") 
    
    for i, v in enumerate(vals):
        # Did we observe this target tonight?
        if v[0] in names:
            # We observed this target, so update the cell in the worksheet
            # update_cell(row, col, val) - col and row are 1 indexed
            otime = times[names.index(v[0])]
            if isinstance(otime,float):
                t = datetime.fromtimestamp(otime)
            else:
                hr, min = otime
                t = datetime(time.year, time.month, time.day, hr, min)
            jd = float(ephem.julian_date(t))
            if jd > float(v[col]):
                ws.update_cell(i+1, col+1, round(jd, 2) )
    apflog( "Updated Googledex",echo=True)

def update_local_googledex(time,googledex_file="googledex.dat", observed_file="observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time. 
    """
    names, times = getObserved(observed_file)

    try:
        g = open(googledex_file, 'r')
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
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
            otime = times[names.index(row[starNameIdx])]
            if isinstance(otime,float):
                t = datetime.fromtimestamp(otime)
            else:
                time = datetime.utcnow()
                hr, min = otime
                t = datetime(time.year, time.month, time.day, hr, min)

            # This keeps the JD precision to one decimal point. There is no real reason for this other than
            # the googledex currently only stores the lastObs field to one decimal precision. Legacy styles FTW.
            jd = round(float(ephem.julian_date(t)), 1) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)),echo=True)
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

def makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y"):

    focval = 0
    if row[DS_APFPRI] > 9.9:
        focval = 2
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
    ret += 'I2=%s ' % (I2)
    # lamp
    ret += 'lamp=none '
    # start time
    ret += 'uth=' + str(t.hour) + ' '
    ret += 'utm=' + str(t.minute) + ' '
    # Exp Count
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
        
    return ret

def getObserved(filename):
    obs = []
    times = []
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename,level="warn",echo=True)
        return obs, times
    else: 
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '#' or line == "":
                    pass
                else:
                    ls = line.split()
                    obs.append(ls[0])
                    if len(ls) > 15:
                        times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
                    else:
                        times.append(float(ls[1]))
            
    obs.reverse()
    times.reverse()
    return obs, times
	
def calculate_ucsc_exposure_time(vmag, precision, elevation, seeing, bmv, decker="W"):
    vmag = np.array(vmag)
    precision = np.array(precision)
    bmv = np.array(bmv)
    precision = np.array(precision)
		
	# Now lets calculate the exposure times
	
	# Desired I2 counts for precision
    i2counts = getI2_K(precision)
    mstars = np.where(bmv > 1.2)
    if len(mstars) > 0:
        i2counts[mstars] = getI2_M(precision[mstars])
	
	# Exposure Meter counts to reach desired I2 counts
	exp_counts = getEXPMeter(i2counts, bmv)
#	exp_counts = 1e9
	# Exposure time to reach desired I2 counts
	exp_time = getEXPTime(i2counts, vmag, bmv, elevation, seeing, decker=decker)
	
	return exp_time, exp_counts, i2counts

def calc_elevations(stars, observer):
    els = []
    for s in stars:
        observer.date = ephem.Date(observer.date)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        els.append(cur_el)
    return np.array(els)

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
    cdate = observer.date
    ret = []
    fin_elevations = []
    start_elevations = []    
    observer.horizon = min_el
    # Now loop over each body to check visibility
    for s, dt in zip(stars, obs_len):
        # Is the target visible at the end of the observations?
        observer.date = ephem.Date(cdate + dt/86400.)
        s.compute(observer)
        fin_el = np.degrees(s.alt)
        fin_elevations.append(fin_el)

        # Is the target visible now?
        observer.date = ephem.Date(cdate)
        s.compute(observer)
        cur_el = np.degrees(s.alt)
        start_elevations.append(cur_el)
        
        if fin_el < min_el or fin_el > max_el:
            ret.append(False)
            continue

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
    return ret, np.array(start_elevations), np.array(fin_elevations)


def smartList(starlist, time, seeing, slowdown):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    # Convert the unix timestamp into a python datetime

    # punt
    dt = datetime.utcnow()

    if type(time) == float:
        dt = datetime.utcfromtimestamp(int(time))
    elif type(time) == datetime:
        dt = time
    elif type(time) == ephem.Date:
        dt = time.datetime()

    observed, _ = getObserved(os.path.join(os.getcwd(),"observed_targets"))

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
        apflog( "No target could be selected because starlist could not be parsed.", level="warn",echo=True)
        return None
    targNum = len(sn)

    # Minimum Brightness based on conditions
    VMAX = 14

    # Distance to stay away from the moon [Between 15 and 25 degrees]

    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  

    available = np.ones(targNum, dtype=bool)
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check

    # If seeing is bad, only observe bright targets ( Large VMAG is dim star )       
    brightenough = np.where(star_table[:, DS_VMAG] < VMAX,True,False)
    available = available & brightenough

    obs_length = star_table[:,DS_EXPT] * star_table[:,DS_NSHOTS] + 45 * (star_table[:,DS_NSHOTS]-1)
    vis, star_elevations, fin_els = is_visible(stars,apf_obs,obs_length,TARGET_ELEVATION_MIN, TARGET_ELEVATION_MAX)
    available = available & vis
        
    done = [ True if n in observed else False for n in sn ]
    availableandnotdone = available & np.logical_not(done)
    
    if not any(availableandnotdone):
        apflog( "All visible targets have been observed",level="warn",echo=True)
        (good,) = np.where(available)
    else:
        (good,) = np.where(availableandnotdone)

    delta = fin_els[good] - star_elevations[good]
    neg   = np.where(delta < 0)
    pos   = np.where(delta >= 0)
    inv_els = fin_els[good]
    inv_els[neg]  = fin_els[good[neg]] - 90
    inv_els[pos] = 90 - fin_els[good[pos]]
    sort_fin_els_idx = np.argsort(inv_els)
    idx = good[sort_fin_els_idx[0]]
                  
    res = dict()

    res['RA']     = stars[idx].a_ra
    res['DEC']    = stars[idx].a_dec
    res['PM_RA']  = star_table[idx, DS_PMRA]
    res['PM_DEC'] = star_table[idx, DS_PMDEC]
    res['VMAG']   = star_table[idx, DS_VMAG]
    res['BV'] = 0.6
    res['PRI'] = 10.
    res['SCORE'] = 1.0
    res['COUNTS'] = star_table[idx, DS_COUNTS]
    res['EXP_TIME'] = star_table[idx, DS_EXPT]
    res['NEXP'] = star_table[idx, DS_NSHOTS]
    res['NAME']   = sn[idx]
    res['SCRIPTOBS'] = lines[idx]
    return res

def format_expmeter(exp_counts, nexp):
    
    exp_counts *= 1.1 
    long_idx = np.where((exp_counts/nexp) > MAX_EXPMETER, True, False)
    exp_counts[long_idx] = MAX_EXPMETER
    nexp[long_idx] = np.ceil((exp_counts[long_idx]/MAX_EXPMETER) + 1)
    return exp_counts, nexp

def format_time(total, i2counts, hitthemall=False):
    total = np.array(total)
    times = np.zeros(len(total))
    exps  = np.zeros(len(total))

    short_idx = np.where(total < MIN_EXPTIME, True, False)
    times[short_idx] = MIN_EXPTIME  # pad out to make it more likely exposure meter threshold sets actual limit
    exps[short_idx] = [ (np.ceil(MIN_EXPTIME/(t+40)) + 1) for t in total[short_idx] ] 

    middle_idx = np.where((total > MIN_EXPTIME ) &(total < MAX_EXPTIME))
    times[middle_idx] = MAX_EXPTIME # pad out to make it more likely exposure meter threshold sets actual limit
    exps[middle_idx] = 1

    max_idx = np.where(total > MAX_EXPTIME, True, False)
    if hitthemall:
        exps[max_idx] = 1
    else:
        exps[max_idx] = np.ceil(total[max_idx]/MAX_EXPTIME)
    times[max_idx] = MAX_EXPTIME



    bright_idx = np.where((i2counts > MAX_I2) & (exps == 1), True, False)
    exps[bright_idx] = [ np.ceil(i/MAX_I2) for i in i2counts[bright_idx] ]
    times[bright_idx] = np.ceil(total[bright_idx]/exps[bright_idx])

    return times, exps


def getNext(time, seeing, slowdown, bstar=False, verbose=False,sheetn="The Googledex"):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """
    # Convert the unix timestamp into a python datetime
    if type(time) == float:
        dt = datetime.utcfromtimestamp(int(time))
    elif type(time) == datetime:
        dt = time
    elif type(time) == ephem.Date:
        dt = time.datetime()
    else:
        dt = datetime.utcnow()
        # punt
        
    if verbose:
        apflog( "getNext(): Finding target for time %s" % (dt),echo=True)

    if slowdown > SLOWDOWN_MAX:
        apflog( "getNext(): Slowndown value of %f exceeds maximum of %f at time %s" % (slowdown,SLOWDOWN_MAX,dt),echo=True)
        return None


    try:
        apfguide = ktl.Service('apfguide')
        ptime = apfguide['midptfin'].read(binary=True)
    except:
        ptime = datetime.utcnow()
    update_local_googledex(ptime,googledex_file=os.path.join(os.getcwd(),"googledex.dat"), observed_file=os.path.join(os.getcwd(),"observed_targets"))

    # List of targets already observed
    observed, _ = getObserved(os.path.join(os.getcwd(),'observed_targets'))
    global last_objs_attempted
    try:
        lastline = ktl.read("apftask","SCRIPTOBS_LINE")
        lastobj = lastline.split()[0]
    except:
        lastobj = None

    if lastobj:
        if lastobj not in observed and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)
            if len(last_objs_attempted) > 5:
                last_objs_attempted = []
                return None
    if len (last_objs_attempted) > 0:
        observed += last_objs_attempted
#    if observed == []:
#        if verbose:
#            apflog( "getNext(): getObserved is empty, setting bstar to true",echo=True)
#        bstar = True

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
        apflog("getNext(): Parsing the Googledex...",echo=True)
    sn, star_table, do_flag, stars = parseGoogledex(sheetn=sheetn)
    sn = np.array(sn)
    targNum = len(sn)

    # Note which of these are B-Stars for later.
    bstars = np.array([ True if 'HR' in n else False for n in sn ], dtype=bool)


    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  


    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    available = np.ones(targNum, dtype=bool)
    totexptimes = np.zeros(targNum, dtype=float)
    cur_elevations = np.zeros(targNum, dtype=float)
    i2cnts = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check
    
    

#    apflog( "Pre loop elevations", echo=True)
#    elstr = "Stars els not behind moon: %s %d" % ( star_elevations[available],len(star_elevations[available]))
#    apflog(elstr, echo=True)

    # We just need a B star, so restrict our math to those
    if bstar:
        available = available & bstars
        f = available
        vis,star_elevations,fin_star_elevations = is_visible([s for s,_ in zip(stars,f) if _ ], apf_obs, [400]*len(bstars[f]), TARGET_ELEVATION_MIN, TARGET_ELEVATION_MAX)
        
        available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]
       		
        star_table[available, DS_COUNTS] = 1e9
        star_table[available, DS_EXPT] = 900
        star_table[available, DS_NSHOTS] = 2
        totexptimes[available] = 400

    # Just need a normal star for observing
    else:
        # Available and not a BStar
        available = np.logical_and(available, np.logical_not(bstars))
        # has the star been observed 
        done = [ True if n in observed else False for n in sn ]
        available = available & np.logical_not(done) # Available and not observed

        mag_limit = 12
        if slowdown > 1.0:
            mag_limit = -2.5*np.log10(slowdown/SLOWDOWN_MIN) + 10
        bright_enough = np.where(star_table[:,DS_VMAG] < mag_limit, True, False)
        available = available & bright_enough
        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function
        f = available
        star_elevations=calc_elevations([s for s,_ in zip(stars,f) if _ ],apf_obs)
        exp_times, exp_counts, i2counts = calculate_ucsc_exposure_time( star_table[f,DS_VMAG], \
                                            star_table[f,DS_ERR], star_elevations, seeing, \
                                            star_table[f,DS_BV])
        
        exp_times = exp_times * slowdown
        totexptimes[f] += exp_times
        i2cnts[f] += i2counts
        star_table[f, DS_EXPT], star_table[f, DS_NSHOTS] = format_time(exp_times,i2counts)
        #        exp_counts /= star_table[f, DS_NSHOTS]
        star_table[f, DS_COUNTS], star_table[f, DS_NSHOTS] = format_expmeter(exp_counts,star_table[f, DS_NSHOTS])
        #        star_table[f, DS_COUNTS] = 1.1*exp_counts / star_table[f, DS_NSHOTS]

        # Is the exposure time too long?
        time_check = np.where( exp_times < TARGET_EXPOSURE_TIME_MAX, True, False)
        
        available[f] = available[f] & time_check
        f = available

        # Is the star currently visible?
        vis,star_elevations,fin_star_elevations = is_visible([s for s,_ in zip(stars,f) if _ ], apf_obs, exp_times, TARGET_ELEVATION_MIN, TARGET_ELEVATION_MAX)
        if vis != []:
            available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]

    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        apflog( "Couldn't find any suitable targets!",level="error",echo=True)
        return None

    cadence_check = (ephem.julian_date(dt) - star_table[:, DS_LAST]) / star_table[:, DS_CAD]
    good_cadence = np.where(cadence_check > 0, True, False)
    good_cadence_available = available & good_cadence

    if len(good_cadence_available):
        try:
            pri = max(star_table[good_cadence_available, DS_APFPRI])
            sort_i = np.where(star_table[good_cadence_available, DS_APFPRI] == pri, True, False)
        except:
            pri = max(star_table[available, DS_APFPRI])
            sort_i = np.where(star_table[available, DS_APFPRI] == pri, True, False)
    elif len(available):
        pri = max(star_table[available, DS_APFPRI])
        sort_i = np.where(star_table[available, DS_APFPRI] == pri, True, False)
    else:
        apflog( "Couldn't find any suitable targets!",level="error",echo=True)
        return None

        
#    print sn[available][sort_i]
#    print star_table[available, DS_APFPRI][sort_i]
    starstr = "star table available: %s" % (sn[good_cadence_available][sort_i]) 
    apflog(starstr,echo=True)

    starstr = "star table available priorities: %s" % (star_table[good_cadence_available, DS_APFPRI][sort_i]) 
    apflog(starstr,echo=True)
     
    if bstar:
        sort_j = cur_elevations[good_cadence_available][sort_i].argsort()[::-1]
    else:
        sort_j = cadence_check[good_cadence_available][sort_i].argsort()[::-1]
        cstr= "cadence check: %s" %( cadence_check[good_cadence_available][sort_i][sort_j][0])
        apflog(cstr,echo=True)
    
    t_n = sn[good_cadence_available][sort_i][sort_j][0]

    elstr= "star elevations %s" % (cur_elevations[good_cadence_available][sort_i][sort_j])
    apflog(elstr,echo=True)

    t_n = sn[good_cadence_available][sort_i][sort_j][0]

    apflog("selected target %s" %( t_n) )

    idx, = np.where(sn == t_n)
    idx = idx[0]

    stars[idx].compute(apf_obs)
    
    res = dict()
    if pri > 9.9:
        star_table[idx, DS_EXPT] = MAX_EXPTIME

        
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
    res['I2CNTS'] = i2cnts[idx]
    res['NAME']   = sn[idx]
    res['SCORE']  = star_table[idx,DS_NSHOTS]
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['SCRIPTOBS'] = makeScriptobsLine(sn[idx], star_table[idx,:], do_flag[idx], dt)
    return res




if __name__ == '__main__':
    # For some test input what would the best target be?
    otfn = "observed_targets"
    ot = open(otfn,"w")
    starttime = time.time()
    result = getNext(starttime, 13.99, 1.8, bstar=True, verbose=True)
    ot.write("%s\n" % (result["SCRIPTOBS"]))
    ot.close()
    starttime += 400
    for i in range(5):
        
        result = getNext(starttime, 13.99, 1.8, bstar=False, verbose=True)
        #result = smartList("tst_targets", time.time(), 13.5, 2.4)

        if result is None:
            print "Get None target"
        else:
            for k in result:
                print k, result[k]
        ot = open(otfn,"a")
        ot.write("%s\n" % (result["SCRIPTOBS"]))
        ot.close()
        starttime += result["EXP_TIME"]
                
#    print "testing googledex updater"
#    update_googledex_lastobs('observed_targets')

    print "Done"
