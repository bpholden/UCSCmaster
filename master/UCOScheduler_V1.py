# UCOScheduler_V1.py

from datetime import datetime, timedelta
import ephem
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
import Visible

# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20 # this elevation is the physical minimum, below this the ADC does not work
TARGET_ELEVATION_HIGH_MIN = 45 # this elevation is the preferred one for stars that will be high in the sky
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX =  1* 60 * 60 # 1 hour
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

# Maximum single exposure time in seconds
MAX_EXPTIME = 1200.
MIN_EXPTIME = 600.
MIN_TOTOBS = 300.
MAX_I2 = 40000
MIN_I2 = 500
MAX_EXPMETER = 2e9

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

def computeMaxTimes(sn,exp_times):
    maxtimes = np.zeros_like(exp_times)
    maxtimes += TARGET_EXPOSURE_TIME_MAX
    
    return maxtimes

def float_keyval(instr):
    try:
        key,val = instr.split('=')
    except:
        return None
    try:
        retval = float(val)
    except:
        return None
    return retval

def parseStarlist(starlist):
    """ Parse a scriptobs-compatible starlist for the scheduler.

    names, star_table, lines, stars = parseStarlist(starlist)
    starlist - a filename

    names - a list of stars in the starlist
    star_table - a numpy array
    lines - a list of strings that can be used for scriptobs input
    stars - a list of pyEphem objects 
    """
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
                star.name = ls[0]
                star._ra = ephem.hours(":".join([ls[1], ls[2], ls[3]]))
                star._dec = ephem.degrees(":".join([ls[4], ls[5], ls[6]]))
                stars.append(star)
            
    return names, np.array(star_table), lines, stars

def findColumns(col_names,req_cols,opt_cols=[]):
    """ findColumns finds the indices for the column names in the list of required columns
    indices = findColumns(col_names, req_cols)
    
    indices - a list of indices, each index maps to where in col_names the column is found and in the order of req_cols
    col_names - list of column names to be searched
    req_cols - list of names that should be in the first list
    """
    idx = []
    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) , level="Alert",echo=True)

    for r in opt_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
            
    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in didx.keys():
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Error",echo=True)

    return didx

def checkflag(key,didx,line,regexp,default):
    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default
    

def parseStarTable(fn="googledex.dat",wdir=None,config={'I2': 'Y', 'decker': 'W' }):
    """ parseStarTable parses the google sheet and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseStarTable(fn="googledex.dat",wdir=None,config={'I2': 'Y', 'decker': 'W' })
    names - a list of stars in the starlist
    star_table - a numpy array
    flags - a dictionary of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    if not outdir :
        outdir = os.getcwd()
    try:
        f = open(os.path.join(wdir,fn),'r')
    except IOError, e:
        print "Cannot open %s: %s" % (os.path.join(wdir,fn),e)
        return [], [], [], []
    full_codex = pickle.load(f)
    f.close()
        
    col_names = full_codex[0]
    codex = full_codex[1:]

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2"
                ]
    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    flags = { "do" : [], "decker" : [], "I2" : [] }
    stars = []
    # Build the star table to return to 
    for ls in codex:
        if ls[0] == '':
            continue
        try:
            apfpri = float(ls[didx["APFpri"]])
        except:
            apfpri = 0.0
        if apfpri < 0.5: continue
        row = []
        # Get the star name
        names.append(ls[didx["Star Name"]])
        # Get the RA
        raval = getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval:
            row.append(raval)
        else:
            row.append(-1.)
        # Get the DEC
        decval = getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
        if decval:
            row.append(decval)
        else:
            row.append(-3.14)

        for coln in ("pmRA", "pmDEC", "Vmag"):
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln == "Vmag":
                    row.append(15.0)
                else:
                    row.append(0.0)
        # For now use the old 1e9 count value
        row.append(1200.0)
        row.append(1.e9)
        for coln in ["APFpri", "APFcad","APFnshots"] :
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln in ("APFpri","APFnshots"):
                    row.append(0)
                else:
                    row.append(1000.0)
                    
        
        for coln in ["lastobs", "B-V", "APF Desired Precision" ]:
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln in ("lastobs", "B-V"):
                    row.append(0.0)
                else:
                    row.append(1000.0)

        check = checkflag("Close Companion",didx,ls,"\A(n|N)","Y")
        if check == "N" or check == "n":
            flags['do'].append("")
        else:
            flags['do'].append(check)
            
        flags['decker'].append(checkflag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        flags['I2'].append(checkflag("I2",didx,ls,"\A(n|N)",config["I2"]))
                                         
            
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]]))
        star._dec = ephem.degrees(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]]))
        stars.append(star)

    return (names, np.array(star_table), flags, stars)
    

def update_local_googledex(intime,googledex_file="googledex.dat", observed_file="observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time.
        update_local_googledex(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """
    names, times = getObserved(observed_file)

    try:
        g = open(googledex_file, 'r')
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
        return names,times

    full_codex = pickle.load(g)
    g.close()

    codex_cols = full_codex[0]

    starNameIdx = codex_cols.index("Star Name")
    lastObsIdx = codex_cols.index("lastobs")
    
    for i in range(1, len(full_codex)):
        row = full_codex[i]
        if row[starNameIdx] in names:
            # We have observed this star, so lets update the last obs field
            obstime = times[names.index(row[starNameIdx])]
            if isinstance(obstime,float):
                t = datetime.utcfromtimestamp(obstime)
            else:
                hr, min = obstime
                if type(intime) != datetime:
                    ctime = datetime.now()
                    td = timedelta(0,3600.*7)
                    intime = ctime + td
                t = datetime(intime.year, intime.month, intime.day, hr, min)

            # This keeps the JD precision to one decimal point. There is no real reason for this other than
            # the googledex currently only stores the lastObs field to one decimal precision. Legacy styles FTW.
            jd = round(float(ephem.julian_date(t)), 2) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)),echo=True)
            row[lastObsIdx] = str(jd)
            full_codex[i] = row

    with open(googledex_file, 'w') as f:
        pickle.dump(full_codex, f)
    f.close()
    
    return names, times


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


def getRARad(hr, mn, sec):
    try:
        ra_hours = float(hr) + float(mn)/60. + float(sec)/3600.
        return ra_hours * 15 * np.pi/180.0
    except:
        return None

def getDECRad(deg, mn, sec, neg=False):
    try:
        deg = float(deg)
        mn = float(mn)
        sec = float(sec)
    except:
        return None
    if deg < 0:
        neg = True
        deg = abs(deg)       
    if  mn < 0:
        neg = True
        mn = abs(mn)
    if sec < 0:
        neg = True
        sec = abs(sec)
    x = deg + mn/60. + sec/3600.
    x = x * np.pi/180.
    if neg:
        return x*-1
    else:
        return x

def getCoordStr(floatval,isRA=False):

    neg = False
    nround = 2
    if isRA:
        floatval /= 15.
        nround = 3
    if floatval < 0:
        neg = True
    floatval = abs(floatval)
    deghrval = int(floatval)
    minval = (floatval % 1) * 60.0 
    secval = round( (minval % 1) *60.0, nround)

    if neg and deghrval != 0:
        ret = "-" + str(deghrval) + ' '
    else:
        ret = str(deghrval) + ' '
    if neg and deghrval == 0 and minval != 0:
        ret += "-" + str(int(minval)) + ' '
    else:
        ret += str(int(minval)) + ' '
    if neg and deghrval == 0 and minval == 0:
        ret += "-" + str(secval)
    else:
        ret += str(secval)
    return ret

        
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

def makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y",owner='Vogt'):
    """ given a name, a row in a star table and a do_flag, will generate a scriptobs line as a string
    line = makeScriptobsLine(name, row, do_flag, t, decker="W",I2="Y")
    name - name of star, first column in line
    row - star_table row for star that begins with name, cotains all of the data needed for the line except
    do_flag - a string for whether or not scriptob needs to do a pointing check before slewing to the target
    t - a datetime object, this is used to fill in the uth and utm fields,
    decker - one character field for the decker, defaults to "W"
    I2 - one character field for whether or not the Iodine cell is in, must be "Y" or "N"
    """
    focval = 0
    """Takes a line from the star table and generates the appropriate line to pass to scriptobs. """
    # Start with the target name
    ret = name + ' '
    # Add the RA as three elements, HR, MIN, SEC
    rastr = getCoordStr(np.degrees(row[DS_RA]),isRA=True)
    ret += rastr + ' '
    # Add the DEC as three elements, DEG, MIN, SEC
    decstr = getCoordStr(np.degrees(row[DS_DEC]))
    ret += decstr + ' '
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
    elif row[DS_EXPT] <= MIN_EXPTIME:
        ret += 'texp=%d ' % (int(MIN_EXPTIME))
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
        
    return ret

def getObserved(filename):
    """ getObserved parses a file to find the object names and times
    names, times = getObserved(filename)
    names - list of names, must be first column of file called filename
    times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line

    """
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
    """ calculate_ucsc_exposure_time uses the recipe from Burt et al. (2015) to compute the exposure time for a target.

    exp_time, exp_counts, i2counts = calculate_ucsc_exposure_time(vmag, precision, elevation, seeing, bmv, decker="W")
    vmag - numpy array of V magnitudes (Johnson filter, Vega mags)
    precision - required precision for the velocity in m/s
    elevation - elevation of the star above the horizon at the start of the exposure
    seeing - FWHM of the seeing in pixels on the guider
    bmv - (B - V) for the star (both Johnson filters, Vega zeropoint)
    decker - apeture 
    
    exp_time - a numpy array of times in seconds, are integer values
    exp_counts - values for the exposure meter, this can be a floating point value
    i2counts - the required number of median Iodine cell counts, this is calculated from the precision and color of the star, this, in effect, sets the exposure time.


    """
    vmag = np.array(vmag)
    precision = np.array(precision)
    bmv = np.array(bmv)

		
	# Now lets calculate the exposure times
	
	# Desired I2 counts for precision
    i2counts = getI2_K(precision)
    mstars = np.where(bmv > 1.2)
    if len(mstars) > 0:
        i2counts[mstars] = getI2_M(precision[mstars])

    # minimum I2 counts so exposures are not rejected by P. Butler's DRP
    mini2_idx = np.where(i2counts < MIN_I2)
    if len(mini2_idx) > 0:
        i2counts[mini2_idx] = MIN_I2
	
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


def smartList(starlist, time, seeing, slowdown,outdir = None):
    """ Determine the best target to observe from the provided scriptobs-compatible starlist.
        Here the best target is defined as an unobserved target (ie not in observed targets )
        that is visible above 30 degrees elevation. Higher elevation targets are prefered,
        but those that rise above 85 degrees will be regected to avoid slewing through the zenith. """
    # Convert the unix timestamp into a python datetime


    if type(time) == float:
        dt = datetime.utcfromtimestamp(int(time))
    elif type(time) == datetime:
        dt = time
    elif type(time) == ephem.Date:
        dt = time.datetime()
    else:
        #punt
        dt = datetime.utcfromtimestamp(int(time.time()))
        
    if not outdir:
        outdir = os.getcwd()
    observed, _ = getObserved(os.path.join(outdir,"observed_targets"))

    # Generate a pyephem observer for the APF
    apf_obs = ephem.Observer()
    apf_obs.lat  = '37:20:33.1'
    apf_obs.long = '-121:38:17.7'
    apf_obs.elevation = 1274
    # Minimum observation to observe things at
    apf_obs.horizon = str(TARGET_ELEVATION_MIN)
    apf_obs.date = dt
    # APF latitude in radians
    apf_lat = apf_obs.lat

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
    vis, star_elevations, fin_els = Visible.is_visible(apf_obs,stars,obs_length)
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

    exps = np.zeros_like(exp_counts)
    exp_counts *= 1.1 
    long_idx = exp_counts > MAX_EXPMETER
    toofew_idx = exp_counts/nexp > MAX_EXPMETER
    nexp[toofew_idx] = np.ceil((exp_counts[toofew_idx]/MAX_EXPMETER) + 1)
    exp_counts[long_idx] = MAX_EXPMETER
    exps[exps < nexp] = nexp[exps < nexp]
    return exp_counts, exps

def format_time(total, i2counts, nexp, hitthemall=False):
    total = np.array(total)
    times = np.zeros(len(total))
    exps  = np.zeros(len(total))

    middle_idx = (total > MIN_TOTOBS ) &(total < MAX_EXPTIME)
    times[middle_idx] = MAX_EXPTIME # pad out to make it more likely exposure meter threshold sets actual limit
    exps[middle_idx] = 1

    max_idx = total > MAX_EXPTIME
    if hitthemall:
        exps[max_idx] = 1
    else:
        exps[max_idx] = np.ceil(total[max_idx]/MAX_EXPTIME)
    times[max_idx] = MAX_EXPTIME

    short_idx = total < MIN_TOTOBS
    times[short_idx] = MIN_EXPTIME  # pad out to make it more likely exposure meter threshold sets actual limit
    exps[short_idx] = [ (np.ceil(MIN_TOTOBS/(t+40)) + 1) for t in total[short_idx] ] 

    exps[exps < nexp] = nexp[exps < nexp]

    return times, exps


def getNext(ctime, seeing, slowdown, bstar=False, verbose=False,sheetn="The Googledex",owner='Vogt',outfn="googledex.dat",outdir=None):
    """ Determine the best target for UCSC team to observe for the given input.
        Takes the time, seeing, and slowdown factor.
        Returns a dict with target RA, DEC, Total Exposure time, and scritobs line
    """

    if not outdir:
        outdir = os.getcwd()
        
    # Convert the unix timestamp into a python datetime
    if type(ctime) == float:
        dt = datetime.utcfromtimestamp(int(ctime))
    elif type(ctime) == datetime:
        dt = ctime
    elif type(ctime) == ephem.Date:
        dt = ctime.datetime()
    else:
        dt = datetime.utcfromtimestamp(int(time.time()))
        # punt

    confg = dict()
    confg['I2'] = 'Y'
    confg['decker']='W'
                
    if verbose:
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
            
    observed, obstimes = update_local_googledex(ptime,googledex_file=os.path.join(outdir,"googledex.dat"), observed_file=os.path.join(outdir,"observed_targets"))

    # List of targets already observed
#    observed, _ = getObserved(os.path.join(os.getcwd(),'observed_targets'))
    global last_objs_attempted
    try:
        lastline = ktl.read("apftask","SCRIPTOBS_LINE")
        if not bstar:             # otherwise from previous night
            lastobj = lastline.split()[0]

        if verbose:
            apflog( "getNext(): Last object attempted %s" % (lastobj),echo=True)
    except:
        lastobj = None

    if lastobj:
        if lastobj not in observed and lastobj not in last_objs_attempted:
            last_objs_attempted.append(lastobj)
            if verbose:
                apflog( "getNext(): Last objects attempted %s" % (last_objs_attempted),echo=True)
            
            if len(last_objs_attempted) > 5:
                apflog( "getNext(): 5 failed acquisition attempts",echo=True)
                last_objs_attempted = []
                return None

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
    sn, star_table, flags, stars = parseGoogledex(sheetn=sheetn,outfn=outfn,outdir=outdir)
    sn = np.array(sn)
    targNum = len(sn)
    if verbose:
        apflog("getNext(): Parsed the Googledex...",echo=True)

    # Note which of these are B-Stars for later.
    bstars = np.array([ True if 'HR' in n else False for n in sn ], dtype=bool)
    if verbose:
        apflog("getNext(): Finding B stars",echo=True)


    # Distance to stay away from the moon
    md = TARGET_MOON_DIST_MAX - TARGET_MOON_DIST_MIN
    minMoonDist = ((moon.phase / 100.) * md) + TARGET_MOON_DIST_MIN  


    moonDist = np.degrees(np.sqrt((moon.ra - star_table[:,DS_RA])**2 + (moon.dec - star_table[:,DS_DEC])**2))

    available = np.ones(targNum, dtype=bool)
    totexptimes = np.zeros(targNum, dtype=float)
    cur_elevations = np.zeros(targNum, dtype=float)
    scaled_elevations = np.zeros(targNum, dtype=float)    
    i2cnts = np.zeros(targNum, dtype=float)

    # Is the target behind the moon?
    if verbose:
        apflog("getNext(): Culling stars behind the moon",echo=True)
    moon_check = np.where(moonDist > minMoonDist, True, False)
    available = available & moon_check

    # We just need a B star, so restrict our math to those
    if bstar:
        if verbose:
            apflog("getNext(): Selecting B stars",echo=True)
        available = available & bstars
        
        f = available
        if verbose:
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
        if verbose:
            apflog("getNext(): Culling B stars",echo=True)
        available = np.logical_and(available, np.logical_not(bstars))
        
        # has the star been observed - commented out as redundant with cadence
        if len(last_objs_attempted)>0:
            for n in last_objs_attempted:
                attempted = (sn == n)
                available = available & np.logical_not(attempted) # Available and not observed

        # Calculate the exposure time for the target
        # Want to pass the entire list of targets to this function
        f = available
        if verbose:
            apflog("getNext(): Computing star elevations",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.is_visible_se(apf_obs, fstars, [0]*len(fstars))
#        vis,star_elevations,fin_star_elevations = Visible.is_visible( apf_obs,fstars,[0]*len(fstars))
        available[f] = available[f] & vis
        f = available
        fstars = [s for s,_ in zip(stars,f) if _ ]
        if verbose:
            apflog("getNext(): Computing exposure times",echo=True)
        exp_times, exp_counts, i2counts = calculate_ucsc_exposure_time( star_table[f,DS_VMAG], \
                                            star_table[f,DS_ERR], star_elevations[np.array(vis)], seeing, \
                                            star_table[f,DS_BV])
        
        exp_times = exp_times * slowdown
        maxtimes = computeMaxTimes(sn[f],exp_times)
        totexptimes[f] += exp_times
        i2cnts[f] += i2counts
        if verbose:
            apflog("getNext(): Formating exposure times",echo=True)
        star_table[f, DS_EXPT], exps = format_time(exp_times,i2counts,star_table[f, DS_NSHOTS])

        if verbose:
            apflog("getNext(): Formating exposure meter",echo=True)
        star_table[f, DS_COUNTS], star_table[f, DS_NSHOTS] = format_expmeter(exp_counts,exps)
        #        star_table[f, DS_COUNTS] = 1.1*exp_counts / star_table[f, DS_NSHOTS]

        # Is the exposure time too long?
        if verbose:
            apflog("getNext(): Removing really long exposures",echo=True)
        time_check = np.where( exp_times < maxtimes, True, False)
        
        available[f] = available[f] & time_check
        f = available

        # Is the star currently visible?
        if verbose:
            apflog("getNext(): Computing stars visibility",echo=True)
        fstars = [s for s,_ in zip(stars,f) if _ ]
        vis,star_elevations,fin_star_elevations, scaled_els = Visible.is_visible_se(apf_obs, fstars, exp_times)
        if vis != []:
            available[f] = available[f] & vis
        cur_elevations[np.where(f)] += star_elevations[np.where(vis)]
        scaled_elevations[np.where(f)] += scaled_els[np.where(vis)]        

    # Now just sort by priority, then cadence. Return top target
    if len(sn[available]) < 1:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

    cadence_check = (ephem.julian_date(dt) - star_table[:, DS_LAST]) / star_table[:, DS_CAD]
    good_cadence = np.where(cadence_check >  1.0, True, False)
    good_cadence_available = available & good_cadence

    if len(good_cadence_available):
        try:
            pri = max(star_table[good_cadence_available, DS_APFPRI])
            sort_i = (star_table[:, DS_APFPRI] == pri) & good_cadence_available
        except:
            pri = max(star_table[available, DS_APFPRI])
            sort_i = (star_table[:, DS_APFPRI] == pri) & available
    elif len(available):
        apflog( "getNext(): No new stars available, going back to the previously observed list.",level="warn",echo=True)
        pri = max(star_table[available, DS_APFPRI])
        sort_i = (star_table[:, DS_APFPRI] == pri) & available
    else:
        apflog( "getNext(): Couldn't find any suitable targets!",level="error",echo=True)
        return None

        
#    print sn[available][sort_i]
#    print star_table[available, DS_APFPRI][sort_i]
    starstr = "getNext(): star table available: %s" % (sn[sort_i]) 
    apflog(starstr,echo=True)

    starstr = "getNext(): star table available priorities: %s" % (star_table[:, DS_APFPRI][sort_i]) 
    apflog(starstr,echo=True)
     
    if bstar:
        sort_j = cur_elevations[sort_i].argsort()[::-1]
    else:
        sort_j = scaled_elevations[sort_i].argsort()[::-1]
#        sort_j = cadence_check[sort_i].argsort()[::-1]
        cstr= "getNext(): cadence check: %s" %( cadence_check[sort_i][sort_j][0])
        apflog(cstr,echo=True)
    
    t_n = sn[sort_i][sort_j][0]

    elstr= "getNext(): star elevations %s" % (cur_elevations[sort_i][sort_j])
    apflog(elstr,echo=True)

    t_n = sn[sort_i][sort_j][0]

    apflog("getNext(): selected target %s" %( t_n) )

    idx, = np.where(sn == t_n)
    idx = idx[0]

    stars[idx].compute(apf_obs)
    cstr= "getNext(): cadence check: %f (%f %f %f)" %( ((ephem.julian_date(dt) - star_table[idx, DS_LAST]) / star_table[idx, DS_CAD]), ephem.julian_date(dt), star_table[idx, DS_LAST], star_table[idx, DS_CAD])
    apflog(cstr,echo=True)
    
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
    res['I2CNTS'] = i2cnts[idx]
    res['NAME']   = sn[idx]
    res['SCORE']  = star_table[idx,DS_NSHOTS]
    res['PRI']    = star_table[idx, DS_APFPRI]
    res['DECKER'] = flags['decker'][idx]
    res['I2'] = confg['I2']    
    res['SCRIPTOBS'] = makeScriptobsLine(sn[idx], star_table[idx,:], flags['do'][idx], dt, decker=flags['decker'][idx], I2=flags['I2'][idx], owner=owner)
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