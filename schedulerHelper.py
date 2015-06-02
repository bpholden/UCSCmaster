#!/usr/bin/env /opt/kroot/bin/kpython
import ephem
import subprocess
import os
from datetime import datetime, timedelta

import ktl

from apflog import *
# Helper function to generate scheduler.inp
# which is the input file for the dynamic scheduler


cwd = r'/u/rjhanson/master/'

# Outline of scheduler.inp
# Number of APF Stars
# UT twilight hr min sec
# UT sunrise hr min sec
# Telescope AZ El
# UT year month day
# UT Hour min sec
# Moon location (ra hr min sec sign dec deg min sec )
# Moon illumination percentage
# Bool - is the moon up?
# Shutter position - fscurpos rscurpos
# Windspeed wind direction

def schedInputs(AZ, EL, FSPOS, RSPOS, WS, WD, StarList):
    # We need the time at the start of an observing day
    # This is noon pacific.
    now = datetime.utcnow()
    if now.hour > 18:
        start = now - timedelta(hours=now.hour-18)
    else:
        start = now - timedelta(hours=now.hour+5)

    # Set up an observer at the APF Dome
    apf = ephem.Observer()
    apf.long = '-121:38:17.7'
    apf.lat  = '37:20:33.1'
    apf.elevation = 1274
    apf.horizon = '-9.0'
    apf.date = start
 
    sun = ephem.Sun()
    sun.compute(apf)
    moon = ephem.Moon()
    moon.compute(apf)

    stars = 0
    with open(StarList,'r') as f:
        for line in f: stars += 1

    with open(cwd + 'scheduler.inp','w') as out:
        # Number of stars
        out.write(repr(stars) + '\n')
        # Twilight
        twi = apf.next_setting(sun).datetime()
        out.write("%d %d %d\n" % (twi.hour, twi.minute, twi.second))
        # Dawn
        dawn = apf.next_rising(sun).datetime()
        out.write("%d %d %d\n" % (dawn.hour, dawn.minute, dawn.second))
        # Tele AZ and EL
        out.write("%f %f\n" % (AZ, EL))
        # UT year day month
        out.write("%d %d %d\n" % (now.year, now.month, now.day))
        # UT Hour minute second
        out.write("%d %d %d\n" % (now.hour, now.minute, now.second))
        # Moon current location (RA and Dec)
        m_ra = str(moon.ra).split(':')
        out.write(' '.join(m_ra))
        m_dec = str(moon.dec).split(':')
        if int(m_dec[0]) < 0:
            m_dec[0] = m_dec[0].lstrip('-')
            out.write(' - ' + ' '.join(m_dec) + '\n')
        else:
            out.write(' + ' + ' '.join(m_dec) + '\n')

        # Moon illumination percentage
        out.write("%f\n" % moon.phase)
        # Is the moon up?
        moonRise = apf.next_rising(moon).datetime()
        moonSet  = apf.next_setting(moon).datetime()
        if moonRise > moonSet:
            if now < moonSet:
                out.write(".true.\n")
            else:
                out.write(".false.\n")
        elif moonRise < now and now < moonSet:
            out.write(".true.\n")
        else:
            out.write(".false.\n")

        # Front and Rear Shutter positions
        out.write("%f %f\n" % (FSPOS, RSPOS))
        # Wind Speed and Wind Direction
        # Need to figure out what the last number is!
        out.write("%f %f 6.0\n" % (WS, WD))

def getObs():
    # Grab the telescope params
    tel = ktl.Service('eostele')
    el = float(tel("AEL").read())
    az = float(tel("AAZ").read())
    dome = ktl.Service('eosdome')
    fspos = float(dome("FSCURPOS").read())
    rspos = float(dome("RSCURPOS").read())
    ckapf = ktl.Service("checkapf")
    ws = ckapf("AVGWSPEED").read(binary=True)
    wd = ckapf("AVGWDIR").read(binary=True)
        
        
    apf_stars = "/u/rjhanson/master/stars_APF"
    schedInputs(az, el, fspos, rspos, ws, wd, apf_stars)
    
    try:
        ret = subprocess.check_call([cwd + 'windshield_all_apf_stars.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    except subprocess.CalledProcessError:
        print "Error with windshield"
        #apflog("Error occured in windshield_all_apf_stars.",echo=True)
        return None

    try:
        ret = subprocess.Popen([cwd + 'apf_scheduler_dyn1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    except subprocess.CalledProcessError:
        print "Error from dyn scheduler"
        #apflog("Error occured in apf_scheduler_dyn1",echo=True)
        return None
    output, err = ret.communicate()
    print output, err

    return cwd + 'apf_sched.txt'

def cleanup():
    """Sets up for the next nights observing."""
    files = ['apf_sched.txt','hit_list.old']
    for f in files:
        try:
            os.remove(cwd + f)
        except OSError:
            apflog("Couldn't delete %s" % f,echo=True)

    try:
        os.rename(cwd + "hit_list", cwd + "hit_list.old")
    except OSError:
        apflog("Couldn't rename hit_list", echo=True)

    # Slgihtly jenky solution to force files created during the night to be in UCSCAPF
    files = ['observed_targets', 'hit_list.old', 'diagnostics_apf.76', 'lastObs.txt', \
             'scheduler.inp', 'slewwind.out']
    for f in files:
        try:
            os.chown(cwd + f, -1, 1018)
            # Set the file permissions for the user and group
            os.chmod(cwd + f, 0664)
        except OSError:
            print "Couldn't reset group on %s" % (cwd + f)

    logpush(cwd + "observed_targets")
    logpush(cwd + "diagnostics_apf.76")

    
    
if __name__ == '__main__':
    print "Generating apf_sched.txt"

    fn = getObs()

    print "Returned: %s" % fn
    
