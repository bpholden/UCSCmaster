#!/opt/kroot/bin/kpython

import pyfits
import numpy as np

import thread
import time
import os
from datetime import datetime

import ktl
import APFTask
from apflog import apflog

from tharfwhm import tharfwhm

import matplotlib
matplotlib.use('Agg')
import pylab as pl


apfmot = ktl.Service('apfmot')
dewfoc = apfmot['DEWARFOCRAW']
ucam = ktl.Service('apfucam')

dewfoc = apfmot['DEWARFOCRAW']
dewfoc.monitor()
startfoc = dewfoc.binary

numsteps = 21
stepsize = 100
tharexp = 10.0
normfoc = 8200

ucscoff = 305.

focuspos = {'HATCHPOS':     'Closed',
            'DECKERNAM':    'Pinhole',
            #'DECKERNAM':    'N (0.50:8.0)',
            'IODINENAM':    'Out',
            'CALMIRRORNAM': 'In',
            'CALSOURCENAM': 'Thorium2',
            'HALOGEN1':     'Off',
            'HALOGEN2':     'Off',
            'THORIUM1':     'Off',
            'THORIUM2':     'On',}

#focuspos = {'HATCHPOS':     'Closed',
#            'DECKERNAM':    'Pinhole',
#            #'DECKERNAM':    'N (0.50:8.0)',
#            'IODINENAM':    'In',
#            'CALMIRRORNAM': 'In',
#            'CALSOURCENAM': 'Halogen Sphere',
#            'HALOGEN1':     'Off',
#            'HALOGEN2':     'On',
#            'THORIUM1':     'Off',
#            'THORIUM2':     'Off',}

obspos = {'HATCHPOS':       'Closed',
            'DECKERNAM':    'W (1.00:3.0)',
            'IODINENAM':    'In',
            'CALMIRRORNAM': 'Out',
            'CALSOURCENAM': 'Halogen Sphere',
            'HALOGEN1':     'Off',
            'HALOGEN2':     'Off',
            'THORIUM1':     'Off',
            'THORIUM2':     'Off',}

def waitwrite(keyword, value, timeout=180):
    apflog("\nSetting and waiting for %s = %s..." % (keyword.name, str(value)),echo=True)
    keyword.write(value)
    result = keyword.waitfor('== "%s"' % str(value),timeout=timeout)
    if result == False:
        raise ktl.ktlError, "Failed to set %s = %s" % (keyword.name, str(value))
    return result

def stagesetup(positions):
    for key,value in positions.items():
        thread.start_new_thread(waitwrite, (apfmot[key], value))

    for key,value in positions.items():
        result = apfmot[key].waitfor('== "%s"' % value, timeout=60)

    apflog("Spectrograph stages are now in position.",echo=True)

def exposure(objnam,exptime, readwait=True):
    result = ucam['EVENT_STR'].waitfor('== "ControllerReady"', timeout=120)
    if result == False:
        apflog("ERROR: Camera failed to become ready in 2 min.", level="error")
        raise ktl.ktlError, "Camera failed to become ready"

    stat = waitwrite(ucam['OBJECT'], objnam)
    stat = waitwrite(ucam['EXPOSURE'], exptime)

    apflog("Starting %3.1f s exposure on target %s..." % (exptime, objnam))
    ucam['START'].write('True')
    ucam['EVENT_STR'].waitfor('== "ExposureBegin"', timeout=10)
    if readwait: ucam['EVENT_STR'].waitfor('== "ControllerReady"')
    else: ucam['EVENT_STR'].waitfor('== "ReadoutBegin"')

def focusloop():
    apflog("Current focus value = %s" % dewfoc)

    focmul = np.array(range(numsteps)) - np.median(range(numsteps))
    focsteps = stepsize*focmul + normfoc
    apflog("Focus trials = %s" % repr(focsteps))
    
    imfiles = []
    for foc in focsteps:
        apflog("Setting dewar focus to %d" % foc)
        dewfoc.write(foc)
        time.sleep(2)
        while abs(dewfoc.binary - foc) > 10:
            time.sleep(2)
        apflog("Actual dewar focus value %d" % dewfoc.binary)
        fname = os.path.join(ucam('OUTDIR').read(), ucam('OUTFILE').read() + ucam('OBSNUM').read() + '.fits')
        imfiles.append(fname)
        exposure('Focus', tharexp, readwait=False)

    apflog("Returning focus to starting position...")
    dewfoc.write(startfoc, wait=False)

    return imfiles

def findfoc(images, doplot=False):


    focus = []
    fwhm = []
    flux = []
    for im in images:
        head = pyfits.getheader(im)
        foc = head['DFOCRAW']
        focus.append(foc)
        w,f = tharfwhm(im, silent=True)
        fwhm.append(w)
        flux.append(f)

    p = np.polyfit(focus, fwhm, 4)
    xmod = np.linspace(min(focus), max(focus), 200)
    yfit = np.polyval(p, xmod)

    bestfoc = xmod[np.argmin(yfit)]
    bestflux = flux[np.argmin(fwhm)]

    bestfoc += ucscoff
    
    offset = startfoc - bestfoc
    
    pl.plot(focus,fwhm, 'ko', markersize=10)
    pl.plot(xmod, yfit, 'b-')
    pl.ylabel('FWHM [px]')
    pl.xlabel('DEWARFOCRAW')
    pl.legend(['Best focus (w/ offset): %4.3f px\n@ %d' % (np.min(yfit),int(round(bestfoc)))], numpoints=1, loc='best')
    pl.savefig('dewar_focus.png')
    if doplot: pl.show()

    apflog("Found best focus at DEWARFOCRAW=%d" % int(round(bestfoc)))
    apflog("FHWM = %4.3f" % np.min(yfit))
    apflog("Max peak counts = %4.0f" % bestflux)
    apflog("Offset from previous value = %+4.0f" % offset)

    import mailer
    msgbody = """Focused the Levy spectograph at %s Pacific local time.
Found best focus (including offset) at DEWARFOCRAW=%d
Change from previous focus = %4.0f
FHWM = %4.3f pixels
Max peak counts = %4.0f ADU\n
""" % (datetime.strftime(datetime.now(), format="%Y-%m-%d %H:%M"), int(round(bestfoc)), int(round(offset)), np.min(yfit), bestflux)
    mailer.sendmail(msgbody, 'Levy Focus', toaddr=mailer.contacts.ucbgroup, attachments=['dewar_focus.png'])

    return int(round(bestfoc))

if __name__ == '__main__':
    
    APFTask.establish('focusinstr',os.getpid())
    try:
        waitwrite(ucam['AUTOSHUT'], True)
        stagesetup(focuspos)
        focimages = focusloop()
        bestfoc = findfoc(focimages)
        if abs(bestfoc - normfoc) < 2500:
            apflog("Setting dewar focus to %d" % bestfoc)    # Disabled updating to best focus, use UCSC value instead
            apfmot['DEWARFOCRAW'].write(bestfoc)
            APFTask.set('focusinstr','LASTFOCUS',bestfoc)
            #pass
        else:
            apflog("ERROR: Best focus value (%d) is very far from nominal!" % bestfoc, level="error")
            APFTask.set('focusinstr','STATUS','Exited/Failure')
        
        stagesetup(obspos)

        APFTask.set('focusinstr','STATUS','Exited/Success')
    except Exception, e:
        apflog("Spectrograph focus failed %s" % (e), level="error")
        APFTask.set('focusinstr','MSG',e)
        APFTask.set('focusinstr','STATUS','Exited/Failure')

    
