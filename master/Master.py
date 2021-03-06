#! @KPYTHON@
from __future__ import print_function

# UCSC script for the master task.
# Monitors and operates the APF for an observing night

import argparse
import atexit
from datetime import datetime, timedelta
import os
import os.path
import shutil
import signal
from select import select
import re
import subprocess
import sys
import threading
import time
import pickle
import functools
import glob

import numpy as np
import astropy.io.ascii

try:
    import ktl
    import APF as APFLib
    import APFTask
except:
    pass

import APFControl 
from apflog import *
import UCOScheduler_V1 as ds
from x_gaussslit import *
import ParseUCOSched
from Observe import Observe

os.umask(0007)

success = False

parent = 'master'

# global
paused = False

def controlWatch(keyword,parent):
    if keyword['populated'] == False:
        return
    try:
        value = keyword['ascii']
        global paused
        if value == "Abort":
            APFTask.set(parent, suffix='STATUS', value='Exited/Failure')
            APF.log("Aborted by APFTask")
            os.kill(os.getpid(), signal.SIGINT)
        elif value == "Pause" and not paused:
            try:
                APFTask.set(parent, suffix='STATUS', value='PAUSED')
                paused = True
            except:
                APF.log("Failure to set STATUS in APFTask", level="error")
                os.kill(os.getpid(),signal.SIGINT)

        else:
            try:
                APFTask.set(parent, suffix='STATUS', value='Running')
                paused = False
            except:
                APF.log("Failure to set STATUS in APFTask", level=error)
                os.kill(os.getpid(), signal.SIGINT)

    except:
        return


def signalShutdown(signal,frame):
    shutdown()
    
def shutdown():
    if success == True:
        status = 'Exited/Success'
        
    else:
        status = 'Exited/Failure'

    try:
        APFTask.set(parent, 'STATUS', status)
        sys.exit()
    except:   
        print('Exited/Failure')
        os._exit(1)
    else:
        print(status)
        os._exit(0)

def getStartTime(hr,mn):
    ct = datetime.now()
    st = datetime(ct.year, ct.month, ct.day+1, hr, mn)
    return float(st.strftime("%s"))

def calcFocusStartTime():
    # computes time 3.25 hours before sunset
    udt = datetime.utcnow()
    dt = datetime.now()    
    time_to_sunset = ds.computeSunset(udt)
    if time_to_sunset > 36000:
        # the sun has already set
        start_time = -1.
    else:
        start_time = time_to_sunset - 3.25*3600.
    t_delta = timedelta(0,start_time)
    sun_delta = timedelta(0,time_to_sunset)
    start_dt = dt + t_delta
    start_str = "%s" % (start_dt)
    sun_dt = dt + sun_delta
    sun_str = "%s" % (sun_dt)    
    return start_time, start_str, sun_str

def args():
    p_c = ["Init", "Focus", "Cal-Pre", "Cal-Post", "Watching"]
    w_c = ["on", "off", "auto"]
    b_c = [1, 2, 4]
    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-n', '--name', type=str,default='apf', help='This values is used as the UCAM observer name, as well as the file prefix.')
    parser.add_argument('-o', '--obsnum', type=int, help='Sets the UCAM observation number to this integer value.')
    parser.add_argument('-b', '--binning', choices=b_c, default=1, type=int, help='Sets the UCAM binning, bins in both pixels, allowed to be 1, 2 or 4.')
    parser.add_argument('-p', '--phase', choices=p_c, help='Specify the starting phase of the watcher. Allows for skipping standard procedures.')
    parser.add_argument('-f', '--fixed', default=None, help='Specify a fixed target list to observe. File will be searched for relative to the current working directory.')
    parser.add_argument('-t', '--test', action='store_true', help="Start the watcher in test mode. No modification to telescope, instrument, or observer settings will be made.")
#    parser.add_argument('-r', '--restart', action='store_true', default=False, help="Restart the specified fixed star list from the beginning. This resets scriptobs_lines_done to 0.") # removed should possible make default True and this option be False
    parser.add_argument('-w', '--windshield', choices=w_c, default='auto', help="Turn windshielding on, off, or let the software decide based on the current average wind speed (Default is auto). Velocity > %.1f mph turns windshielding on." % (APFControl.WINDSHIELD_LIMIT))
    parser.add_argument('-c', '--calibrate', default='uco', type=str, help="Specify the calibrate script to use. Specify string to be used in calibrate 'arg' pre/post")
    parser.add_argument('-l', '--line', type=int, help="If a fixed starlist is given, starts the list at line N.")
    parser.add_argument('-s', '--start', default=None, type=str, help="When specified with a fixed starlist, this option starts that list at that time.")
    parser.add_argument('--raster', action='store_true', default=False, help="If a fixed starlist is given, use it for a raster scan.")

    parser.add_argument('--rank_table', default='2020B_ranks', help="Optional name for table of sheet ranks")
    parser.add_argument('--sheet',default="RECUR_A100",help="Optional name for a Google spreadsheet")
    parser.add_argument('--owner',default='public',help="Optional name for file owners")    
    parser.add_argument('--frac_table',default='2020B_frac',help="Table of fractions of the night")
    parser.add_argument('--time_left',default=None,help="File name for remaining time in programs")
    
    opt = parser.parse_args()

    if opt.start != None:

        mtch = re.search("(\d+)\:(\d+)",opt.start)
        if mtch:
            hr = int(mtch.group(1))
            mn = int(mtch.group(2))
            opt.start = getStartTime(hr,mn)
        else:
            print("Start time %s does not match required format hours:minutes where both the hours and the minutes are integers")
            shutdown()


    if opt.sheet != None:
        opt.sheet = opt.sheet.split(",")

        
    return opt


def findObsNum(last):

    last += 100 - (last % 100)
    if last % 1000 > 700:
        last += 1000 - (last % 1000)
    if last >= 20000:
        last = 10000

    return last

def setObsDefaults(opt):
    opt.owner = 'public'
    if opt.name == None:
        opt.name = 'apf'
    if opt.obsnum == None:
        apflog("Figuring out what the observation number should be.",echo=False)
        opt.obsnum = findObsNum(int(ktl.read('apftask','MASTER_LAST_OBS_UCSC',binary=True)))
    else:
        opt.obsnum = int(opt.obsnum)

    return opt


def downloadFiles(opt):

    try:
        star_table,stars = ParseUCOSched.parseUCOSched(sheetns=opt.sheet,outfn='googledex.dat',outdir=os.getcwd())
    except Exception as e:
        apflog("Error: Cannot download googledex?! %s" % (e),level="error")
        # goto backup
        if os.path.exists("googledex.dat.1"):
            shutil.copyfile("googledex.dat.1","googledex.dat")

    try:
        rank_table = ds.makeRankTable(sheet_table_name=opt.rank_table,outdir=os.getcwd())
    except Exception as e:
        apflog("Error: Cannot download rank_table?! %s" % (e),level="error")
        # goto backup
        if os.path.exists("rank_table.1"):
            shutil.copyfile("rank_table.1","rank_table")


    if opt.time_left is None:
        hour_constraints=None
    else:
        if os.path.exists(opt.time_left):
            try:
                hour_constraints = astropy.io.ascii.read(opt.time_left)
            except Exception as e:
                hour_constraints = None
                apflog("Error: Cannot read file of time left %s : %s" % (opt.time_left,e))
            
    try:
        hour_table = ds.makeHourTable(opt.frac_table,datetime.now(),outdir=os.getcwd(),hour_constraints=hour_constraints)
    except Exception as e:
        apflog("Error: Cannot download frac_table?! %s" % (e),level="error")

if __name__ == '__main__':

    apflog("Starting master...")

    # Parse the command line arguments
    opt = args()

    # Register the atexit function after parsing the command line arguments
    # This prevents printing the help followed by exited/failure message
    atexit.register(shutdown)
    signal.signal(signal.SIGINT,  signalShutdown)
    signal.signal(signal.SIGTERM, signalShutdown)

    # Log the Command line arguments
    apflog("Command Line Args:")
    od = vars(opt)
    for o in od:
        if od[o] is not None:
            apflog(o + ': ' +str(od[o]))


    if opt.test:
        debug = True
        parent = 'example'
        apflog("master starting in test mode.")
    else:
        debug = False
        parent = 'master'


    apftask = ktl.Service("apftask")        
    # Establish this as the only running master script ( Or example task in test mode )
    try:
        apflog("Attempting to establish apftask as %s" % parent)
        APFTask.establish(parent, os.getpid())
    except Exception as e:
        apflog("Cannot establish as name %s: %s." % (parent,e), echo=True)
        shutdown()
    else:
        # Set up monitoring of the current master phase
        phase = apftask("%s_PHASE" % parent)
        phase.monitor()

    # Set preliminary signal and tripwire conditions
    apflog("Setting APFTask signal and tripwire.")
    APFTask.set(parent, "SIGNAL", "TERM")
    APFTask.set(parent, "TRIPWIRE", "TASK_ABORT")

    control = apftask[parent + '_CONTROL']
    cw = functools.partial(controlWatch, parent=parent)
    control.callback(cw)
    control.monitor(wait=False)

    apflog("Master initiallizing APF monitors.", echo=True)

    # Aquire an instance of the APF class, which holds wrapper functions for controlling the telescope
    apf = APFControl.APF(task=parent, test=debug)
    APFTask.waitFor(parent, True, timeout=5)
    apf.initGuideCam()
    
    # All the phase options that this script uses. This allows us to check if we exited out of the script early.
    possible_phases = ["Init", "Focus", "Cal-Pre", "Watching", "Cal-Post", "Focus-Post"]
    phase_index = 0
    # If a command line phase was specified, use that.
    if opt.phase != None:
        APFTask.phase(parent, opt.phase)
        phase_index = possible_phases.index(str(opt.phase).strip())
        phase.poll()
    elif phase.read() != "":
        # start where the master script left off
        phase_index = possible_phases.index(str(phase.read()).strip())
        phase.poll()
        
    # If the phase isn't a valid option, (say the watchdog was run last)
    # then assume we are starting a fresh night and start from setting the observer information.
    if str(phase).strip() not in possible_phases:
        apflog("Starting phase is not valid. Phase being set to Init", echo=True)
        APFTask.phase(parent, "Init")
    apflog("Phase at start is: %s" % phase, echo=True)
    APFLib.write(apf.robot["SCRIPTOBS_MESSAGE"], '')

    # Start the actual operations
    # Goes through 5 steps:
    # 1) Set the observer information
    # 2) Run focuscube
    # 3) Run calibrate ucsc pre
    # 4) Start the main watcher
    # 5) Run calibrate ucsc post
    # Specifying a phase jumps straight to that point, and continues from there.


    # Make sure that the command line arguments are respected.
    # Regardless of phase, if a name, obsnum, or reset was commanded, make sure we perform these operations.
    
    if "Init" == str(phase).strip():
        apflog("Setting the task step to 0")
        APFTask.step(parent, 0)
        apflog("Setting Observer Information", echo=True)
        opt = setObsDefaults(opt)
        apflog("Using %s for name and %s for obs number." % (opt.name, repr(opt.obsnum)), echo=True)
        apf.setObserverInfo(num=opt.obsnum, name=opt.name, owner=opt.owner)


        APFLib.write(apf.robot["SCRIPTOBS_MESSAGE"], "Setting defaults for observing.")
        if opt.fixed:
            if not os.path.exists(opt.fixed):
                errmsg = "starlist %s does not exist" % (opt.fixed)
                apflog(errmsg, level="error", echo=True)
                shutdown()
            if not debug:
                APFTask.set(parent, "STARLIST", opt.fixed)
        else:
            if not debug:
                APFTask.set(parent, "STARLIST", "")
        
        if os.path.exists(os.path.join(os.getcwd(),"googledex.dat")):
            logpush(os.path.join(os.getcwd(),"googledex.dat"))
        
        if os.path.exists(os.path.join(os.getcwd(),"robot.log")):
            logpush(os.path.join(os.getcwd(),"robot.log"))
            
        if os.path.exists(os.path.join(os.getcwd(),"googledex.dat")):
            logpush(os.path.join(os.getcwd(),"googledex.dat"))

        apflog("Setting SCRIPTOBS_LINES_DONE to 0")
        APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
        APFLib.write(apf.robot["MASTER_OBSBSTAR"], True,binary=True)        
        apflog("Initialization finished")
        
        stime, s_str, sun_str = calcFocusStartTime()
        waitstr = "Will now wait %.1f seconds before starting focusinstr" % (stime)
        apflog(waitstr, echo=True)
        APFTask.set(parent, suffix="MESSAGE", value=waitstr, wait=False)
        startstr = "Estimating focusinstr will start at %s for a sunset at %s" % (s_str,sun_str)
        apflog(startstr, echo=True)
        APFTask.set(parent, suffix="MESSAGE", value=startstr, wait=False)
        
        APFTask.wait(parent, True, timeout=stime)
        phase_index += 1
        APFTask.phase(parent, possible_phases[phase_index])
        apflog("Phase is now %s" % phase,echo=True)
        
    # 2) Run autofocus cube
    if "Focus" == str(phase).strip():
        
        apflog("Starting focusinstr script.", level='Info', echo=True)

        result = apf.ucam_status()
        if result is False:
            apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
            shutdown()
            
        result = apf.focusinstr()
        
        apflog("Focus has finished. Setting phase to Cal-Pre")
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())
        phase_index += 1
        APFTask.phase(parent, possible_phases[phase_index])
        apflog("Phase now %s" % phase)

    # 3) Run pre calibrations
    if 'Cal-Pre' == str(phase).strip():
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

        apflog("Starting calibrate pre script.", level='Info', echo=True)
        apf.instrPermit()

        result = apf.ucam_status()
        if result is False:
            apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
            shutdown()
        
        result = apf.calibrate(script = opt.calibrate, time = 'pre')
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

        if result == False:
            apflog("Calibrate Pre has failed. Trying again",level='warn',echo=True)
            apf.instrPermit()
            result = apf.calibrate(script = opt.calibrate, time = 'pre')
            if not result:
                apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                apf.turnOffLamps()
                sys.exit(2)

        phase_index += 1

        APFTask.phase(parent, possible_phases[phase_index])
        apflog("Calibrate Pre has finished. Setting phase to %s." % phase)


    # 4) Start the main watcher thread
    observe = Observe(apf,opt)
    if 'Watching' == str(phase).strip():
        apf.instrPermit()
        apflog("Starting the main watcher." ,echo=True)

        downloadFiles(opt)
        
        bstr = "%d,%d" % (opt.binning,opt.binning)
        apf.ucam['BINNING'].write(bstr)

        if opt.binning > 1:
            apfmon = ktl.Service('apfmon')
            d = time.time() + 22.*3600
            try:
                apfmon['BINNINGDIS'].write(d,binary=True,timeout=2)
            except:
                apflog("cannot write apfmon7 keyword" % (e),level='Alert',echo=True)
            
        if opt.fixed != None:
            lastList = apf.robot["MASTER_STARLIST"].read()
            # This resets lines done if this is a new target list
            if opt.fixed != lastList:
                APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
                APFLib.write(apf.robot["MASTER_STARLIST"], opt.fixed)
            if opt.line != None:
                APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], int(opt.line))
                apflog("Will be starting star list %s at line %d" % (opt.fixed, int(opt.line)),echo=True)
            else:
                apflog("Will be starting star list %s at line 0" % opt.fixed,echo=True)
        else:
            apflog("Starting dynamic scheduler", echo=True)
        observe.task = parent

        observe.start()
    else:
        observe.signal = False

    while observe.signal:
        # Observe is running, check for keyboard interupt
        try:
            currTime = datetime.now()
            # Check if it is after ~8:00AM.
            # If it is, something must be hung up, so lets force
            #  a closeup and run post cals. 
            #  Log that we force a close so we can look into why this happened.
            if currTime.hour == 8:
                # Its 8 AM. Lets closeup
                observe.stop()
                observe.signal=False
                apflog("Observe was still running at 8AM. It was stopped and post calibrations will be attempted.", level='warn')
                break
            if observe.isAlive() is False:
                # this means that the thread died, usually because of an uncaught exception
                # observe.signal will still be True
                # overwrite observe with a new instance and restart the thread
                apflog("Observe thread has died. Will attempt to restart.", level='error', echo=True)
                apf.killrobot()
                observe = Observe(apf,opt)
                observe.task = parent

                observe.start()
                
            if debug:
                print('Observe is running.')
                print(str(apf))
            APFTask.waitFor(parent, True, timeout=30)
        except KeyboardInterrupt:
            apflog("master has been killed by user.", echo=True)
            observe.stop()
            shutdown()
        except:
            apflog("master killed by unknown.", echo=True)
            observe.stop()
            shutdown()

    # Check if the observe left us an exit message.
    # If so, something strange likely happened so log it.
    try:
        msg = observe.exitMessage
    except AttributeError:
        pass
    else:
        apflog(msg, level='Info', echo=True)

    # We have finished taking data, and presumably it is the morning.
    apf.setTeqMode('Morning')
    apf.close(force=True)

    # Keep a copy of observed_targets around for a bit just in case
    if os.path.exists(os.path.join(os.getcwd(),"observed_targets")):
        try:
            apflog("Updating the online googledex with the observed times", level='Info', echo=True)
            ParseUCOSched.updateSheetLastobs(os.path.join(os.getcwd(),"observed_targets"),sheetns=opt.sheet,outfn="googledex.dat")


        except Exception as e:
            apflog("Error: Updating the online googledex has failed: %s" % (e), level="error")
        logpush(os.path.join(os.getcwd(),"observed_targets"))

    if os.path.exists(os.path.join(os.getcwd(),"googledex.dat")):
        logpush(os.path.join(os.getcwd(),"googledex.dat"))
        
    if os.path.exists(os.path.join(os.getcwd(),"rank_table")):
        logpush(os.path.join(os.getcwd(),"rank_table"))
        
    if os.path.exists(os.path.join(os.getcwd(),"hour_table")):
        logpush(os.path.join(os.getcwd(),"hour_table"))
        
    if os.path.exists(os.path.join(os.getcwd(),"frac_table")):
        os.remove(os.path.join(os.getcwd(),"frac_table"))
        
    if os.path.exists(os.path.join(os.getcwd(),"robot.log")):
        logpush(os.path.join(os.getcwd(),"robot.log"))
        

    apfmon = ktl.Service('apfmon')
    try:
        if apfmon['BINNINGDIS'].read(binary=True,timeout=2) > 0:
            apfmon['BINNINGDIS'].write(0,binary=True)
    except:
        apflog("cannot read apfmon7 keyword" % (e),level='Alert',echo=True)

    try:
        apf.ok2open.monitor(start=False)
    except Exception, e:
        apflog("Note: Cannot stop monitoring ok2open. %s" % (e), level="warn", echo=True)
    # 5) Take morning calibrations
    phase_index += 1
    APFTask.phase(parent, possible_phases[phase_index])
    apf.instrPermit()
    result = apf.calibrate(script=opt.calibrate, time='post')
    if not result:
        apflog("Calibrate Post has failed.", level='warn',echo=True)
        result = apf.calibrate(script=opt.calibrate, time='post')
        if not result:
            apflog("Calibrate Post has failed twice.", level='error',echo=True)
            APFTask.set(parent, suffix="MESSAGE",value="Calibrate Post failed twice",wait=False)

    if not debug:
        APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

    bstr = "%d,%d" % (1,1)
    apf.ucam['BINNING'].write(bstr) 
            
    # Focus the instrument once more
    phase_index += 1
    APFTask.phase(parent, possible_phases[phase_index])
    apflog("Running Focus-Post", echo=True)
    result = apf.focusinstr()
    if not result:
        apflog("Focus-post has failed", level='error', echo=True)
    else:
        apflog("Focus-post has finished successfully.", echo=True)

    # We have done everything we needed to, so leave
    # the telescope in day mode to allow it to start thermalizing to the next night.
    apf.setTeqMode('Day')

    # Update the last observation number to account for the morning calibration shots.

    if not debug:
        APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())
        APFTask.set(parent, suffix="MESSAGE",value="Updating last observation number to %s" % (apf.ucam["OBSNUM"].read()),wait=False)


    rv = apf.turnOffInst()
    APFTask.set(parent, suffix="MESSAGE",value="Turning off the motors",wait=False)

    # All Done!
    APFTask.phase(parent, "Finished")

    success = True
    shutdown()


