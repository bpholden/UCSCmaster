#!/usr/bin/env  /opt/kroot/bin/kpython
# Heimdallr.py
# UCSC script for the master task.
# Monitors and operates the APF for an observing night

import argparse
import atexit
from datetime import datetime, timedelta
import numpy as np
import os
import os.path
from select import select
import subprocess
import sys
import threading
import time
import pickle
import functools
import glob

import ktl
import APF as APFLib
import APFTask

import APFControl as ad
from apflog import *
import UCSCScheduler_V2 as ds
from x_gaussslit import *
import ExposureCalculations

os.umask(0007)

success = False

parent = 'master'

SUNEL_ENDLIM = -10.0
SUNEL_STARTLIM = -9.0
SUNEL_HOR = -3.2
FOCUSTIME = 3600.
AVERAGE_INSTRFOC = 8522
DMLIM = 1140
# global
paused = False

def sunel_startlim():
    return SUNEL_STARTLIM

def control_watch(keyword,parent):
    if keyword['populated'] == False:
        return
    try:
        value = keyword['ascii']
        global paused
        if value == "Abort":
            APFTask.set(parent,suffix='STATUS',value='Exited/Failure')
            APF.log("Aborted by APFTask")
            os.kill(os.getpid(),signal.SIGINT)
        elif value == "Pause":
            try:
                APFTask.set(parent,suffix='STATUS',value='PAUSED')
                paused = True
            except:
                APF.log("Failure to set STATUS in APFTask",level=error)
                os.kill(os.getpid(),signal.SIGINT)

        else:
            try:
                APFTask.set(parent,suffix='STATUS',value='Running')
                paused = False
            except:
                APF.log("Failure to set STATUS in APFTask",level=error)
                os.kill(os.getpid(),signal.SIGINT)

    except:
        return


def shutdown():
    if success == True:
        status = 'Exited/Success'
        
    else:
        status = 'Exited/Failure'

    try:
        APFTask.set(parent, 'STATUS', status)
    except:   
        print 'Exited/Failure'
        os._exit(1)
    else:
        print status
        os._exit(0)


def args():
    p_c = ["ObsInfo", "Focus", "Cal-Pre", "Cal-Post", "Watching"]
    w_c = ["on", "off", "auto"]
    b_c = [1,2,4]
    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-n', '--name', type=str, help='This values is used as the UCAM observer name, as well as the file prefix.')
    parser.add_argument('-o', '--obsnum', type=int, help='Sets the UCAM observation number to this integer value.')
    parser.add_argument('-b', '--binning', choices=b_c, default=1, type=int, help='Sets the UCAM binning, bins in both pixels, allowed to be 1, 2 or 4.')
    parser.add_argument('-p', '--phase', choices=p_c, help='Specify the starting phase of the watcher. Allows for skipping standard procedures.')
    parser.add_argument('-f', '--fixed', help='Specify a fixed target list to observe. File will be searched for relative to the current working directory.')
    parser.add_argument('-t', '--test', action='store_true', help="Start the watcher in test mode. No modification to telescope, instrument, or observer settings will be made.")
#    parser.add_argument('-r', '--restart', action='store_true', default=False, help="Restart the specified fixed star list from the beginning. This resets scriptobs_lines_done to 0.") # removed should possible make default True and this option be False
    parser.add_argument('-w', '--windshield', choices=w_c, default='auto', help="Turn windshielding on, off, or let the software decide based on the current average wind speed (Default is auto). Velocity > %.1f mph turns windshielding on." % (ad.WINDSHIELD_LIMIT))
    parser.add_argument('-c', '--calibrate', default='ucsc', type=str, help="Specify the calibrate script to use. Specify string to be used in calibrate 'arg' pre/post")
    parser.add_argument('-l', '--line', type=int, help="If a fixed starlist is given, starts the list at line N.")
    parser.add_argument('-s', '--smartObs', default=False, action='store_true', help="When specified with a fixed starlist, this option will pass the starlist to a selection algorithm to choose the optimal unobserved target, regardless of starlist ordering.")
    parser.add_argument('--sheet',default=None,help="Optional name for a Google spreadsheet")
    parser.add_argument('--owner',default='Vogt',help="Optional name for file owners")    
    
    opt = parser.parse_args()
    return opt


def findObsNum(apf):

    obsNum = int(apf.robot["MASTER_LAST_OBS_UCSC"].read().strip())

    last_times = int(float(apf.robot["MASTER_VAR_2"].read()))
    deltat = time.time() - last_times
    deltat /= (24*3600)
    ndays = int(deltat+0.5)

    obsNum += 100 - (obsNum % 100)
    obsNum += ndays*200

    if obsNum % 10000 > 9700:
        obsNum += 10000 - (obsNum % 10000)

    return obsNum

def getTotalLines(filename):
    tot = 0
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            elif line.strip()[0] == '#':
                continue
            else:
                tot += 1
    return tot

               

class Master(threading.Thread):
    def __init__(self, apf, user='ucsc',sheetn="The Googledex",owner='Vogt'):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.APF = apf
        self.user = user
        self.owner = owner
        self.name = 'Heimdallr'
        self.signal = True
        self.windshield = 'auto'
        self.scriptobs = None
        self.BV = None
        self.VMAG = None
        self.decker = "W"
        self.obsBstar = False
        self.fixedList = None
        self.sheetn = sheetn
        self.targetlogname = os.path.join(os.getcwd(),"targetlog.txt")
        self.targetlog = None

        self.nighttargetlogname = os.path.join(os.getcwd(),"nighttargetlog.txt")
        self.nighttargetlog = None

    def set_autofocval(self):
        APF = self.APF
        # check last telescope focus
        lastfoc = APF.robot['FOCUSTEL_LAST_SUCCESS'].read(binary=True)
        if time.time() - lastfoc > FOCUSTIME :
            APF.autofoc.write("robot_autofocus_enable")
        else:
            APF.autofoc.write("robot_autofocus_disable")


    def run(self):
        APF = self.APF

        # This is called when an observation finishes, and selects the next target
        def getTarget():
            apflog("getTarget(): Scriptobs phase is input, determining next target.",echo=True)
            APFLib.write(APF.ucam["RECORD"], "Yes") # safe / sorry

            try:
                self.obsBstar = bool(ktl.read("apftask","master_var_3"))
            except:
                self.obsBstar = True
                
            if self.scriptobs is None:
                apflog("Called getTarget, but there is not instance of scriptobs associated with Heimdallr. This is an error condition.", echo=True)
                return None
            
            # Calculate the slowdown factor.
            if self.BV is None:
                apflog("Warning!: Ended up in getTarget() with no B Magnitude value, slowdown can't be computed.", echo=True)
                self.BV = 0.6 # use a default average

            if self.VMAG is None:
                apflog("Warning!: Ended up in getTarget() with no V Magnitude value, slowdown can't be computed.", echo=True)
                slowdown = 5
            elif APF.avg_fwhm < 1.0:
                apflog("Warning!: AVG_FWHM = %4.2f. By Odin's beard that seems low." % APF.avg_fwhm, echo=True)
                slowdown = 5
            else:
                apflog("Calculating expected counts")
                apflog("self.VMAG [%4.2f] - self.BV [%4.2f] - APF.ael [%4.2f]" % (self.VMAG, self.BV, APF.ael))
                exp_cnts_sec = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, APF.ael,APF.avg_fwhm, decker=self.decker)
                try:
                    slowdown = exp_cnts_sec / APF.countrate
                    if slowdown < 0:
                        slowdown = 1
                        apflog("Countrate non-sensical %g" % (APF.countrate), echo=True, level='warn')
                        APF.counts.monitor(start=False)
                        APF.counts.monitor(start=True)
                        APF.counts.callback(ad.countmon)
                        # yes this happened.
                    if slowdown < ds.SLOWDOWN_MIN:
                        slowdown = ds.SLOWDOWN_MIN
                        apflog("slowdown too low, countrate= %g" % (APF.countrate), echo=True, level='debug')
                        # yes this happened.
                    if slowdown > ds.SLOWDOWN_MAX:
                        slowdown = ds.SLOWDOWN_MAX
                        apflog("slowdown too high, countrate= %g" % (APF.countrate), echo=True, level='debug')
                except ZeroDivisionError:
                    apflog("Current countrate was 0. Slowdown will be set to 1.", echo=True)
                    slowdown = 1
            apflog("getTarget(): slowdown factor = %4.2f" % slowdown, echo=True)
            APFLib.write(apf.robot["MASTER_VAR_1"], slowdown)
            apflog("getTarget(): countrate = %.2f" % APF.countrate)

            # Check for a valid seeing measurment. If there isn't one, use a default
            if APF.avg_fwhm == 0.:
                apflog("getTarget(): Warning AVG_FWHM is 0. A default value of 15 will be used in its place.",echo=True)
                seeing = 15
            else:
                seeing = float(APF.avg_fwhm)
                apflog("getTarget(): Current AVG_FWHM = %4.2f" % seeing)
            
            if self.fixedList is None:
                # Pull from the dynamic scheduler
                target = ds.getNext(time.time(), seeing, slowdown, bstar=self.obsBstar, verbose=True, sheetn=self.sheetn, owner=self.owner)
            else:
                # Get the best target from the star list
                if os.path.exists(self.fixedList):
                    target = ds.smartList(self.fixedList, time.time(), seeing, slowdown)
                else:
                    apflog("Error: starlist %s does not exist" % (self.fixedList), level="error",echo=True)
                    self.fixedList=None

            self.set_autofocval()
            if target is None:
                apflog("No acceptable target was found. Since there does not seem to be anything to observe, Heimdallr will now shut down.", echo=True)
                # Send scriptobs EOF to finish execution - wouldn't want to leave a zombie scriptobs running
                self.scriptobs.stdin.close()
                APF.close()
                apf.countrate = -1.0
                # sleep for a half hour to see if the clouds blow by
                APFTask.waitfor(self.task, True, timeout=60*45)
                return
            else:
                apflog("Observing target: %s" % target['NAME'], echo=True)
                self.scriptobs.stdin.write(target["SCRIPTOBS"] + '\n')
            # Set the Vmag and B-V mag of the latest target
            self.VMAG = target["VMAG"]
            self.BV   = target["BV"]
            self.decker = target["DECKER"]
            apflog("getTarget(): V=%.2f  B-V=%.2f Pri=%.2f " % (self.VMAG, self.BV, target["PRI"]))
            apflog("getTarget(): FWHM=%.2f  Slowdown=%.2f  Countrate=%.2f" % (APF.avg_fwhm, slowdown, APF.countrate))

            apflog("getTarget(): Target= %s" % target["NAME"])
            apflog("getTarget(): Counts=%.2f  EXPTime=%.2f  Nexp=%d" % (target["COUNTS"], target["EXP_TIME"], target["NEXP"]))
            APF.updateLastObs()

        def opening(sunel,sunset=False):
            when = "night"
            if sunset:
                when = "sunset"
                
            apflog("Running open at %s as sunel = %4.2f" % (when, float(sunel)))
            (apfopen,what) =APF.isOpen()
            if apfopen:
                APF.DMReset()
            else:
                APF.DMZero()

            result = APF.openat(sunset=sunset)
            apflog("openatsunset completed with result %s" % (result), echo=True)
            if result == False:
                apflog("openatsunset hasn't successfully opened. Current sunel = %4.2f" % ( float(sunel)), level='warn', echo=True)
                if ( float(sunel) < SUNEL_ENDLIM):
                    result = APF.openat(sunset=sunset)
                    if not result:
                        apflog("Error: openatsunset has failed twice.", level='error', echo=True)
                        APF.close()


            if datetime.now().strftime("%p") == 'PM':
                setting = True
            else:
                setting = False
            APF.DMReset()
            chk_done = "$eostele.SUNEL < %f" % (SUNEL_STARTLIM*np.pi/180.0)
            result = False
            while float(sunel.read()) > SUNEL_STARTLIM and setting:
                outstr = "Sun setting is %s and sun at elevation of %.3f" % (setting, float(sunel.read()))
                apflog(outstr,level='info', echo=True)
                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=60)
                APF.DMReset()
            return

        def startScriptobs():
            # Update the last obs file and hitlist if needed

            APF.updateLastObs()
            APF.updateWindshield(self.windshield)
            apflog("Starting an instance of scriptobs",echo=True)
            ripd, running = APF.findRobot()
            if running:
                apflog("Scriptobs is already running yet startScriptobs was called",level="warn",echo=True)
                return

            if self.fixedList is not None and self.smartObs == False:
                # We wish to observe a fixed target list, in it's original order
                if not os.path.exists(self.fixedList):
                    apflog("Error: starlist %s does not exist" % (self.fixedList), level="error")
                    self.fixedList=None
                    return
                tot = getTotalLines(self.fixedList)
                apflog("%d total starlist lines and %d lines done." % (tot, APF.ldone)) 
                if APF.ldone == tot and APF.user != "ucsc":
                    APF.close()
                    APF.updateLastObs()
                    self.exitMessage = "Fixed list is finished. Exiting the watcher."
                    self.stop()
                    # The fixed list has been completely observed so nothing left to do
                elif APF.ldone == tot and APF.user == "ucsc":
                    self.fixedList = None
                    self.smartObs = False
                    APFTask.set(self.task,suffix="STARLIST",value="")
                    ripd, running = APF.findRobot()
                    if running:
                        APF.killRobot(now=True)
                    apflog("Finished fixed list on line %d, will start dynamic scheduler" % int(APF.ldone), echo=True)
                    expression="($apftask.SCRIPTOBS_STATUS != 0) and ($apftask.SCRIPTOBS_STATUS != 1) "
                    if APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                        apflog("Starting an instance of scriptobs for dynamic observing.",echo=True)
                        self.scriptobs = APF.startRobot()
                else:
                    apflog("Found Fixed list %s" % self.fixedList, echo=True)
                    apflog("Starting fixed list on line %d" % int(APF.ldone), echo=True)
                    APF.observe(str(self.fixedList))
                    APFTask.waitFor(self.task, True, timeout=10)
            else:
                if self.BV is None:
                    apflog("No B-V value at the moment", echo=True)
                    #self.BV = 0.028
                if self.VMAG is None:
                    apflog("No VMag value at the moment", echo=True)
                    #self.VMAG = None
                # We wish to observe either a starlist with intelligent ordering, or the dynamic scheduler
                apflog("Starting an instance of scriptobs for dynamic observing.",echo=True)
                self.scriptobs = APF.startRobot()
                            
                # Don't let the watcher run over the robot starting up
                APFTask.waitFor(self.task, True, timeout=10)
            
            return

       
            

        ###############################

        # Actual Watching loop
        apflog("Beginning observing process....",echo=True)
        APF.DMZero()
        haveobserved = False
        failstart = 0
        while self.signal:
            # Check on everything
            if datetime.now().strftime("%p") == 'AM':
                rising = True
                sunel_lim = SUNEL_ENDLIM
            else:
                rising = False
                sunel_lim = sunel_startlim()
            wind_vel = APF.wvel
            ripd, running = APF.findRobot()
            sunel = APF.sunel
            sunel.monitor()


            # if paused:
            #     apflog("Pausing because of apftask request",level='warn',echo=True)
            #     APFTask.waitfor(self.task, True, timeout=60)
            
            # Check and close for weather
            if APF.isOpen()[0] and not APF.openOK:
                closetime = datetime.now()
                APFTask.set(parent,suffix="MESSAGE",value="Closing for weather",wait=False)
                apflog("No longer ok to open.", echo=True)
                apflog("OPREASON: " + APF.checkapf["OPREASON"].read(), echo=True)
                apflog("WEATHER: " + APF.checkapf['WEATHER'].read(), echo=True)
                if running:
                    APF.killRobot(now=True)

                APF.close()
                APF.updateLastObs()

            # Check the slowdown factor to close for clouds
            if self.VMAG is not None and self.BV is not None and False:
                exp_cntrate = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, APF.ael,APF.avg_fwhm)
                try:
                    slow = exp_cntrate / APF.countrate
                    if slow < 0:
                        slow = 5
                        apflog("Countrate non-sensical %g" % APF.countrate, echo=True)
                except ZeroDivisionError:
                    apflog("No current countrate", echo=True)
                    slow = 5
                APFTask.set(parent,suffix="MESSAGE",value="FWHM = %.2f and slowdown %.2f" % (APF.avg_fwhm,slow),wait=False)
                if slow > 16:
                    # The slowdown is too high, we should close up and wait.
                    APFTask.set(parent,suffix="MESSAGE",value="Closing for clouds",wait=False)
                    apflog("Slowdown factor of %.2f is too high. Waiting 30 min to check again." % slow, echo=True)
                    APF.killRobot(now=True)
                    APF.close()
                    APFTask.waitfor(self.task, True, timeout=60*30)
                    self.VMAG=None
                    self.BV=None
                    APF.countrate = 0

            
            # If scriptobs is running and waiting for input, give it a target
            if running == True and float(sunel) < sunel_lim and APF.sop.read().strip() == 'Input':
                if self.fixedList is None:
                    APFTask.set(parent,suffix="MESSAGE",value="Calling getTarget",wait=False)
                    apflog("Scriptobs phase is input ( dynamic scheduler ), calling getTarget.")
                    getTarget()
                    apflog("Observing target")
                    APFTask.set(parent,suffix="MESSAGE",value="Observing Target",wait=False)
                    APFTask.waitfor(self.task, True, timeout=15)
                    haveobserved = True                    
                    if self.obsBstar:
                        self.obsBstar = False
                    try:
                        s=""
                        if self.obsBstar:
                            s="True"
                        APFTask.set(parent,suffix="VAR_3",value=s,wait=False)
                    except:
                        apflog("Error: Cannot communicate with apftask",level="error")


                elif self.smartObs == True:
                    apflog("Scriptobs phase is input ( smartlist ), calling getTarget.")
                    APFTask.set(parent,suffix="MESSAGE",value="Calling getTarget for a smartlist",wait=False)
                    getTarget()
                    APFTask.waitfor(self.task, True, timeout=15)
                    apflog("Observing target")
                    APFTask.set(parent,suffix="MESSAGE",value="Observing Target",wait=False)
                    haveobserved = True
                                                           
            # check last telescope focus
            lastfoc = APF.robot['FOCUSTEL_LAST_SUCCESS'].read(binary=True)
            if time.time() - lastfoc > FOCUSTIME and running and float(sunel) <= sunel_lim and haveobserved and APF.sop.read().strip() == 'Input':
                APFTask.set(parent,suffix="MESSAGE",value="More than %.1f hours since telescope focus, now focusing" % (FOCUSTIME/3600.),wait=False)
#                APF.focusTel()
                haveobserved = False
                
            # If the sun is rising and we are finishing an observation
            # Send scriptobs EOF. This will shut it down after the observation
            if float(sunel) >= sunel_lim and running == True:
                APFTask.set(parent,suffix="MESSAGE",value="Last call",wait=False)
                if self.scriptobs is None:
                    apflog("Robot claims to be running, but no self.scriptobs instance can be found. Instead calling killRobot().", echo=True)
                    APF.killRobot()
                else:
                    self.scriptobs.stdin.close()
                    APF.killRobot()
            
            # If the sun is rising and scriptobs has stopped, run closeup
            if float(sunel) > sunel_lim and running == False and rising == True:
                apflog("Closing due to sun elevation. Sunel = % 4.2f" % float(sunel), echo=True)
                APFTask.set(parent,suffix="MESSAGE",value="Closing, sun is rising",wait=False)
                if APF.isOpen()[0]:
                    msg = "APF is open, closing due to sun elevation = %4.2f" % float(sunel)
                    APF.close()
                else:
                    msg = "Telescope was already closed when sun got to %4.2f" % float(sunel)
                
                if APF.isOpen()[0]:
                    apflog("Error: Closeup did not succeed", level='error', echo=True)
                APF.updateLastObs()
                self.exitMessage = msg
                self.stop()

            # If we can open, try to set stuff up so the vent doors can be controlled by apfteq
            if APF.openOK and not rising and not APF.isOpen()[0]:
                APFTask.set(parent,suffix="MESSAGE",value="Powering up for APFTeq",wait=False)                    
                if APF.clearestop():
                    try:
                        APFLib.write(APF.dome['AZENABLE'],'enable',timeout=10)
                    except:
                        apflog("Error: Cannot enable AZ drive, exiting",level="error")
                        return
                    apf.setTeqMode('Evening')
                else:
                    apflog("Error: Cannot clear emergency stop, sleeping for 600 seconds",level="error")
                    APFTask.waitFor(parent,True,timeout=600)
                
            # Open at sunset
            sun_between_limits = float(sunel) < SUNEL_HOR and float(sunel) > sunel_lim 
            if not APF.isReadyForObserving()[0] and float(sunel) < SUNEL_HOR and float(sunel) > sunel_lim and APF.openOK and not rising:
                APFTask.set(parent,suffix="MESSAGE",value="Open at sunset",wait=False)                    
                success = opening( sunel, sunset=True)
                if success == False:
                    apflog("Error: Cannot open the dome",echo=True,level='error')
                    APF.close()
                    os._exit()
                
            # Open at night
            if not APF.isReadyForObserving()[0]  and float(sunel) < sunel_lim and APF.openOK:
                if not rising or (rising and float(sunel) < (sunel_lim - 5)):
                    APFTask.set(parent,suffix="MESSAGE",value="Open at night",wait=False)                    
                    success = opening( sunel)
                    if success == False:
                        apflog("Error: Cannot open the dome",echo=True,level='error')
                        APF.close()
                        os._exit()

            # Check for servo errors
            if APF.isOpen()[0] and APF.slew_allowed == False:
                APFTask.set(parent,suffix="MESSAGE",value="Servo failure.",wait=False)                    
                
                apflog("Error: APF is open, and slew_allowed is false. Likely an amplifier fault.", level="error", echo=True)
#                apflog("Forcing checkapf to close the dome. Heimdallr will then exit.", echo=True)
                chk_done = "$checkapf.MOVE_PERM == true"
#                APFLib.write(APF.dmtimer, 0)
                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=600)
                if not result and "DomeShutter" in APF.isOpen()[1]:
                    apflog("Error: After 10 min move permission did not return, and the dome is still open.", level='error', echo=True)

                APF.close(force=True)
                if not APF.power_down_telescope():
                    apflog("Error: Cannot reset telescope after servo failure",level="error", echo=True)
                    os._exit(1)
                
            # If we are open and scriptobs isn't running, start it up
            if APF.isReadyForObserving()[0] and not running and float(sunel) <= sunel_lim:
                APFTask.set(parent,suffix="MESSAGE",value="Starting scriptobs",wait=False)
                startScriptobs()
                if not APFTask.waitFor(self.task,True,expression="$apftask.SCRIPTOBS_STATUS == 'Running'",timeout=10):
                    if failstart % 11 == 0 and failstart > 0:
                        lvl = "error"
                    else:
                        lvl = "warn"
                    apflog("scriptobs is not running just after being started!", level=lvl, echo=True)
                    APFTask.set(parent,suffix="MESSAGE",value="scriptobs is not running just after being started!",wait=False)                    

                
            # Keep an eye on the deadman timer if we are open 
            if APF.isOpen()[0] and APF.dmtime <= DMLIM:
                APFTask.set(parent,suffix="MESSAGE",value="Reseting DM timer",wait=False)                    
                APF.DMReset()
#                apflog("The APF is open, the DM timer is clicking down, and scriptobs is %s." % ( str(running)),level="debug")

            if not APF.isOpen()[0] and not rising:
                APFTask.set(parent,suffix="MESSAGE",value="Waiting for sunset",wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
            if  APF.isOpen()[0] and float(sunel) > sunel_lim:
                APFTask.set(parent,suffix="MESSAGE",value="Waiting for the end of twilight",wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
            

    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)


if __name__ == '__main__':

    apflog("Starting Heimdallr...")

    # Parse the command line arguments
    opt = args()

    # Register the atexit function after parsing the command line arguments
    # This prevents printing the help followed by exited/failure message
    atexit.register (shutdown)

    # Log the Command line arguments
    apflog("Command Line Args:")
    od = vars(opt)
    for o in od:
        if od[o] is not None:
            apflog(o + ': ' +str(od[o]))


    if opt.test:
        debug = True
        parent = 'example'
        apflog("Heimdallr starting in test mode.")
    else:
        debug = False
        parent = 'master'


    if not opt.sheet:
        opt.sheet = "The Googledex"
    apftask = ktl.Service("apftask")        
    # Establish this as the only running master script ( Or example task in test mode )
    try:
        apflog("Attempting to establish apftask as %s" % parent)
        APFTask.establish(parent, os.getpid())
    except Exception as e:
        print e
        apflog("Task is already running with name %s." % parent, echo=True)
        sys.exit("Couldn't establish APFTask %s" % parent)
    else:
        # Set up monitoring of the current master phase
        phase = apftask("%s_PHASE" % parent)
        phase.monitor()

    # Set preliminary signal and tripwire conditions
    apflog("Setting APFTask signal and tripwire.")
    APFTask.set(parent, "SIGNAL", "TERM")
    APFTask.set(parent, "TRIPWIRE", "TASK_ABORT")

    control = apftask[parent + '_CONTROL']
    cw = functools.partial(control_watch,parent=parent)
    control.monitor()
    control.callback(cw)

    apflog("Master initiallizing APF monitors.", echo=True)

    # Aquire an instance of the APF class, which holds wrapper functions for controlling the telescope
    apf = ad.APF(task=parent, test=debug)
    APFTask.waitFor(parent, True, timeout=5)

    if apf.checkapf['USERKIND'].read(binary=True) != 3:
        apflog("checkapf not in robotic mode, exiting",level="error",echo=True)
        sys.exit()
    
    # Check to see if the instrument has been released
    # if not debug:
    #     if apf.checkapf['INSTRELE'].read().strip().lower() != 'yes':
    #         apflog("The instrument has not been released. Check that Observer Location has been submitted.", echo=True, level='error')
    #         sys.exit(1)
        

    # All the phase options that this script uses. This allows us to check if we exited out of the script early.
    possible_phases = ["ObsInfo", "Focus", "Cal-Pre", "Cal-Post", "Watching"]

    # If a command line phase was specified, use that.
    if opt.phase != None:
        APFTask.phase(parent, opt.phase)
        phase.poll()
        
    # If the phase isn't a valid option, (say the watchdog was run last)
    # then assume we are starting a fresh night and start from setting the observer information.
    apflog("Phase at start is: %s" % phase, echo=True)
    if str(phase).strip() not in possible_phases:
        apflog("Starting phase is not valid. Phase being set to ObsInfo", echo=True)
        APFTask.phase(parent, "ObsInfo")

    # Make sure that the command line arguments are respected.
    # Regardless of phase, if a name, obsnum, or reset was commanded, make sure we perform these operations.
    apflog("Setting scriptobs_lines_done=0")
    APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
    if not opt.fixed:
        APFTask.set(parent,"STARLIST","")
    else:
        if not os.path.exists(opt.fixed):
            errmsg = "starlist %s does not exist" % (opt.fixed)
            apflog(errmsg,level="error",echo=True)
            sys.exit(errmsg)
        APFTask.set(parent,"STARLIST",opt.fixed)

    if str(phase).strip() != "ObsInfo":
        if opt.obsnum:
            apflog("option -o specified. Setting UCAM OBSNUM to %d." % int(opt.obsnum)) 
            APFLib.write(apf.ucam["OBSNUM"], int(opt.obsnum))
        if opt.name:
            apflog("option -n specified. Setting UCAM Observer to %s." % opt.name)
            APFLib.write(apf.ucam["OBSERVER"], opt.name)
            APFLib.write(apf.ucam["OUTFILE"], opt.name)

    # Start the actual operations
    # Goes through 5 steps:
    # 1) Set the observer information
    # 2) Run focuscube
    # 3) Run calibrate ucsc pre
    # 4) Start the main watcher
    # 5) Run calibrate ucsc post
    # Specifying a phase jumps straight to that point, and continues from there.


    # 1) Setting the observer information.
    # Sets the Observation number, observer name, file name, and file directory
    if "ObsInfo" == str(phase).strip():
        apflog("Setting the task step to 0")
        APFTask.step(parent,0)
        if opt.obsnum == None:
            apflog("Figuring out what the observation number should be.",echo=False)
            obsNum = findObsNum(apf)
        else:
            obsNum = int(opt.obsnum)

        apflog("Using %s for obs number." % repr(obsNum),echo=True)
        apflog("Setting Observer Information", echo=True)
        if opt.name is None:
            opt.owner = 'Vogt'
            opt.name = 'ucsc'
            apf.setObserverInfo(num=obsNum, name='ucsc',owner=opt.owner)
        else:
            if opt.owner == None:
                opt.owner = opt.name
            apf.setObserverInfo(num=obsNum, name=opt.name, owner=opt.owner)
                
        apflog("Setting ObsInfo finished. Setting phase to Focus.")
        apflog("Setting SCRIPTOBS_LINES_DONE to 0")
        APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
        APFLib.write(apf.robot["MASTER_VAR_2"], time.time())
        APFTask.phase(parent, "Focus")
        apflog("Phase is now %s" % phase)

    # Run autofocus cube
    if "Focus" == str(phase).strip():
        apflog("Starting focusinstr script.", level='Info', echo=True)
        instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        while not instr_perm:
            apflog("Waiting for instrument permission to be true")
            APFTask.waitfor(parent,True,expression="$checkapf.INSTR_PERM = true",timeout=600)
            instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        result = apf.focus(user='ucsc')
        if not result:
            apflog("Focusinstr has failed. Observer is exiting.",level='error',echo=True)
            sys.exit(1)
        apflog("Focus has finished. Setting phase to Cal-Pre")
        apf.updateLastObs()

        APFTask.phase(parent, "Cal-Pre")
        apflog("Phase now %s" % phase)

    # Run pre calibrations
    if 'Cal-Pre' == str(phase).strip():
        try:
            APFTask.set(parent,suffix="VAR_3",value="True")
        except:
            apflog("Error: Cannot communicate with apftask",level="error")

        try:
            if opt.fixed == None:
                names,star_table,do_flags,stars = ds.parseGoogledex(sheetn=opt.sheet)
        except Exception as e:
            apflog("Error: Cannot download googledex?! %s" % (e),level="error")
            
        APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
        # try:
        #     APFLib.write("apfmot.DEWARFOCRAW", AVERAGE_INSTRFOC,timeout=60)
        #     APFLib.write("apftask.FOCUSINSTR_LASTFOCUS", AVERAGE_INSTRFOC,timeout=60)
        #     apflog("Moved instrument focus to %d" % (AVERAGE_INSTRFOC),echo=True)            
        # except:
        #     apflog("Cannot move instrument focus to %d" % (AVERAGE_INSTRFOC),level="error",echo=True)

        if apf.ucam['OUTFILE'].read() == 'ucsc':
            APFTask.set(parent,suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())
            

        apflog("Starting calibrate pre script.", level='Info', echo=True)
        instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        while not instr_perm:
            apflog("Waiting for instrument permission to be true")
            APFTask.waitfor(parent,True,expression="$checkapf.INSTR_PERM = true",timeout=600)
            instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)

        result = apf.calibrate(script = opt.calibrate, time = 'pre')
        apf.updateLastObs()
        if result == False:
            apflog("Calibrate Pre has failed. Trying again",level='warn',echo=True)
            result = apf.calibrate(script = opt.calibrate, time = 'pre')
            if not result:
                apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                sys.exit(2)
        apflog("Calibrate Pre has finished. Setting phase to Watching.")
        APFTask.phase(parent, "Watching")
        apflog("Phase is now %s" % phase)
        if apf.ucam['OUTFILE'].read() == 'ucsc':
            APFTask.set(parent,suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())


    # Start the main watcher thread
    master = Master(apf,user=opt.name,sheetn=opt.sheet,owner=opt.owner)
    try:
        if opt.name == "ucsc":
                names,star_table,do_flags,stars = ds.parseGoogledex(sheetn=opt.sheet)
    except Exception as e:
        apflog("Error: Cannot download googledex?!  %s" % (e),level="error")
    
    if 'Watching' == str(phase).strip():
        apflog("Starting the main watcher." ,echo=True)

        bstr = "%d,%d" % (opt.binning,opt.binning)
        apf.ucam['BINNING'].write(bstr)

        if opt.binning > 1:
            apfmon = ktl.Service('apfmon')
            d = time.time() + 22.*3600
            apfmon['BINNINGDIS'].write(d,binary=True)
            
        if opt.fixed != None:
            lastList = apf.robot["MASTER_STARLIST"].read()
            # This resets lines done if this is a new target list
            if opt.fixed != lastList:
                APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
                APFLib.write(apf.robot["MASTER_STARLIST"], opt.fixed)
            if opt.line != None:
                APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], int(opt.line))
                apflog("Starting star list %s at line %d" % (opt.fixed, int(opt.line)),echo=True)
            else:
                apflog("Starting star list %s" % opt.fixed,echo=True)
        else:
            apflog("Starting dynamic scheduler", echo=True)
        master.fixedList = opt.fixed
        master.smartObs = opt.smartObs
        master.task = parent
        master.windsheild = opt.windshield
        master.start()
    else:
        master.signal = False

    while master.signal:
        # Master is running, check for keyboard interupt
        try:
            currTime = datetime.now()
            # Check if it is after ~9:00AM.
            # If it is, something must be hung up, so lets force
            #  a closeup and run post cals. 
            #  Log that we force a close so we can look into why this happened.
            if currTime.hour == 9:
                # Its 9 AM. Lets closeup
                master.stop()
                apflog("Master was still running at 9AM. It was stopped and post calibrations will be attempted.", level='warn')
                break

            if debug:
                print 'Master is running.'
                print str(apf)
            APFTask.waitFor(parent, True, timeout=30)
        except KeyboardInterrupt:
            apflog("Heimdallr has been killed by user.", echo=True)
            master.stop()
            sys.exit()
        except:
            apflog("Heimdallr killed by unknown.", echo=True)
            master.stop()
            sys.exit()

    # Check if the master left us an exit message.
    # If so, something strange likely happened so log it.
    try:
        msg = master.exitMessage
    except AttributeError:
        pass
    else:
        apflog(msg, level='Info', echo=True)

    # We have finished taking data, and presumably it is the morning.
    apf.setTeqMode('Morning')


    # Keep a copy of observed_targets around for a bit just in case
    if os.path.exists(os.path.join(os.getcwd(),"observed_targets")):
        if opt.name == "ucsc":
            try:
                apflog("Updating the online googledex with the observed times", level='Info', echo=True)
                ds.update_googledex_lastobs(os.path.join(os.getcwd(),"observed_targets"),sheetn=master.sheetn)
            except:
                apflog("Error: Updating the online googledex has failed.", level="error")
        logpush(os.path.join(os.getcwd(),"observed_targets"))
        
    if os.path.exists(os.path.join(os.getcwd(),"robot.log")):
        logpush(os.path.join(os.getcwd(),"robot.log"))
        
    if master.nighttargetlog:
        try:
            logpush(master.nighttargetlogname)
        except:
            apflog("cannot roll %s" % (master.nighttargetlogname))

    # If there is a copy of the googledex laying around, remove it so it gets re-downloaded.
    try:
        os.remove(os.path.join(os.getcwd(),"googledex.dat"))
    except OSError:
        apflog("Note: There was no googledex save file to delete today.", echo=True)

    apfmon = ktl.Service('apfmon')
    if apfmon['BINNINGDIS'].read(binary=True) > 0:
        apfmon['BINNINGDIS'].write(0,binary=True)
              
    # Take morning calibrations
    APFTask.phase(parent, "Cal-Post")
    result = apf.calibrate(script=opt.calibrate, time='post')
    if not result:
        apflog("Calibrate Post has failed.", level='warn',echo=True)
        result = apf.calibrate(script=opt.calibrate, time='post')
        if not result:
            apflog("Calibrate Post has failed twice.", level='error',echo=True)
            APFTask.set(parent,suffix="MESSAGE",value="Calibrate Post failed twice",wait=False)
    apf.updateLastObs()

    bstr = "%d,%d" % (1,1)
    apf.ucam['BINNING'].write(bstr) 
            
    # Focus the instrument once more
    APFTask.phase(parent, "Focus")
    apflog("Running Focus Post", echo=True)
    result = apf.focus(user='ucsc')
    if not result:
        apflog("Focus cube post has failed", level='error', echo=True)
    else:
        apflog("Focus post has finished successfully.", echo=True)

    # We have done everything we needed to, so leave
    # the telescope in day mode to allow it to start thermalizing to the next night.
    apf.setTeqMode('Day')

    # Update the last observation number to account for the morning calibration shots.
    apf.updateLastObs()
    APFTask.set(parent,suffix="MESSAGE",value="Updating last observation number",wait=False)

    # All Done!
    APFTask.phase(parent, "Finished")

    success = True
    sys.exit()


