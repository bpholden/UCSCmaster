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

SUNEL_LIM = -10.0
SUNEL_HOR = -3.2

def control_watch(keyword,parent):
    if keyword['populated'] == False:
        return
    try:
        value = keyword['ascii']
        if value == "Abort":
            APF.log("Aborted by APFTask")
            sys.exit("Aborted by APFTask")
        elif value == "Pause":
            try:
                APFTask.set(parent,suffix='STATUS',value='PAUSED')
            except:
                APF.log("Failure to set STATUS in APFTask",level=error)
                sys.exit("Failure to communicate with APFTask")

        else:
            try:
                APFTask.set(parent,suffix='STATUS',value='Running')
            except:
                APF.log("Failure to set STATUS in APFTask",level=error)
                sys.exit("Failure to communicate with APFTask")

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
    parser = argparse.ArgumentParser(description="Set default options")
    parser.add_argument('-n', '--name', type=str, help='This values is used as the UCAM observer name, as well as the file prefix.')
    parser.add_argument('-o', '--obsnum', type=int, help='Sets the UCAM observation number to this integer value.')
    parser.add_argument('-p', '--phase', choices=p_c, help='Specify the starting phase of the watcher. Allows for skipping standard procedures.')
    parser.add_argument('-f', '--fixed', help='Specify a fixed target list to observe. File will be searched for relative to the current working directory.')
    parser.add_argument('-t', '--test', action='store_true', help="Start the watcher in test mode. No modification to telescope, instrument, or observer settings will be made.")
    parser.add_argument('-r', '--restart', action='store_true', default=False, help="Restart the specified fixed star list from the begining. This resets scriptobs_lines_done to 0.")
    parser.add_argument('-w', '--windshield', choices=w_c, default='auto', help="Turn windshielding on, off, or let the software decide based on the current average wind speed (Default is auto). Velocity > %.1f mph turns windshielding on." % (ad.WINDSHIELD_LIMIT))
    parser.add_argument('-c', '--calibrate', default='ucsc', type=str, help="Specify the calibrate script to use. Specify string to be used in calibrate 'arg' pre/post")
    parser.add_argument('-l', '--line', type=int, help="If a fixed starlist is given, starts the list at line N.")
    parser.add_argument('-s', '--smartObs', default=False, action='store_true', help="When specified with a fixed starlist, this option will pass the starlist to a selection algorithm to choose the optimal unobserved target, regardless of starlist ordering.")
    parser.add_argument('--sheet',default=None,help="Optional name for a Google spreadsheet")
    
    opt = parser.parse_args()
    return opt


def findObsNum(apf):

    obsNum = int(apf.robot["MASTER_LAST_OBS_UCSC"].read().strip())
    try:
        backup_obsNum = int(apf.robot["MASTER_VAR_2"].read(binary=True))
    except:
        backup_obsNum = -1
    if  backup_obsNum > obsNum:
        obsNum = backup_obsNum

    obsNum += 100 - (obsNum % 100)

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
    def __init__(self, apf, user='ucsc',sheetn="The Googledex"):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.APF = apf
        self.user = user
        self.name = 'Heimdallr'
        self.signal = True
        self.windshield = 'auto'
        self.scriptobs = None
        self.BV = None
        self.VMAG = None
        self.obsBstar = False
        self.fixedList = None
        self.sheetn = sheetn
        self.targetlogname = os.path.join(os.getcwd(),"targetlog.txt")
#        try:
#            self.targetlog = open(self.targetlogname,"w+")
#        except Exception, e:
        self.targetlog = None
#            apflog("cannot open %s: %s" % (self.targetlogname,e),level="error")

        self.nighttargetlogname = os.path.join(os.getcwd(),"nighttargetlog.txt")
#        try:
#            self.nighttargetlog = open(self.nighttargetlogname,"w+")
#        except Exception, e:
        self.nighttargetlog = None
#            apflog("cannot open %s: %s" % (self.nighttargetlogname,e),level="error")

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
                self.BV = 0.8 # use a default average
                slowdown = 5
            elif self.VMAG is None:
                apflog("Warning!: Ended up in getTarget() with no V Magnitude value, slowdown can't be computed.", echo=True)
                slowdown = 5
            elif APF.avg_fwhm < 1.0:
                apflog("Warning!: AVG_FWHM = %4.2f. By Odin's beard that seems low." % APF.avg_fwhm, echo=True)
                slowdown = 5
            else:
                apflog("Calculating expected counts")
                apflog("self.VMAG [%4.2f] - self.BV [%4.2f] - APF.ael [%4.2f]" % (self.VMAG, self.BV, APF.ael))
                exp_cnts_sec = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, APF.ael,APF.avg_fwhm)
                try:
                    slowdown = exp_cnts_sec / APF.countrate
                    if slowdown < 0:
                        slowdown = 1
                        apflog("Countrate non-sensical %g" % APF.countrate, echo=True)
                        # yes this happened.
                except ZeroDivisionError:
                    apflog("Current countrate was 0. Slowdown will be set to 5.", echo=True)
                    slowdown = 1
            apflog("getTarget(): slowdown factor = %4.2f" % slowdown, echo=True)
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
                target = ds.getNext(time.time(), seeing, slowdown, bstar=self.obsBstar, verbose=True,sheetn=self.sheetn)
            else:
                # Get the best target from the star list
                target = ds.smartList(self.fixedList, seeing, slowdown, APF.aaz, APF.ael)

            if target is None:
                apflog("No acceptable target was found. Since there does not seem to be anything to observe, Heimdallr will now shut down.", echo=True)
                # Send scriptobs EOF to finish execution - wouldn't want to leave a zombie scriptobs running
                self.scriptobs.stdin.close()
                APF.close()
                apf.countrate = -1.0
                # sleep for a half hour to see if the clouds blow by
                APFTask.waitfor(self.task, True, timeout=60*30)
                return
            else:
                apflog("Observing target: %s" % target['NAME'], echo=True)
                APFLib.write(self.APF.robot["SCRIPTOBS_AUTOFOC"], "robot_autofocus_enable")
                self.scriptobs.stdin.write(target["SCRIPTOBS"] + '\n')
            # Set the Vmag and B-V mag of the latest target
            self.VMAG = target["VMAG"]
            self.BV   = target["BV"]
            apflog("getTarget(): V=%.2f  B-V=%.2f Pri=%.2f " % (self.VMAG, self.BV, target["PRI"]))
            apflog("getTarget(): FWHM=%.2f  Slowdown=%.2f  Countrate=%.2f" % (APF.avg_fwhm, slowdown, APF.countrate))

            apflog("getTarget(): Target= %s" % target["NAME"])
            apflog("getTarget(): Counts=%.2f  EXPTime=%.2f  Nexp=%d" % (target["COUNTS"], target["EXP_TIME"], target["SCORE"]))


        def opening(sunset=False):
            when = "night"
            if sunset:
                when = "sunset"
                
            apflog("Running open at %s as sunel = %4.2f" % (when,el))
            (apfopen,what) =APF.isOpen()
            if apfopen:
                APF.DMReset()
            else:
                APF.DMZero()

            result = APF.openat(sunset=sunset)
            if not result:
                apflog("After two tries openatsunset hasn't successfully opened. \
                        Emailing for help and exiting.", level='error', echo=True)
                APF.close()
                os._exit(1)

            APF.DMReset()
            while float(APF.sunel) > SUNEL_LIM:
                chk_done = "$eostele.SUNEL < %f" % (SUNEL_LIM)
                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=60)
                APF.DMReset()
            return

        def startScriptobs():
            # Update the last obs file and hitlist if needed

            APF.updateLastObs()
            APF.updateWindshield(self.windshield)
            apflog("Starting an instance of scriptobs",echo=True)

            if self.fixedList is not None and self.smartObs == False:
                # We wish to observe a fixed target list, in it's original order
                tot = getTotalLines(self.fixedList)
                apflog("%d total starlist lines and %d lines done." % (tot, APF.ldone)) 
                if APF.ldone == tot and APF.user != "ucsc":
                    APF.close()
                    APF.updateLastObs()
                    self.exitMessage = "Fixed list is finished. Exiting the watcher."
                    self.stop()
                    # The fixed list has been completely observed so nothing left to do
                else:
                    apflog("Found Fixed list %s" % self.fixedList, echo=True)
                    apflog("Starting fixed list on line %d" % int(APF.ldone), echo=True)
                    APF.observe(str(self.fixedList), skip=int(APF.ldone))
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
                APFTask.waitFor(self.task, True, timeout=5)
            
            return

       
            

        ###############################

        # Actual Watching loop
        apflog("Beginning observing process....",echo=True)                
        while self.signal:
            # Check on everything
            if datetime.now().strftime("%p") == 'AM':
                rising = True
            else:
                rising = False
            wind_vel = APF.wvel
            ripd, running = APF.findRobot()
            el = float(APF.sunel)          

            # Check and close for weather
            if APF.isOpen()[0] and not APF.openOK:
                closetime = datetime.now()
                APFTask.set(parent,suffix="VAR_1",value="Closing for weather",wait=False)
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
                APFTask.set(parent,suffix="VAR_1",value="FWHM = %.2f and slowdown %.2f" % (APF.avg_fwhm,slow),wait=False)
                if slow > 16:
                    # The slowdown is too high, we should close up and wait.
                    APFTask.set(parent,suffix="VAR_1",value="Closing for clouds",wait=False)
                    apflog("Slowdown factor of %.2f is too high. Waiting 30 min to check again." % slow, echo=True)
                    APF.killRobot(now=True)
                    APF.close()
                    APFTask.waitfor(self.task, True, timeout=60*30)
                    self.VMAG=None
                    self.BV=None
                    APF.countrate = 0

            # If scriptobs is running and waiting for input, give it a target
            if running == True and el < SUNEL_LIM and APF.sop.read().strip() == 'Input':
                if self.fixedList is None:
                    APFTask.set(parent,suffix="VAR_1",value="Calling getTarget",wait=False)
                    apflog("Scriptobs phase is input ( dynamic scheduler ), calling getTarget.")
                    getTarget()
                    apflog("Observed target")
                    APFTask.set(parent,suffix="VAR_1",value="Observed Target",wait=False)
                    APFTask.waitfor(self.task, True, timeout=15)
                    
                    if self.obsBstar:
                        self.obsBstar = False
                    try:
                        s=""
                        if self.obsBstar:
                            s="True"
                        APFTask.set(parent,suffix="VAR_3",value=s,wait=False)
                    except:
                        apflog("Cannot communicate with apftask",level="error")


                elif self.smartObs == True:
                    apflog("Scriptobs phase is input ( smartlist ), calling getTarget.")
                    APFTask.set(parent,suffix="VAR_1",value="Calling getTarget for a smartlist",wait=False)
                    getTarget()
                    APFTask.waitfor(self.task, True, timeout=15)
                    
            # If the sun is rising and we are finishing an observation
            # Send scriptobs EOF. This will shut it down after the observation
            if el >= SUNEL_LIM and running == True:
                APFTask.set(parent,suffix="VAR_1",value="Last call",wait=False)
                if self.scriptobs is None:
                    apflog("Robot claims to be running, but no self.scriptobs instance can be found. Instead calling killRobot().", echo=True)
                    APF.killRobot()
                else:
                    self.scriptobs.stdin.close()
            
            # If the sun is rising and scriptobs has stopped, run closeup
            if el > SUNEL_LIM and running == False and rising == True:
                apflog("Closing due to sun elevation. El = % 4.2f" % el, echo=True)
                APFTask.set(parent,suffix="VAR_1",value="Closing, sun is rising",wait=False)
                if APF.isOpen()[0]:
                    msg = "APF is open, closing due to sun elevation = %4.2f" % el
                    APF.close()
                else:
                    msg = "Telescope was already closed when sun got to %4.2f" % el
                
                if APF.isOpen()[0]:
                    apflog("Closeup did not succeed", level='Error', echo=True)
                APF.updateLastObs()
                self.exitMessage = msg
                self.stop()

            # If we can open, try to set stuff up so the vent doors can be controlled by apfteq
            if APF.openOK and not rising and not APF.isOpen()[0]:
                APFTask.set(parent,suffix="VAR_1",value="Powering up for APFTeq",wait=False)                    
                APF.clearestop()
                try:
                    APFLib.write(APF.dome['AZENABLE'],'enable',timeout=10)
                except:
                    apflog("cannot enable AZ drive",level="error")
                apf.setTeqMode('Evening')
                
            # Open at sunset
            sun_between_limits = el < SUNEL_HOR and el > SUNEL_LIM 
            if not APF.isReadyForObserving()[0] and el < SUNEL_HOR and el > SUNEL_LIM and APF.openOK and not rising:
                APFTask.set(parent,suffix="VAR_1",value="Open at sunset",wait=False)                    
                opening(sunset=True)
                
            # Open at night
            if not APF.isReadyForObserving()[0]  and el < SUNEL_LIM and APF.openOK:
                APFTask.set(parent,suffix="VAR_1",value="Open at night",wait=False)                    
                opening()

            # Check for servo errors
            if APF.isOpen()[0] and APF.slew_allowed == False:
                APFTask.set(parent,suffix="VAR_1",value="Servo failure.",wait=False)                    
                
                apflog("APF is open, and slew_allowed is false. Likely an amplifier fault.", level="error", echo=True)
                apflog("Forcing checkapf to close the dome. Heimdallr will then exit.", echo=True)
                chk_done = "$checkapf.MOVE_PERM == true"
                APFLib.write(APF.dmtimer, 0)
                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=600)
                if not result and "DomeShutter" in APF.isOpen()[1]:
                    apflog("After 10 min move permission did not return, and the dome is still open.", level='error', echo=True)

                APF.closeup(force=True)
                if not APF.power_down_telescope():
                    apflog("Cannot reset telescope after servo failure",level="error", echo=True)
                    os._exit(1)
                

            # If we are open and scriptobs isn't running, start it up
            if APF.isReadyForObserving()[0] and not running and el <= SUNEL_LIM:
                APFTask.set(parent,suffix="VAR_1",value="Starting scriptobs",wait=False)                    
                startScriptobs()
                
                    
                
            # Keep an eye on the deadman timer if we are open 
            if APF.isOpen()[0] and APF.dmtime <= 1140:
                APFTask.set(parent,suffix="VAR_1",value="Reseting DM timer",wait=False)                    
                APF.DMReset()
                apflog("The APF is open, the DM timer is clicking down, and scriptobs is %s." % ( str(running)))

            if not APF.isOpen()[0]:
                APFTask.set(parent,suffix="VAR_1",value="Waiting for sunset",wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
            if  APF.isOpen()[0] and el > SUNEL_LIM:
                APFTask.set(parent,suffix="VAR_1",value="Waiting for the end of twilight",wait=False)
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
        apftask = ktl.Service("apftask")
        phase = apftask("%s_PHASE" % parent)
        phase.monitor()

    # Set preliminary signal and tripwire conditions
    apflog("Setting APFTask signal and tripwire.")
    APFTask.set(parent, "SIGNAL", "TERM")
    APFTask.set(parent, "TRIPWIRE", "TASK_ABORT")

    control = apftask[parent + '_CONTROL']
    cw = functools.partial(control_watch,parent=parent)
    control.callback(cw)
    control.read()
    control.monitor()

    apflog("Master initiallizing APF monitors.", echo=True)

    # Aquire an instance of the APF class, which holds wrapper functions for controlling the telescope
    apf = ad.APF(task=parent, test=debug)
    APFTask.waitFor(parent, True, timeout=5)

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
    if opt.restart:
        apflog("Restart specified. Setting scriptobs_lines_done=0")
        APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
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

        print "Welcome! I think the starting observation number should be:"
        print repr(obsNum)
        print ''
        print "If you believe this number is an error, please enter the correct number within the next 15 seconds..."
        rlist, _, _ = select([sys.stdin], [], [], 15)
        if rlist:
            s = sys.stdin.readline()
            while True:
                try:
                    v = int(s.strip())
                except ValueError:
                    print "Couldn't turn %s into an integer." % s
                else:
                    break
                s = raw_input("Enter Obs Number:")
            
            obsNum = v

        apflog("Using %s for obs number." % repr(obsNum),echo=True)
        apflog("Setting Observer Information", echo=True)
        if opt.name is None:
            apf.setObserverInfo(num=obsNum, name='ucsc')
        else:
            apf.setObserverInfo(num=obsNum, name=opt.name)
        apflog("Setting ObsInfo finished. Setting phase to Focus.")
        apflog("Setting SCRIPTOBS_LINES_DONE to 0")
        APFLib.write(apf.robot["SCRIPTOBS_LINES_DONE"], 0)
        APFTask.phase(parent, "Focus")
        apflog("Phase is now %s" % phase)

    # Run autofocus cube
    if "Focus" == str(phase).strip():
        apflog("Starting focusinstr script.", level='Info', echo=True)
        instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        while not instr_perm:
            apflog("Waiting for instrument permission to be true")
            APFTask.waitfor(parent,True,expression="$checkapf.INSTR_PERM = true",timeout=600)
        result = apf.focus(user='ucsc')
        if not result:
            apflog("Focusinstr has failed. Observer is exiting.",level='error',echo=True)
            sys.exit(1)
        apflog("Focus has finished. Setting phase to Cal-Pre")
        APFTask.phase(parent, "Cal-Pre")
        apflog("Phase now %s" % phase)

    # Run pre calibrations
    if 'Cal-Pre' == str(phase).strip():
        try:
            APFTask.set(parent,suffix="VAR_3",value="True")
        except:
            apflog("Cannot communicate with apftask",level="Error")

        try:
            if opt.fixed == None:
                (names,) = ds.parseGoogledex(sheetn=opt.sheetn)
        except:
            apflog("Cannot download googledex?!",level="Error")

        apflog("Starting calibrate pre script.", level='Info', echo=True)
        instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        while not instr_perm:
            apflog("Waiting for instrument permission to be true")
            APFTask.waitfor(parent,True,expression="$checkapf.INSTR_PERM = true",timeout=600)
        result = apf.calibrate(script = opt.calibrate, time = 'pre')
        if not result:
            apflog("Calibrate Pre has failed. Observer is exiting.",level='error',echo=True)
            sys.exit(2)
        apflog("Calibrate Pre has finished. Setting phase to Watching.")
        APFTask.phase(parent, "Watching")
        apflog("Phase is now %s" % phase)


    # Start the main watcher thread
    master = Master(apf,user=opt.name)
    try:
        if opt.fixed == None:
                (names,) = ds.parseGoogledex(sheetn=opt.sheetn)
    except:
        apflog("Cannot download googledex?!",level="Error")
    
    if 'Watching' == str(phase).strip():
        apflog("Starting the main watcher." ,echo=True)
        if master.nighttargetlog != None:
            try:
                apfguide = ktl.Service('apfguide')
                midpt = apfguide['MIDPTFIN']
                try:
                    midpt.callback(functools.partial(ad.midptmon,outputfile=master.nighttargetlog,permoutfile=master.targetlog))
                except Exception, e:
                    apflog("Cannot setup midpoint monitor: %s" % (e),level="warn")
            except Exception, e:
                apflog("Cannot setup midpoint keyword: %s" % (e),level="warn")
            
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
                apflog("Master was still running at 9AM. It was stopped and post calibrations will be attempted.", level='Warn')
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

    if master.fixedList is None:
        try:
            if os.path.exists(master.nighttargetlogname):
                master.nighttargetlog.close()
                master.targetlog.close()
                ds.update_googledex_lastobs(master.nighttargetlogname,sheetn=master.sheetn)
            else:
                ds.update_googledex_lastobs(os.path.join(os.getcwd(),"observed_targets"),sheetn=master.sheetn)
        except:
            apflog("Updating the online googledex has failed.", level="Error")

    # Keep a copy of observed_targets around for a bit just in case
    logpush(os.path.join(os.getcwd(),"observed_targets"))
    logpush(os.path.join(os.getcwd(),"robot.log"))
    if master.nighttargetlog:
        try:
            logpush(master.nighttargetlogname)
        except:
            apflog("cannot roll %s" % (master.nighttargetlogname))

    # If there is a copy of the googledex laying around, remove it so it gets re-downloaded tomorrow.
    try:
        os.remove(os.path.join(os.getcwd(),"googledex.dat"))
    except OSError:
        apflog("Note: There was no googledex save file to delete today.", echo=True)

              
    # Take morning calibrations
    APFTask.phase(parent, "Cal-Post")
    result = apf.calibrate(script=opt.calibrate, time='post')
    if not result:
        apflog("Calibrate Post has failed.", level='error',echo=True)

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

    # All Done!
    APFTask.phase(parent, "Finished")

    success = True
    sys.exit()


