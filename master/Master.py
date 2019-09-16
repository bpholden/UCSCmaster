from datetime import datetime, timedelta
import os
import os.path
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

try:
    import ktl
    import APF as APFLib
    import APFTask
except:
    pass

import APFControl
from apflog import *
import UCSCScheduler_V2 as ds
from x_gaussslit import *
import ExposureCalculations
import ParseGoogledex
import SchedulerConsts

SUNEL_ENDLIM = -10.0
SUNEL_STARTLIM = -9.0
SUNEL_HOR = -3.2
AVERAGE_INSTRFOC = 8522
DMLIM = 1140
FOCUSTIME = 3600.

class Master(threading.Thread):
    def __init__(self, apf, opt,totTemps=4):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.APF = apf
        self.task = 'master'
        self.user = opt.name
        self.owner = opt.owner
        self.name = 'Heimdallr'
        self.signal = True
        self.windshield = opt.windshield
        self.scriptobs = None
        self.BV = None
        self.VMAG = None
        self.decker = "W"
        self.obsBstar = True
        self.lastObsSuccess = True
        self.lastObsFinished = True
        self.fixedList = opt.fixed
        self.sheetn = opt.sheet
        self.targetlogname = os.path.join(os.getcwd(),"targetlog.txt")
        self.targetlog = None
        self.starttime = opt.start
        self.raster = opt.raster
        self.debug = opt.test
        self.doTemp = True
        self.nTemps = 0
        self.focval = 0
        self.totTemps = totTemps        

        self.target = None
        
        self.apftask = ktl.Service('apftask')
        self.lineresult = self.apftask['SCRIPTOBS_LINE_RESULT']
        self.lineresult.monitor()
        self.observed = self.apftask['SCRIPTOBS_OBSERVED']
        self.observed.monitor()
        
        self.nighttargetlogname = os.path.join(os.getcwd(),"nighttargetlog.txt")
        self.nighttargetlog = None



    def set_autofocval(self):
        """ Master.set_autofocval()
            tests when the last time the telescope was focused, if more than FOCUSTIME enable focus check
        """
        APF = self.APF
        # check last telescope focus
        lastfoc = APF.robot['FOCUSTEL_LAST_SUCCESS'].read(binary=True)
        current_val = APF.autofoc.read()
        rising = APF.rising
        cur_sunel = APF.sunel.read(binary=True)
        too_close = rising and (cur_sunel > -20)
        if time.time() - lastfoc > FOCUSTIME:
            if current_val != "robot_autofocus_enable" and not too_close:
                APF.autofoc.write("robot_autofocus_enable")
                self.focval=1
                APFTask.set(parent, suffix="MESSAGE", value="More than %.1f hours since telescope focus" % (FOCUSTIME/3600.), wait=False)            
        else:
            if current_val == "robot_autofocus_enable":
                APF.autofoc.write("robot_autofocus_disable")
                self.focval=0


    def checkObsSuccess(self):
        """ Master.checkObsSuccess() 
            checks the value of SCRIPTOBS_LINE_RESULT to see if the last observation suceeded.
        """
        retval = False
        
        if self.lineresult.read(binary=True) == 3:
            retval = True
        return retval

    def checkObsFinished(self):
        """ Master.checkObsFinished() 
            checks the value of SCRIPTOBS_LINE to see if we are on the last line of the block
            checks SCRIPTOBS_LINE_RESULT and SCRIPTOBS_OBSERVED to see if the last line is done
        """
        retval = False

        mtch = re.search("end\Z",apf.line.read())
        if apf.ldone.read(binary=True) == 0 or mtch:
            retval = True
        return retval


    def checkBstar(self,haveobserved):
        """ Master.obsBstar(haveobserved) 
            if observing has begun, and the last observation was a success, set Master.obsBstar to false, writes master_obsbstar to
            the current value of obsBstar
            The variable OBSBSTAR still overrides
        """
        self.obsBstar = ktl.read('apftask','MASTER_OBSBSTAR',binary=True)
        
        if haveobserved and self.lastObsSuccess:
            self.obsBstar = False
            try:
                ktl.write('apftask','MASTER_OBSBSTAR',self.obsBstar,binary=True)
            except Exception, e:
                apflog("Error: Cannot communicate with apftask: %s" % (e),level="error")
            
    def readStarlistFile(filename):
        tot = 0
        self.target = dict()
        self.target["SCRIPTOBS"] = []
        with open(filename, 'r') as f:
            for line in f:
                sline = line.strip()
                if sline == '':
                    continue
                elif sline[0] == '#':
                    continue
                else:
                    tot += 1
                    self.target["SCRIPTOBS"].append(sline)
        return tot

    def shouldStartList(self):
        """ Master.shouldStartList()
            should we start a fixed observing list or not? true if start time is None or if w/in + 1 hour - 0.5 hours of start time
        """
        if self.starttime == None:
            return True
        ct = time.time()
        if ct - self.starttime < 3600 or self.starttime - ct < 1800:
            return True
        return False


    ####
    # run is the main event loop, for historical reasons it has its own functions that are local in scope
    ####

    
    def run(self):
        """ Master.run() - runs the observing
        """
        APF = self.APF

        def calcSlowdown():
        
            if self.BV is None:
                apflog("Warning!: Ended up in getTarget() with no B Magnitude value, slowdown can't be computed.", echo=True)
                self.BV = 0.6 # use a default average

            if self.VMAG is None:
                apflog("Warning!: Ended up in getTarget() with no V Magnitude value, slowdown can't be computed.", echo=True)
                return 5
                
            if APF.avg_fwhm < 1.0:
                apflog("Warning!: AVG_FWHM = %4.2f. By Odin's beard that seems low." % APF.avg_fwhm, echo=True)
                return 5
            slowdown = 1
            apflog("Calculating expected counts")
            apflog("self.VMAG [%4.2f] - self.BV [%4.2f] - APF.ael [%4.2f]" % (self.VMAG, self.BV, APF.ael))
            exp_cnts_sec = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, APF.ael, APF.avg_fwhm, self.decker)
            try:
                if APF.countrate <= 0:
                    try:
                        APF.countrate = APF.ccountrate
                    except:
                        APF.countrate = -1.0
                if APF.countrate*10 <  APF.ccountrate:
                    APF.countrate = APF.ccountrate
                slowdown = exp_cnts_sec / APF.countrate
                if slowdown < 0:
                    slowdown = 1
                    apflog("Countrate non-sensical %g" % (APF.countrate), echo=True, level='warn')
                    APF.counts.monitor(start=False)
                    APF.counts.monitor(start=True)
                    APF.counts.callback(APFControl.countmon)
                    # yes this happened.
                if slowdown < SchedulerConsts.SLOWDOWN_MIN:
                    slowdown = SchedulerConsts.SLOWDOWN_MIN
                    apflog("slowdown too low, countrate= %g" % (APF.countrate), echo=True, level='debug')
                    # yes this happened.
                if slowdown > SchedulerConsts.SLOWDOWN_MAX:
                    slowdown = SchedulerConsts.SLOWDOWN_MAX
                    apflog("slowdown too high, countrate= %g" % (APF.countrate), echo=True, level='debug')
            except ZeroDivisionError:
                apflog("Current countrate was 0. Slowdown will be set to 1.", echo=True)
                slowdown = 1

            return slowdown
        
        # This is called when an observation finishes, and selects the next target
        def getTarget():
            APFLib.write(APF.ucam["RECORD"], "Yes") # safe / sorry

            if APF.nerase != 2:
                APF.killRobot()
                rv = APF.ucam_reboot()
                self.scriptobs = APF.startRobot()
            
            if self.target is not None and 'SCRIPTOBS' in self.target.keys():
                if len(self.target["SCRIPTOBS"]) > 0:
                    # just keep going with last block
                    apflog("getTarget(): Going through remaining target queue.",echo=True)
                    try:
                        curstr = self.target["SCRIPTOBS"].pop() + '\n'
                        self.scriptobs.stdin.write(curstr)
                    except Exception, e:
                        apflog("Failure in getTarget: %s" % (e),level='error',echo=True)
                        pass
                    return
            if self.checkObsFinished():
                apflog("getTarget(): Scriptobs phase is input, determining next target.",echo=True)
            else:
                apflog("getTarget(): Not at end of block but out of targets.",echo=True)


            try:
                self.obsBstar = ktl.read("apftask", "MASTER_OBSBSTAR",binary=True)
                apflog("getTarget(): Setting obsBstar to %s" % (str(self.obsBstar)),echo=True)

                
            except Exception, e:
                apflog("getTarget(): Cannot change obsBstar: %s" % (e),level='error',echo=True)
                self.obsBstar = True

            if self.obsBstar:
                APF.autofoc.write("robot_autofocus_enable")
                
            if self.scriptobs is None:
                apflog("Called getTarget, but there is not instance of scriptobs associated with Heimdallr. This is an error condition.", echo=True)
                return 
            
            # Calculate the slowdown factor.
            slowdown = calcSlowdown()
            apflog("getTarget(): slowdown factor = %4.2f" % slowdown, echo=True)
            APFLib.write(apf.robot["MASTER_SLOWDOWN"], slowdown)
            apflog("getTarget(): countrate = %.2f, ccountrate = %.2f" % (APF.countrate, APF.ccountrate))

            # Check for a valid seeing measurment. If there isn't one, use a default
            if APF.avg_fwhm == 0.:
                apflog("getTarget(): Warning AVG_FWHM is 0. A default value of 15 will be used in its place.", echo=True)
                seeing = 15
            else:
                seeing = float(APF.avg_fwhm)
                apflog("getTarget(): Current AVG_FWHM = %4.2f" % seeing)
            
            self.target = ds.getNext(time.time(), seeing, slowdown, bstar=self.obsBstar,sheetns=self.sheetn, owner=self.owner, template=self.doTemp,focval=self.focval)

            self.set_autofocval()
            if self.target is None:
                apflog("No acceptable target was found. Since there does not seem to be anything to observe, Heimdallr will now shut down.", echo=True)
                # Send scriptobs EOF to finish execution - wouldn't want to leave a zombie scriptobs running
                self.scriptobs.stdin.close()
                APF.close()
                if self.fixedList is None:
                    APFLib.write(apf.ldone, 0)
                apf.countrate = -1.0
                # sleep for a half hour to see if the clouds blow by
                APFTask.waitfor(self.task, True, timeout=60*30)
                return
            else:
                apflog("Observing target: %s" % self.target['NAME'], echo=True)
                APFTask.set(parent, suffix="MESSAGE", value="Observing target: %s"  % self.target['NAME'], wait=False)
                
                self.scriptobs.stdin.write(self.target["SCRIPTOBS"].pop() + '\n')
                
            # Set the Vmag and B-V mag of the latest target
            self.VMAG = self.target["VMAG"]
            self.BV   = self.target["BV"]
            self.decker = self.target["DECKER"]
            istemp = str(self.target['isTemp'])
            
            apflog("getTarget(): V=%.2f  B-V=%.2f Pri=%.2f " % (self.VMAG, self.BV, self.target["PRI"]))
            apflog("getTarget(): FWHM=%.2f  Slowdown=%.2f  Countrate=%.2f" % (APF.avg_fwhm, slowdown, APF.countrate))

            apflog("getTarget(): Target= %s Temp=%s" % (self.target["NAME"], istemp))
            apflog("getTarget(): Counts=%.2f  EXPTime=%.2f  Nexp=%d" % (self.target["COUNTS"], self.target["EXP_TIME"], self.target["NEXP"]))
            if self.target['isTemp']:
                self.nTemps += 1
                if self.nTemps >= self.totTemps:
                    self.doTemp = False

            return 

        # opens the dome & telescope, if sunset is True calls open at sunset, else open at night
        def opening(sunel,sunset=False):
            when = "night"
            if sunset:
                when = "sunset"

            result = APF.ucam_status()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
                os._exit()

            result = APF.enable_obs_inst()
            if result == False:
                apflog("Cannot enable instrument", level='warn', echo=True)
                result = APF.enable_obs_inst()
                if not result:
                    apflog("Error: cannot enable instrument twice.", level='alert', echo=True)
                    return result
                
            apflog("Running open at %s as sunel = %4.2f" % (when, float(sunel)),echo=True)
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
                        return result

                    
            if datetime.now().strftime("%p") == 'PM':
                setting = True
            else:
                setting = False
            APF.DMReset()

            return result


        # closing
        def closing(force=False):
            if running:
                APF.killRobot(now=True)

            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

            rv = APF.close(force=force)
            if rv:
                return
            if APF.power_down_telescope() is False:
                apflog("Error: Cannot close and power off telescope ", level="alert", echo=True)
                
            rv = APF.disable_inst()
            
            return

        def checkTelState():
            slewing     = '$eostele.AZSSTATE == Slewing  or  $eostele.ELSSTATE == Slewing'
            tracking    = '$eostele.AZSSTATE == Tracking and $eostele.ELSSTATE == Tracking'
            istracking = ktl.Expression(tracking)
            isslewing = ktl.Expression(slewing)

            if istracking.evaluate() or isslewing.evaluate():
                rv = True
            else:
                rv = False
            return rv

        def startTelescope():
            '''This starts up the telescope if the Az drive is disabled or the E-Stop State is True
            If the telescope is just disabled, the start up procedure for a new version of scriptobs should clear that state.
            '''
            rv = False

            eosdome = ktl.Service('eosdome')
            isenabled = eosdome['AZDRVENA'].read(binary=True)
            isstopped = eosdome['ESTOPST'].read(binary=True)
            fullstop = eosdome['SWESTOP'].read(binary=True)
            if fullstop:
                rv = False
                # cannot start the telescope
            else:
                # we can!
                if isstopped:
                    eosdome['ESTOPCMD'].write('ResetEStop')
                if isenabled is False:
                    eosdome['AZENABLE'].write('Enable')
                isenabled = eosdome['AZDRVENA'].read(binary=True)
                isstopped = eosdome['ESTOPST'].read(binary=True)
                if isenabled and isstopped is False:
                    rv = True
                else:
                    rv = False
            
            return rv
        
        # starts an instance of scriptobs 
        def startScriptobs():
            # Update the last obs file and hitlist if needed

            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

            APF.updateWindshield(self.windshield)
            ripd, running = APF.findRobot()
            if running:
                apflog("Scriptobs is already running yet startScriptobs was called",level="warn",echo=True)
                return
            message = apf.message.read()
            mtch = re.search("ERR/UCAM",message)
            if mtch:
                # uh oh
                apflog("scriptobs has failed post UCAM recovery",level="error",echo=True)
                # reboot warsaw
                APF.ucam_restart()
            
            apflog("Starting an instance of scriptobs",echo=True)
            if self.fixedList is not None and self.shouldStartList():
                # We wish to observe a fixed target list, in it's original order
                if not os.path.exists(self.fixedList):
                    apflog("Error: starlist %s does not exist" % (self.fixedList), level="error")
                    self.fixedList=None
                    self.starttime = None
                # this reads in the list and appends it to self.target so getTarget will automatically use it
                tot = self.readStarlistFile(self.fixedList)
                if tot == 0:
                    apflog("Error: starlist %s is empty" % (self.fixedList), level="error")
                    self.fixedList=None
                    self.starttime = None
                    self.target=None
                else:
                    apflog("%d total starlist lines and %d lines done." % (tot, APF.ldone)) 
                if APF.ldone == tot :
                    self.fixedList = None
                    self.starttime = None
                    self.target = None
                    if not APF.test:
                        APFTask.set(self.task,suffix="STARLIST",value="")
                    apflog("Finished fixed list on line %d, will start dynamic scheduler" % int(APF.ldone), echo=True)
                else:
                    apflog("Found Fixed list %s" % self.fixedList, echo=True)
                    apflog("Starting fixed list on line %d" % int(APF.ldone), echo=True)
            else:
                if self.BV is None:
                    apflog("No B-V value at the moment", echo=True)
                    #self.BV = 0.028
                if self.VMAG is None:
                    apflog("No VMag value at the moment", echo=True)
                    #self.VMAG = None
                # We wish to observe with the dynamic scheduler
            ripd, running = APF.findRobot()
            if running is False:
                apflog("Starting an instance of scriptobs for dynamic observing.", echo=True)
                self.scriptobs = APF.startRobot()
                # Don't let the watcher run over the robot starting up
                APFTask.waitFor(self.task, True, timeout=10)
            
            return

       
            

        ###############################

        # Actual Watching loop
        apflog("Beginning observing process....", echo=True)
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
                sunel_lim = SUNEL_STARTLIM
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
                APFTask.set(parent, suffix="MESSAGE", value="Closing for weather", wait=False)
                apflog("No longer ok to open.", echo=True)
                apflog("OPREASON: " + APF.checkapf["OPREASON"].read(), echo=True)
                apflog("WEATHER: " + APF.checkapf['WEATHER'].read(), echo=True)
                closing()

            # Check the slowdown factor to close for clouds
            if self.VMAG is not None and self.BV is not None and False:
                exp_cntrate = ExposureCalculations.getEXPMeter_Rate(self.VMAG, self.BV, APF.ael, APF.avg_fwhm,["W"])
                try:
                    slow = exp_cntrate / APF.countrate
                    if slow < 0:
                        slow = 5
                        apflog("Countrate non-sensical %g" % APF.countrate, echo=True)
                except ZeroDivisionError:
                    apflog("No current countrate", echo=True)
                    slow = 5
                APFTask.set(parent, suffix="MESSAGE",value="FWHM = %.2f and slowdown %.2f" % (APF.avg_fwhm,slow), wait=False)
                if slow > 16:
                    # The slowdown is too high, we should close up and wait.
                    APFTask.set(parent, suffix="MESSAGE",value="Closing for clouds", wait=False)
                    apflog("Slowdown factor of %.2f is too high. Waiting 30 min to check again." % slow, echo=True)
                    closing()
                    APFTask.waitfor(self.task, True, timeout=60*30)
                    self.VMAG=None
                    self.BV=None
                    APF.countrate = 0

            
            # If scriptobs is running and waiting for input, give it a target
            if running == True and float(sunel) < sunel_lim and APF.sop.read().strip() == 'Input':
                if self.fixedList is None or self.shouldStartList() == False:
                    self.lastObsSuccess = self.checkObsSuccess()
                    self.obsBstar = self.checkBstar(haveobserved)
                    
                    APFTask.set(parent, suffix="MESSAGE", value="Calling getTarget", wait=False)
                    apflog("Scriptobs phase is input ( dynamic scheduler ), calling getTarget.")
                    getTarget()
                    APFTask.waitfor(self.task, True, timeout=15)
                    
                    haveobserved = True                    
                elif self.starttime != None and self.shouldStartList():
                    APF.killRobot()

            # check last telescope focus
            if running and  float(sunel) <= sunel_lim:
                self.set_autofocval()
                
            # If the sun is rising and we are finishing an observation
            # Send scriptobs EOF. This will shut it down after the observation
            if float(sunel) >= sunel_lim and running == True:
                APFTask.set(parent, suffix="MESSAGE", value="Last call", wait=False)
                if self.scriptobs is None:
                    apflog("Robot claims to be running, but no self.scriptobs instance can be found. Instead calling killRobot().", echo=True)
                    APF.killRobot()
                else:
                    self.scriptobs.stdin.close()
                    APF.killRobot()
            
            # If the sun is rising and scriptobs has stopped, run closeup
            if float(sunel) > sunel_lim and running == False and rising == True:
                apflog("Closing due to sun elevation. Sunel = % 4.2f" % float(sunel), echo=True)
                APFTask.set(parent, suffix="MESSAGE", value="Closing, sun is rising", wait=False)
                if APF.isOpen()[0]:
                    msg = "APF is open, closing due to sun elevation = %4.2f" % float(sunel)
                    closing()
                else:
                    msg = "Telescope was already closed when sun got to %4.2f" % float(sunel)
                
                if APF.isOpen()[0]:
                    apflog("Error: Closeup did not succeed", level='error', echo=True)

                self.exitMessage = msg
                self.stop()

            # If we can open, try to set stuff up so the vent doors can be controlled by apfteq
            if APF.openOK and not rising and not APF.isOpen()[0]:
                APFTask.set(parent, suffix="MESSAGE", value="Powering up for APFTeq", wait=False)                    
                if APF.clearestop():
                    try:
                        APFLib.write(APF.dome['AZENABLE'], 'enable', timeout=10)
                    except:
                        apflog("Error: Cannot enable AZ drive", level="error")

                    apf.setTeqMode('Evening')
                    vent_open = "$eosdome.VD4STATE = VENT_OPENED"
                    result = APFTask.waitfor(self.task, True, expression=vent_open, timeout=180)
                    if result:
                        try:
                            APFLib.write(APF.dome['AZENABLE'], 'disable', timeout=10)
                        except:
                            apflog("Error: Cannot disable AZ drive", level="error",echo=True)

                    else:
                        apflog("Error: Vent doors did not open, is apfteq and eosdome running correctly?", level='error',echo=True)
                else:
                    apflog("Error: Cannot clear emergency stop, sleeping for 600 seconds", level="error")
                    APFTask.waitFor(parent, True, timeout=600)
                
            # Open 

            if not APF.isReadyForObserving()[0] and float(sunel) < SUNEL_HOR and APF.openOK:
                if float(sunel) > sunel_lim and not rising:
                    APFTask.set(parent, suffix="MESSAGE", value="Open at sunset", wait=False)                    
                    success = opening(sunel, sunset=True)
                    if success is False:
                        if APF.openOK:
                            apflog("Error: Cannot open the dome", level="alert",echo=True)
                            os._exit()
                        else:
                            # lost permision during opening, happens more often than you think
                            apflog("Error: No longer have opening permission", level="error",echo=True)
                            
                    else:
                        rv = APF.evening_star()
                        if not rv:
                            apflog("evening star targeting and telescope focus did not work",level='warn', echo=True)
            
                        chk_done = "$eostele.SUNEL < %f" % (SUNEL_STARTLIM*np.pi/180.0)
                        while float(sunel.read()) > SUNEL_STARTLIM and not rising:
                            outstr = "Sun is setting and sun at elevation of %.3f" % (float(sunel.read()))
                            apflog(outstr,level='info', echo=True)
                            result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=60)
                            APF.DMReset()
                    
                elif not rising or (rising and float(sunel) < (sunel_lim - 5)):
                    APFTask.set(parent, suffix="MESSAGE", value="Open at night", wait=False)                    
                    success = opening(sunel)
                else:
                    success = True
                if success == False:
                    apflog("Error: Cannot open the dome", echo=True, level='error')
                    APF.close()
                    os._exit()


            # Check for servo errors
            if APF.slew_allowed.read(binary=True) is False and APF.isReadyForObserving()[0]:
                APFTask.set(parent, suffix="MESSAGE", value="Likely servo failure.", wait=False)                    
                if running:
                    APF.killRobot(now=True)
                
                apflog("Error: APF is open, and slew_allowed is false. Likely a servo error/amplifier fault.", level="error", echo=True)
                chk_done = "$checkapf.MOVE_PERM == true"
                result = APFTask.waitFor(self.task, True, expression=chk_done, timeout=600)
                if result:
                    rv = APF.power_down_telescope()
                    if rv:
                        apflog("APF power cycled.", echo=True)
                    else:
                        apflog("Error: APF power cycle failed.", level="error", echo=True)
                        closing(force=True)
                elif result is False and "DomeShutter" in APF.isOpen()[1]:
                    apflog("Error: After 10 min move permission did not return, and the dome is still open.", level='error', echo=True)

                    closing(force=True)
                
            # If we are open and scriptobs isn't running, start it up
            if APF.isReadyForObserving()[0] and not running and float(sunel) <= sunel_lim:
                APFTask.set(parent, suffix="MESSAGE", value="Starting scriptobs", wait=False)
                rv = checkTelState()
                if rv is False:
                    # this means that the telescope is not slewing and is not tracking
                    rv = startTelescope()
                    # if needed, will power up the Az drive and clear the estop state
                    if rv == False:
                        apflog("Telescope stopped and cannot be restarted", level='Alert', echo=True)
                        closing(force=True)

                startScriptobs()
                if not APFTask.waitFor(self.task,True,expression="$apftask.SCRIPTOBS_STATUS == 'Running'",timeout=10):
                    if failstart % 11 == 0 and failstart > 0:
                        lvl = "error"
                    else:
                        lvl = "warn"
                    apflog("scriptobs is not running just after being started!", level=lvl, echo=True)
                    APFTask.set(parent, suffix="MESSAGE",value="scriptobs is not running just after being started!",wait=False)                    

                
            # Keep an eye on the deadman timer if we are open 
            if APF.isOpen()[0] and APF.dmtime <= DMLIM:
                APFTask.set(parent, suffix="MESSAGE", value="Reseting DM timer", wait=False)                    
                APF.DMReset()
#                apflog("The APF is open, the DM timer is clicking down, and scriptobs is %s." % ( str(running)),level="debug")

            current_msg = APFTask.get("master", ["message"])
            if not APF.isOpen()[0] and not rising:
                omsg = "Waiting for sunset"
                if current_msg['message'] != omsg:
                    APFTask.set(parent, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
            if  APF.isOpen()[0] and float(sunel) > sunel_lim and not rising:
                omsg = "Waiting for sunset"
                if current_msg['message'] != omsg:
                    APFTask.set(parent, suffix="MESSAGE", value=omsg, wait=False)
                APFTask.waitFor(self.task, True, timeout=5)
            

    def stop(self):
        self.signal = False
        threading.Thread._Thread__stop(self)

