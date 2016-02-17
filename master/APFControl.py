#!/usr/bin/env  /opt/kroot/bin/kpython

# Class definition for an APF object which tracks the state of the telescope.

import ktl
import APF as APFLib
import APFTask

import subprocess
import time
import os
import os.path
import math
from datetime import datetime, timedelta

import numpy as np

from apflog import *

m1 = 22.8
windlim = 40.0
slowlim = 100
WINDSHIELD_LIMIT = 10.
wxtimeout = timedelta(seconds=1800)
SUNEL_HOR = -3.2

#ScriptDir = '@LROOT@/bin/robot/'
ScriptDir = '/usr/local/lick/bin/robot/'

deckscale = {'M': 1.0, 'W':1.0, 'N': 3.0, 'B': 0.5, 'S':2.0, 'P':1.0}


# Aquire the ktl services and associated keywords
tel        = ktl.Service('eostele')
sunelServ  = tel('SUNEL')
checkapf   = ktl.Service('checkapf')
apfmet     = ktl.Service('apfmet')
ok2open    = ktl.cache('checkapf','OPEN_OK')
dmtimer    = ktl.cache('checkapf','DMTIME')
wx         = ktl.cache('apfmet','M5WIND')
robot      = ktl.Service('apftask')
vmag       = robot['scriptobs_vmag']
ucam       = ktl.Service('apfucam')
apfteq     = ktl.Service('apfteq')
teqmode    = apfteq['MODE']
guide      = ktl.Service('apfguide')
counts     = ktl.cache('apfguide','COUNTS')
countrate  = guide['countrate']
thresh     = guide['xpose_thresh']

elapsed    = ucam['elapsed']
motor      = ktl.Service('apfmot')
decker     = motor['DECKERNAM']


def cmdexec(cmd, debug=False, cwd='./'):
    args = cmd.split()
    apflog("Executing Command: %s" % repr(cmd), echo=True)
    
    p = subprocess.Popen(args, stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=cwd)
    
    while p.poll() is None:
        l = p.stdout.readline().rstrip('\n')
        if debug: apflog(l, echo=debug)

    out, err = p.communicate()
    if debug: apflog(out, echo=debug)
    if len(err): apflog(err, echo=debug)
    ret_code = p.returncode
    if ret_code == 0:
        return True, ret_code
    else:
        return False, ret_code



def countmon(counts):
    APF.countrate = -1.0
    if counts['populated'] == False:
        return

    try:
        cnts = float(counts)
        time = float(elapsed.read(binary=True))
        APF.countrate = cnts/time
    except ZeroDivisionError:
        return
    except:
        apflog("Cannot read apfguide.counts or apfucam.elapsed")
        return


# Callback for ok2open permission
# -- Check that if we fall down a logic hole we don't error out
def okmon(ok2open):
    if ok2open['populated'] == False:
        return
    try:
        ok = ok2open # historical
    except Exception, e:
        apflog("Exception in okmon: %s" % (e), level='warn')
        return
    try:
        if not checkapf['MOVE_PERM'].read(binary=False):
            ok = False
    except Exception, e:
        apflog("Exception in okmon: %s" % (e), level='warn')
        return
    try:
        if not checkapf['USERKIND'].read(binary=True) == 3:
            ok = False
    except Exception, e:
        apflog("Exception in okmon: %s" % (e), level='warn')
        return
    APF.openOK = ok
    return

# Callback for the windspeed
def windmon(wx):
    if wx['populated'] == False:
        return
    windshield = robot["scriptobs_windshield"].read()
    try:
        wvel = float(wx)
    except Exception, e:
        apflog("Exception in windmon: %s" % (e), level='warn')
        return
        
    if APF.wslist == []:
        APF.wslist = [wvel]*20

    else:
        APF.wslist.append(wvel)
        APF.wslist = APF.wslist[-20:]

    APF.wvel = np.median(APF.wslist)


# Callback for Deadman timer
def dmtimemon(dmtime):
    if dmtime['populated'] == False:
        return
    try:
        APF.dmtime = dmtime
    except Exception, e:
        apflog("Exception in dmtimemon: %s" % (e), level='warn')


class APF:
    """ Class which creates a monitored state object to track the condition of the APF telescope. """

    # Initial seeing conditions
    seeinglist = []
    speedlist  = []
    conditions = 'bad'
    cwd        = os.getcwd()
    slowdown   = 0.0 

    # Initial Wind conditions
    wslist = []
    wdlist = []

    # KTL Services and Keywords
    tel        = ktl.Service('eostele')
    sunel      = tel('SUNEL')
    ael        = tel('AEL')
    aaz        = tel('AAZ')
    aafocus    = tel('AAFOCUS')
    dome       = ktl.Service('eosdome')
    rspos      = dome('RSCURPOS')
    fspos      = dome('FSCURPOS')
    checkapf   = ktl.Service('checkapf')
    apfmet     = ktl.Service('apfmet')
    ok2open    = checkapf('OPEN_OK')
    dmtimer    = checkapf('DMTIME')
    wx         = apfmet('M5WIND')
    mv_perm    = checkapf('MOVE_PERM')
    chk_close  = checkapf('CHK_CLOSE')
    robot      = ktl.Service('apftask')
    vmag       = robot['scriptobs_vmag']
    ldone      = robot['scriptobs_lines_done']
    sop        = robot['scriptobs_phase']
    autofoc    = robot["SCRIPTOBS_AUTOFOC"]
    slew_allowed = robot['slew_allowed']
    ucam       = ktl.Service('apfucam')
    user       = ucam['OUTFILE']
    elapsed    = ucam['elapsed']
    apfteq     = ktl.Service('apfteq')
    teqmode    = apfteq['MODE']
    guide      = ktl.Service('apfguide')
    counts     = guide['COUNTS']
    countrate  = guide['countrate']
    thresh     = guide['xpose_thresh']
    avg_fwhm   = guide['AVG_FWHM']
    motor      = ktl.Service('apfmot')
    decker     = motor['DECKERNAM']
    dewarfoc   = motor["DEWARFOCRAW"]

    def __init__(self, task="example", test=False):
        """ Initilize the current state of APF. Setup the callbacks and monitors necessary for automated telescope operation."""
        # Set up the calling task that set up the monitor and if this is a test instance
        self.test = test
        self.task = task
        
        self.cloudObsNum = 1
  
        # Set the callbacks and monitors
        self.wx.monitor()
        self.wx.callback(windmon)

        self.ok2open.monitor()
        self.ok2open.callback(okmon)

        self.dmtimer.monitor()
        self.dmtimer.callback(dmtimemon)

        self.counts.monitor()
        self.counts.callback(countmon)

        self.teqmode.monitor()
        self.vmag.monitor()
        self.autofoc.monitor()
        self.ldone.monitor()
        self.counts.monitor()
        self.decker.monitor()
        self.avg_fwhm.monitor()
        self.dewarfoc.monitor()
        self.mv_perm.monitor()
        self.chk_close.monitor()
        self.slew_allowed.monitor()

        self.sunel.monitor()
        self.aaz.monitor()
        self.ael.monitor()
        self.fspos.monitor()
        self.rspos.monitor()
        self.aafocus.monitor()

        # Grab some initial values for the state of the telescope
        
        self.wx.poll()
        self.counts.poll()
        self.ok2open.poll()

    def __str__(self):
        # Determine if the sun rising / setting check is working
        now = datetime.now()
        if now.strftime("%p") == 'AM':
            rising = True
        else:
            rising = False
        s = ''
        s += "At %s state of telescope is:\n" % str(now)
        s += "Sun elevation = %4.2f %s\n" % (self.sunel, "Rising" if rising else "Setting")
        s += "Telescope -- AZ=%4.2f  EL=%4.2f \n" % (self.aaz, self.ael)
        s += "Front/Rear Shutter=%4.2f / %4.2f\n"%(self.fspos, self.rspos)
        s += "Wind = %3.1f mph \n" % (self.wvel)
        s += "Seeing %4.2f arcsec\n" % self.seeing
        s += "Slowdown = %5.2f x\n" % self.slowdown
        #s += "Conditions are - %s\n" % self.conditions
        s += "Teq Mode - %s\n" % self.teqmode
        s += "M2 Focus Value = % 4.3f\n" % self.aafocus
        s += "Okay to open = %s -- %s\n" % (repr(self.openOK), self.checkapf['OPREASON'].read() )
        s += "Current Weather = %s\n" % self.checkapf['WEATHER'].read()
        isopen, what = self.isOpen()
        if isopen:
            s += "Currently open: %s\n" % what
        else:
            s += "Not currently open\n"
        ripd, rr = self.findRobot()
        if rr:
            s += "Robot is running\n"
        else:
            s += "Robot is not running\n"

        return s


    # Fucntion for checking what is currently open on the telescope
    def isOpen(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        what = self.checkapf("WHATSOPN").read()
        if "DomeShutter" in what or "MirrorCover" in what or "Vents" in what:
            return True, what
        else:
            return False, ''

    # Fucntion for checking what is currently open on the telescope
    def isReadyForObserving(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        what = self.checkapf("WHATSOPN").read()
        try:
            if "DomeShutter" in what and "MirrorCover" in what:
                return True, what
            else:
                return False, ''
        except Exception, e:
            apflog("Exception in isReadyForObserving: %s (what = %s)" % (e,what), level='error')
            return False, ''

    def setObserverInfo(self, num=10000, name='Robot'):
        if self.test: return
        apflog("Setting science camera parameters.")
        self.ucam('OBSERVER').write(name)
        self.ucam('OBSNUM').write(str(num))
        self.ucam('OUTDIR').write('/data/apf/')
        self.ucam('OUTFILE').write(name)

        apflog("Upadted science camera parameters:")
        apflog("Observer = %s" % self.ucam('OBSERVER').read(),echo=True)
        apflog("Output directory = %s" % self.ucam('OUTDIR').read(),echo=True)
        apflog("Observation number = %s" % self.ucam('OBSNUM').read(), echo=True)
        apflog("File prefix = %s" % self.ucam('OUTFILE').read(), echo=True)

        

    def calibrate(self, script, time):
        s_calibrate = os.path.join(ScriptDir,"calibrate")
        if self.test: 
            print "Test Mode: calibrate %s %s." % (script, time)
            APFTask.waitFor(self.task, True, timeout=10)
            return True
        if time == 'pre' or 'post':
            try:
                APFLib.write("apfmot.DEWARFOCRAW",ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True))
            except:
                apflog("Cannot read the last best fitting focus value or write the dewar focus value", level='error')
            apflog("Running calibrate %s %s" % (script, time), level = 'info')
            cmd = '%s %s %s' % (s_calibrate,script, time)
            result, code = cmdexec(cmd)
            if not result:
                apflog("%s %s failed with return code %d" % (s_calibrate, script, code),echo=True)
            expression="($apftask.CALIBRATE_STATUS != 0) and ($apftask.CALIBRATE_STATUS != 1) "
            if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                apflog("%s %s failed to exit" % (s_calibrate,script),echo=True)
                
            return result
        else:
            print "Couldn't understand argument %s, nothing was done." % time

    def focus(self, user='ucsc'):
        """Runs the focus routine appropriate for the style string."""
        if user == 'ucsc':
            if self.test: 
                APFTask.waitFor(self.task, True, timeout=10)
                print "Test Mode: Would be running focusinstr."
                return True
            else:
                apflog("Running focusinstr routine.",echo=True)
#                cmd = '/u/user/devel_scripts/ucscapf/auto_focuscube.sh pre t'
                cmdpath = '/usr/local/lick/bin/robot/'
                cmd = os.path.join(cmdpath,'focusinstr -b')
                result, code = cmdexec(cmd,cwd=os.path.curdir)
                if not result:
                    apflog("focusinstr failed with code %d" % code, echo=True)
                expression="($apftask.FOCUSINSTR_STATUS != 0) and ($apftask.FOCUSINSTR_STATUS != 1) "
                if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                    apflog("focusinstr failed to exit" ,echo=True)
                return result
        else:
            print "Don't recognize user %s. Nothing was done." % style

    def setTeqMode(self, mode):
        apflog("Setting TEQMode to %s" % mode)
        if self.test: 
            print "Would be setting TEQMode to %s" % mode
            return
        self.teqmode.write(mode)
        result = self.teqmode.waitfor('== %s' % mode, timeout=60)
        if not result:
            apflog("Error setting the TEQMODE.")
            raise RuntimeError, "Couldn't set TEQ mode"


    def clearestop(self):
        if self.test: return True
        if self.mv_perm.binary == False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        cmd = '/usr/local/lick/bin/robot/clear_estop'
        result, code = cmdexec(cmd)
        if result:
            return True
        else:
            return False


    def openat(self, sunset=False):
        """Function to ready the APF for observing. Calls either openatsunset or openatnight.
           This function will attempt to open successfully twice. If both attempts
           fail, then it will return false, allowing the master to register the error
           and behave accodingly. Otherwise it will return True. """
        # If this is a test run, just return True
        if self.test: return True

        if not self.ok2open:
            # This should really never happen. In case of a temporary condition, we give
            # a short waitfor rather than immediatly exiting.
            chk_open = "$checkapf.OPEN_OK == true"
            result = APFLib.waitFor(self.task, False, chk_open, timeout=30) 
            if not result:
                apflog("Tried calling openat with OPEN_OK = False. Can't open.", echo=True)
                apflog(self.checkapf["OPREASON"].read(), echo=True)
                return False

        if float(self.sunel) > SUNEL_HOR:
            apflog("Sun is still up. Current sunel = %4.2f. Can't open." % self.sunel, echo=True)
            return False
        
        if self.mv_perm.binary == False:
            apflog("Waiting for permission to move...", echo=True)
            chk_move = "$checkapf.MOVE_PERM == true"
            result = APFTask.waitFor(self.task, False, chk_move, timeout=600)
            if not result:
                apflog("Can't open. No move permission.",echo=True)
                return False

        # Everything seems acceptable, so lets try opening
        if sunset:
            cmd = '/usr/local/lick/bin/robot/openatsunset'
        else:
            cmd = '/usr/local/lick/bin/robot/openatnight'

        # Make two tries at opening. If they both fail return False so the caller can act
        # accordingly.
        result, code = cmdexec(cmd)
        if not result:
            apflog("First openup attempt has failed. Exit code = %d. After a pause, will make one more attempt." % code,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            result, code = cmdexec(cmd)
            if result:
                try:
                    APFLib.write("eostele.FOCUS",ktl.read("apftask","FOCUSTEL_LASTFOCUS",binary=True))
                except:
                    apflog("Cannot move secondary focus.",level="error")
                return True
            else:
                apflog("Second openup attempt also failed. Exit code %d. Giving up." % code,echo=True)
                return False
        else:
            try:
                APFLib.write("eostele.FOCUS",ktl.read("apftask","FOCUSTEL_LASTFOCUS",binary=True))
            except:
                apflog("Cannot move secondary focus.",level="error")

            return True

    def power_down_telescope(self):
        """Checks that we have the proper permission and dome is closed, then resets telescope power."""
        if self.test: return True
        cmd = "/usr/local/lick/bin/robot/power_down_telescope"
        self.DMReset()
        if self.mv_perm.binary == False:
                apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
        if not result:
            apflog("Didn't have move permission after 5 minutes.", echo=True) 
            return False
        apflog("Running power_down_telescope script")
        result, code = cmdexec(cmd)
        if not result:
            apflog("power_down_telescope has failed. Human intervention likely required.", level='error', echo=True)
            APFTask.waitFor(self.task, True, timeout=30)
        else:
            pass
        if result:    
            return True
        else:
            return False

            
    def close(self, force=False):
        """Checks that we have the proper permission, then runs the closeup script."""
        if self.test: return True
        cmd = "/usr/local/lick/bin/robot/closeup"
        if force:
            apflog("Calling a single instance of closeup. Will return regardless of result.", echo=True)
            result, code = cmdexec(cmd)
            return result
        if self.mv_perm.binary == False:
            if self.chk_close.binary == True:
                apflog("Waiting for checkapf to close up")
            else:
                apflog("Waiting for permission to move")
        chk_mv = '$checkapf.MOVE_PERM == true'
        result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
        if not result:
            apflog("Didn't have move permission after 5 minutes. Going ahead with closeup.", echo=True) 
        apflog("Running closeup script")
        attempts = 0
        close_start = datetime.now()
        while (datetime.now() - close_start).seconds < 1800:
            attempts += 1
            result, code = cmdexec(cmd)
            if not result:
                apflog("Closeup failed with exit code %d" % code, echo=True)
                if attempts == 3:
                    apflog("Closeup has failed 3 times consecutively. Human intervention likely required.", level='error', echo=True)
                APFTask.waitFor(self.task, True, timeout=30)
            else:
                break
        if result:    
            return True
        else:
            apflog("After 30 minutes of trying, closeup could not successfully complete.")

    def focusTel(self):
        """Slew the telescope to a bright star, open the shutters, and call measure_focus."""
        # Short plan
        # Check if we are tracking a star. ( apftask slew_phase? )
        # If not, getNext(Bstar) -> self.target()
        # .call(["focus_telescope", "-o"])

        if self.robot["SLEW_PHASE"].read() != "Tracking":
            apflog("We don't seem to be tracking a star. This is required to focus the telescope.", echo=True)
            apflog("focusTel() will return without doing anything.", echo=True)
            return

        apflog("Attempting to run a full focus sweep.", echo=True)
        args = ["./focus_telescope", "-o"]

        PIPE = subprocess.PIPE

        ret = subprocess.call(args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        if ret != 0:
            apflog("Focus sweep has failed.", echo=True)
        else:
            apflog("Focus finished successfully")



    def updateLastObs(self):
        """ If the last observation was a success, this function updates the file storing the last observation number and the hit_list which is required by the dynamic scheduler."""
        if self.ucam('OUTFILE').read() != 'ucsc':
            apflog("UCAM Observer name is not ucsc, so the lastObs file will not be updated.")
            apflog("Number of last observation is %s" % self.ucam('OBSNUM').read())
            return
        APFLib.write(self.robot["MASTER_LAST_OBS_UCSC"], self.ucam["OBSNUM"].read())
        with open(os.path.join(self.cwd,'lastObs.txt'),'w') as f:
                f.write("%s\n" % self.ucam('OBSNUM').read())
                apflog("Recording last ObsNum as %d" % int(self.ucam["OBSNUM"].read()))


    def updateWindshield(self, state):
        """Checks the current windshielding mode, and depending on the input and wind speed measurements makes sure it is set properly."""
        currState = self.robot["SCRIPTOBS_WINDSHIELD"].read().strip().lower()
        if state == 'on':
            if currState != 'enable':
                apflog("Setting scriptobs_windshield to Enable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Enable")
        elif state == 'off':
            if currState != 'disable':
                apflog("Setting scriptobs_windshield to Disable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Disable")
        else:
            # State must be auto, so check wind
            if currState == 'enable' and self.wvel <= WINDSHIELD_LIMIT:
                apflog("Setting scriptobs_windshield to Disable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Disable")
            if currState == 'disable' and self.wvel > WINDSHIELD_LIMIT:
                apflog("Setting scriptobs_windshield to Enable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Enable")

    def target(self, name, ra, dec, pmra, pmdec):
        """Aim the APF at the desired target. This calles prep-obs, slewlock, and focus-telescope. A workaround to relying on scriptobs."""
        if self.isOpen()[0] == False:
            apflog("APF is not open. Can't target a star while closed.",echo=True)
            return
        apflog("Targeting telescope on %s" % name, echo=True)
        # Move the shutters to fully open ( This might have already been done by open at sunset. In this case we rely on shutters to realize this and simply return)
        apflog("Fully opening shutters.")
        result = subprocess.call(["shutters","-o"])
        if result != 0:
            apflog("Shutters returned error code %d. Targeting object %s has failed." % (result, name),level='error',echo=True)
            return
        # Call prep-obs
        apflog("Calling prep-obs.")
        result = subprocess.call(['prep-obs'])
        if result != 0:
            apflog("Prep-obs returned error code %d. Targeting object %s has failed." % (result, name),level='error',echo=True)
            return
        # Slew to the specified RA and DEC, set guide camera settings, and centerup( Slewlock )
        apflog("Calling slewlock.")
        result = subprocess.call(["slewlock","target",name,str(ra),str(dec),str(pmra),str(pmdec),'210'])
        if result != 0:
            apflog("Slewlock returned error code %d. Targeting object %s has failed." % (result, name), level='error',echo=True)
            return
        # Focus the telescope?
            

    def checkClouds(self, target):
        """
        This function will take a test exposure of a B-Star. By using scriptobs to take this exposure,
        the avg_fwhm and countrate keywords will be updated.
        The resulting file will be stored as Heimdallr1.fits?
        :return:
        """

        apflog("checkClouds(): starting transparency check", echo=True)
        
        # Things to reset after this
        apflog("checkClouds(): Storing current UCAM keywords for later...", echo=True)
        # UCAM outfile
        obs_file = self.ucam["OUTFILE"].read()
        # UCAM Obsnumber
        obs_num = self.ucam["OBSNUM"].read()
        robot['MASTER_VAR_2'].write(obs_num)
        # scriptobs_lines_done
        lines_done = int(self.robot["SCRIPTOBS_LINES_DONE"])
        APFLib.write(self.autofoc, "robot_autofocus_enable")

        apflog("checkClouds(): File name=%s - Number=%d - Lines Done=%d" % (obs_file, int(obs_num), lines_done) )

        APFLib.write(self.ucam["OUTFILE"], "heimdallr")
        APFLib.write(self.ucam["OBSNUM"], self.cloudObsNum)
        APFLib.write(self.ucam["RECORD"], "No")


        cmdexec('/usr/local/lick/bin/robot/prep-obs')
        args = ['/usr/local/lick/bin/robot/scriptobs', '-dir', os.path.join(self.cwd,'checkClouds/')]

        apflog("checkClouds(): Observing %s as a test star to check the clouds/seeing/transparency" % target["NAME"], echo=True)
        outfile = open('robot.log', 'a')
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=outfile, stderr=outfile)
        foo, bar = p.communicate(target["SCRIPTOBS"])
        outfile.close()
        ret = p.returncode
        if ret != 0:
            apflog("checkClouds(): The cloud cover test exposure has failed.", echo=True)
        else:
            apflog("checkClouds(): Cloud check was successful", echo=True)
            self.cloudObsNum += 2

        apflog("checkClouds(): Resetting UCAM keywords to previous values.", echo=True)
        APFLib.write(self.ucam["RECORD"], "Yes")
        APFLib.write(self.ucam["OUTFILE"], obs_file)
        APFLib.write(self.ucam["OBSNUM"], obs_num)
        APFLib.write(self.robot["SCRIPTOBS_LINES_DONE"], lines_done)
	APFTask.waitFor(self.task, True, timeout=5)
        apflog("checkClouds(): Keywords successfully written.", echo=True)
	apflog("checkClouds(): New values")
	apflog("checkClouds(): File=%s - Num=%d - LinesDone=%d" % (self.ucam["OUTFILE"].read(), int(self.ucam["OBSNUM"].read()), int(self.robot["SCRIPTOBS_LINES_DONE"].read()) ) )

    def observe(self, observation, skip=0):
        """ Currently: Takes a string which is the filename of a properly formatted star list. """

        if self.test:
            apflog("Would be taking observation in starlist %s" % observation)
            APFTask.waitFor(self.task, True, timeout=300)
            return
        # Make sure the telescope autofocus is enabled 
        APFLib.write(self.autofoc, "robot_autofocus_enable")
        APFLib.write(self.ucam["RECORD"], "Yes")
        chk_foc = '$apftask.SCRIPTOBS_AUTOFOC == robot_autofocus_enable'
        result = APFTask.waitFor(self.task, False, chk_foc, timeout=60)
        if not result:
            apflog("Error setting scriptobs_autofoc", echo=True)
            return
        # Make sure APFTEQ is in night mode for observations
        if self.teqmode.read() != 'Night':
            self.setTeqMode('Night')
        # Check the instrument focus for a reasonable value
        if self.dewarfoc > 8600 or self.dewarfoc < 8400:
            lastfit_dewarfoc = ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True)
            apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values. Resetting to last derived value %d" % (self.dewarfoc,lastfit_dewarfoc), level = "error", echo=True)
            APFLib.write("apfmot.DEWARFOCRAW",lastfit_dewarfoc)
            
        # Check Telescope M2 Focus
        
        robotdir = "/usr/local/lick/bin/robot/"
        infile = open(observation,'r')
        outfile = open("robot.log", 'a')
        
        if skip != 0:
            args = ['/usr/local/lick/bin/robot/scriptobs', '-dir', os.getcwd(),'-skip', str(skip)]
        else:
            args = ['/usr/local/lick/bin/robot/scriptobs', '-dir', os.getcwd()]

        p = subprocess.Popen(args,stdin=infile, stdout=outfile,stderr = subprocess.PIPE, cwd=robotdir)

    def DMReset(self):
        try:
            APFLib.write(self.checkapf['ROBOSTATE'], "master operating")
        except Exception, e:
            try:
                ukind = self.checkapf['USERKIND'].read()
            except:
                ukind = "Unknown"
            ostr = "Error: Cannot write to ROBOSTATE, USERKIND = %s, reason: %s" % (ukind,e)
            apflog(ostr,level='error',echo=True)

    def DMZero(self):
        try:
            if self.checkapf['DMTIME'].read(binary=True) < 1:
                APFLib.write(self.checkapf['DMTIME'], -1)
        except Exception, e:
            ostr = "Warning: cannot touch DM Timer: %s " %( e)
            apflog(ostr,level='warn',echo=True)

    def findRobot(self):
        """Trys to find a running instance of robot.csh. Returns the PID along with a boolean representing if the robot was succesfully found."""
        rpid = self.robot['SCRIPTOBS_PID'].read(binary=True)
        if rpid == '' or rpid == -1:
            return rpid, False
        else:
            return rpid, True

    def startRobot(self):
        """Start an instance of scriptobs. Returns the result from subprocess.Popen()."""
        # For running in test mode
        if self.test:
            apflog("Would be taking observation in starlist %s" % observation)
            APFTask.waitFor(self.task, True, timeout=300)
            return
        
        # Make sure the telescope autofocus is enabled 
        APFLib.write(self.autofoc, "robot_autofocus_enable")
        chk_foc = '$apftask.SCRIPTOBS_AUTOFOC == robot_autofocus_enable'
        result = APFTask.waitFor(self.task, False, chk_foc, timeout=60)
        if not result:
            apflog("Error setting scriptobs_autofoc", echo=True)
            return
        # Make sure APFTEQ is in night mode for observations
        if self.teqmode.read() != 'Night':
            self.setTeqMode('Night')
        # Check the instrument focus for a reasonable value
        if self.dewarfoc > 8600 or self.dewarfoc < 8400:
            lastfit_dewarfoc = ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True)
            apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values. Resetting to last derived value %d" % (self.dewarfoc,lastfit_dewarfoc), level = "error", echo=True)
            APFLib.write("apfmot.DEWARFOCRAW",lastfit_dewarfoc)

        robotdir = "/usr/local/lick/bin/robot/"

        rv, retc = cmdexec(os.path.join(robotdir,"slew --hold"))
        if not rv:
            return rv
        rv, retc = cmdexec(os.path.join(robotdir,"prep-obs"))
        if not rv:
            return rv
        # Start scriptobs
        outfile = open("robot.log", 'a')
        args = ['/usr/local/lick/bin/robot/scriptobs', '-dir', os.getcwd()]

        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=outfile, stderr=outfile)
        
        return p
        

    def killRobot(self, now=False):
        """ In case during an exposure there is a need to stop the robot and close up."""
        apflog("Terminating Robot.csh")
        if now:
            apflog("Abort exposure, terminating robot now.")
        else:
            if not ucam['EVENT_STR'].read() == "ControllerReady":
                apflog("Waiting for current exposure to finish.")
                ucam['EVENT_STR'].waitfor(" = ReadoutBegin", timeout=1200)
        apflog("Killing Robot.")
        ripd, running = self.findRobot()
        if running:
            try:
                APFLib.write(self.robot['scriptobs_control'], "abort")
            except Exception, e:
                errstr = "Cannot abort scriptobs: %s" % (e)
                apflog(errstr,level="Warn",echo=True)


if __name__ == '__main__':
    print "Testing telescope monitors, grabbing and printing out current state."

    task = 'example'

    APFTask.establish(task, os.getpid())
    apf = APF(task=task,test=False)

    # Give the monitors some time to start up
    APFTask.waitFor(task, True,timeout=10)

    
    print str(apf)

    while True:
        try:
            if raw_input("Print Telescope State? (y/n): ") != 'y':
                break
        except KeyboardInterrupt:
            break
        else:
            print str(apf)


        






