# Class definition for an APF object which tracks the state of the telescope.

import subprocess
import time
import os
import os.path
import math
from datetime import datetime, timedelta

import numpy as np

try:
    from apflog import *
    import ktl
    import APF as APFLib
    import APFTask
except:
    from fake_apflog import *


windlim = 40.0
slowlim = 100
WINDSHIELD_LIMIT = 10. # mph at the APF
TEMP_LIMIT = 35. # deg F at the APF
wxtimeout = timedelta(seconds=1800)
SUNEL_HOR = -3.2
DEWARMAX = 8600
DEWARMIN = 8400

#ScriptDir = '@LROOT@/bin/robot/'
ScriptDir = '/usr/local/lick/bin/robot/'

# Aquire the ktl services and associated keywords
tel        = ktl.Service('eostele')
sunelServ  = tel('SUNEL')
apfmet     = ktl.Service('met3apf')
checkapf   = ktl.Service('checkapf')
ok2open    = ktl.cache('checkapf','OPEN_OK')
dmtimer    = ktl.cache('checkapf','DMTIME')
wx         = ktl.cache('met3apf','M5WIND')

robot      = ktl.Service('apftask')
vmag       = robot['SCRIPTOBS_VMAG']

ucam       = ktl.Service('apfucam')
apfteq     = ktl.Service('apfteq')
teqmode    = apfteq['MODE']
guide      = ktl.Service('apfguide')
counts     = ktl.cache('apfguide','COUNTS')
kcountrate     = ktl.cache('apfguide','COUNTRATE')
thresh     = guide['XPOSE_THRESH']
elapsed    = ktl.cache('apfucam','ELAPSED')
motor      = ktl.Service('apfmot')
decker     = motor['DECKERNAM']


def cmdexec(cmd, debug=False, cwd='./'):
    args = ["apftask","do"]
    args = args + cmd.split()
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
    if counts['populated'] == False:
        return
    try:
        cnts = float(counts.read(binary=True))
    except:
        return
    try:
        time = float(elapsed.read(binary=True))
    except:
        return

    try:
        APF.ccountrate = cnts/time
    except:
        return


def countratemon(kcountrate):
     if kcountrate['populated'] == False:
         return

     try:
         ctr = float(kcountrate['binary'])
     except:
         apflog("Cannot read apfguide.COUNTRATE",level='warn',echo=True)
         return
     APF.countrate *=  (1.0*APF.ncountrate)/(APF.ncountrate+1)
     APF.countrate += ctr/(APF.ncountrate+1)
     APF.ncountrate += 1
     return

def eventmon(event):
    if event['populated'] == False:
        return

    try:
        eventval = event.read(binary=True)
        if eventval == 0 or eventval == 7 :
            APF.ncountrate = 0
    except:
        return


    try:
        cnts = float(counts.read(binary=True))
        time = float(elapsed.read(binary=True))
    except:
        return
    try:
        APF.ccountrate = cnts/time
    except:
        return

# Callback for ok2open permission
# -- Check that if we fall down a logic hole we don't error out
def okmon(ok2open):
    if ok2open['populated'] == False:
        return
    try:
        ok = ok2open # historical
    except Exception, e:
        apflog("Exception in okmon for checkapf.OPEN_OK: %s" % (e), level='error')
        return
    try:
        if checkapf['MOVE_PERM'].read(binary=False) == False:
            ok = False
    except Exception, e:
        apflog("Exception in okmon for checkapf.MOVE_PERM: %s" % (e), level='error')
        return
    try:
        if not checkapf['USERKIND'].read(binary=True) == 3:
            ok = False
    except Exception, e:
        apflog("Exception in okmon checkapf.USERKIND: %s" % (e), level='error')
        return
    APF.openOK = ok
    return



# Callback for the windspeed
def windmon(wx):
    if wx['populated'] == False:
        return
    try:
        wvel = float(wx)
    except Exception, e:
        apflog("Exception in windmon: %s" % (e), level='error')
        return
        
    if APF.wslist == []:
        APF.wslist = [wvel]*20

    else:
        APF.wslist.append(wvel)
        APF.wslist = APF.wslist[-20:]

    APF.wvel = np.median(APF.wslist)

def altwindmon(wx):
    if wx['populated'] == False:
        return
    try:
        downval = APF.down.read(binary=True)
    except Exception, e:
        apflog("Exception in altwindmon: %s" % (e), level='error')
        return
    if downval == 0:
        return 
    try:
        wvel = float(wx)
    except Exception, e:
        apflog("Exception in altwindmon: %s" % (e), level='error')
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
        apflog("Exception in dmtimemon: %s" % (e), level='error')


class APF:
    """ Class which creates a monitored state object to track the condition of the APF telescope. """

    # Initial seeing conditions
    seeinglist = []
    speedlist  = []
    cwd        = os.getcwd()
    slowdown   = 0.0 
    ncountrate = 0
    countrate = 0.0
    ccountrate = 0.0        

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
    shclosed   = dome('SHCLOSED')

    eostdio    = ktl.Service('eostdio')
    mcopen     = eostdio('MCOPEN')
    
    checkapf   = ktl.Service('checkapf')
    ok2open    = checkapf('OPEN_OK')
    dmtimer    = checkapf('DMTIME')
    whatsopn   = checkapf('WHATSOPN')  
    mv_perm    = checkapf('MOVE_PERM')
    chk_close  = checkapf('CHK_CLOSE')

    apfmet     = ktl.Service('met3apf')
    wx         = apfmet('M5WIND')
    down       = apfmet('M5DOWN')
    altwx      = apfmet('M3WIND')
    temp       = apfmet('M5OUTEMP')

    robot        = ktl.Service('apftask')
    vmag         = robot['SCRIPTOBS_VMAG']
    ldone        = robot['SCRIPTOBS_LINES_DONE']
    line         = robot['SCRIPTOBS_LINE']
    sop          = robot['SCRIPTOBS_PHASE']
    message      = robot['SCRIPTOBS_MESSAGE']  
    autofoc      = robot["SCRIPTOBS_AUTOFOC"]
    slew_allowed = robot['SLEW_ALLOWED']
    observed     = robot['SCRIPTOBS_OBSERVED']

    ucam       = ktl.Service('apfucam')
    user       = ucam['OUTFILE']
    elapsed    = ucam['ELAPSED']
    obsnum     = ucam['OBSNUM']
    event      = ucam['EVENT']
    combo_ps   = ucam['COMBO_PS']
    nerase     = ucam['NERASE']

    apfschedule= ktl.Service('apfschedule')
    
    apfteq     = ktl.Service('apfteq')
    teqmode    = apfteq['MODE']

    guide      = ktl.Service('apfguide')
    counts     = guide['COUNTS']
    kcountrate     = guide['COUNTRATE']
    thresh     = guide['xpose_thresh']
    avg_fwhm   = guide['AVG_FWHM']

    motor      = ktl.Service('apfmot')
    decker     = motor['DECKERNAM']
    dewarfoc   = motor["DEWARFOCRAW"]

    eosgcam    = ktl.Service('eosgcam')
    fits3pre   = eosgcam('FITS3PRE')
    save3d     = eosgcam('SAVE3D')

    apfmon     = ktl.Service('apfmon')

    def __init__(self, task="example", test=False):
        """ Initilize the current state of APF. Setup the callbacks and monitors necessary for automated telescope operation."""
        # Set up the calling task that set up the monitor and if this is a test instance
        self.test = test
        self.task = task
        
        self.cloudObsNum = 1

        self.rising = self.sunRising()
        
        # Set the callbacks and monitors
        self.wx.monitor()
        self.wx.callback(windmon)

        self.altwx.monitor()
        self.altwx.callback(altwindmon)

        self.ok2open.monitor()
        self.ok2open.callback(okmon)

        self.dmtimer.monitor()
        self.dmtimer.callback(dmtimemon)

        self.kcountrate.monitor()
        self.kcountrate.callback(countratemon)

        self.elapsed.monitor()

        self.obsnum.monitor()
        self.obsnum.callback(self.updateLastObs)

        self.event.monitor()
        self.event.callback(eventmon)

        self.nerase.monitor()

        self.down.monitor()
        self.temp.monitor()
        self.whatsopn.monitor()

        self.counts.monitor()
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
        self.temp.poll()
        self.counts.poll()
        self.ok2open.poll()

    def __str__(self):
        # Determine if the sun rising / setting check is working
        now = datetime.now()
        s = ''
        s += "At %s state of telescope is:\n" % str(now)
        s += "Sun elevation = %4.2f %s\n" % (self.sunel, "Rising" if self.rising else "Setting")
        s += "Telescope -- AZ=%4.2f  EL=%4.2f \n" % (self.aaz, self.ael)
        s += "Front/Rear Shutter=%4.2f / %4.2f\n"%(self.fspos, self.rspos)
        s += "Wind = %3.1f mph \n" % (self.wvel)
        s += "Slowdown = %5.2f x\n" % self.slowdown
        s += "countrate = %5.2g cts/s\n" % self.countrate
        s += "kcountrate = %5.2g cts/s\n" % self.kcountrate
        s += "ncountrate = %d frames \n" % self.ncountrate
        s += "elapsed = %5.2f sec \n" % self.elapsed
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


    def sunRising(self):
        now = datetime.now()
        if now.strftime("%p") == 'AM':
            return True
        else:
            return False

    def initGuidecam(self):
        self.save3d.write(False,binary=True)
        self.fits3pre.write('')
        return
    
    # Fucntion for checking what is currently open on the telescope
    def isOpen(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        try:
            whatstr = str(self.whatsopn)
            what = whatstr.split()
        except:
            apflog("checkapf.WHATSOPN returned a value that str.split cannot split",level='warn',echo=True)
            return False, ''
        if hasattr(what,'__iter__'):
            if "DomeShutter" in what or "MirrorCover" in what or "Vents" in what:
                return True, what
            else:
                return False, ''
        else:
            return False, ''

    # Fucntion for checking what is currently open on the telescope

    def isReadyForObservingDirect(self):
        what = ''
        rv = False
        try:
            ismcopen = self.mcopen.read(binary=True)
        except:
            return False, ''
        try:
            isshutterclosed = self.shclosed.read(binary=True)
        except:
            return False, ''
            
        if ismcopen:
            what = what + "MirrorCover"
            rv = True
        if isshutterclosed == False:
            if len(what) > 0 :
                what = what + " DomeShutter"
            else:
                what = "DomeShutter"
            if rv:
                rv = True
        else:
            rv = False
        return rv, what
        
    def isReadyForObserving(self):
        """Returns the state of checkapf.WHATSOPN as a tuple (bool, str)."""
        try:
            whatstr = str(self.whatsopn)
            what = whatstr.split()
        except:
            apflog("checkapf.WHATSOPN returned a value that str.split cannot split",level='warn',echo=True)
            return self.isReadyForObservingDirect()

        
        if hasattr(what,'__iter__'):
            if "DomeShutter" in what and "MirrorCover" in what:
                return True, what
            else:
                return False, ''
        else:
            return self.isReadyForObservingDirect()

    def setObserverInfo(self, num=10000, name='Robot', owner='public'):
        if self.test: return
        apflog("Setting science camera parameters.")
        self.ucam('OBSERVER').write(name)
        self.apfschedule('OWNRHINT').write(owner)        
        self.ucam('OUTFILE').write(name)
        self.ucam('OUTDIR').write('/data/apf/')
        self.ucam('OBSNUM').write(str(num))
        self.robot['UCAMLAUNCHER_UCAM_PCC'].write(0)

        apflog("Updated science camera parameters:")
        apflog("Observer = %s" % self.ucam('OBSERVER').read(),echo=True)
        apflog("Ownrhint = %s" % self.apfschedule('OWNRHINT').read(),echo=True)        
        apflog("Output directory = %s" % self.ucam('OUTDIR').read(),echo=True)
        apflog("File prefix = %s" % self.user.read(), echo=True)
        apflog("Observation number = %s" % self.obsnum.read(), echo=True)

        return

    def instrPermit(self):
        instr_perm = ktl.read("checkapf","INSTR_PERM",binary=True)
        userkind = ktl.read("checkapf","USERKIND",binary=True)
        while not instr_perm or userkind != 3:
            apflog("Waiting for instrument permission to be true and userkind to be robotic")
            APFTask.waitfor(self.task, True, expression="$checkapf.INSTR_PERM = true", timeout=600)
            APFTask.waitfor(self.task, True, expression="$checkapf.USERKIND = robotic", timeout=600)
            instr_perm = ktl.read("checkapf", "INSTR_PERM", binary=True)
            userkind = ktl.read("checkapf", "USERKIND", binary=True)

        return True

    def turnOffLamps(self):
        for lamp in ("HALOGEN2","HALOGEN1","THORIUM1","THORIUM2"):
            try:
                rv = ktl.write("apfmot",lamp,"Off",wait=False)
            except Exception, e:
                apflog("Exception: %s" % (e),echo=True,level="alert")
                rv = False
            if rv is False:
                apflog("Cannot turn off lamp %s" % (lamp),echo=True,level="alert")
        return rv
    

    def writeStages(self,stagelist,component,state):
        rv = True
        for stage in stagelist:
            curkwd = stage + component
            try:
                crv = ktl.write("apfmot",curkwd,state,wait=True,timeout=10)
                if crv is False:
                    rv = False
            except:
                rv = False

    def enableObsInst(self):

        rv = True

        stagelist = ['CALMIRROR','CALSOURCE','IODINE','GUIDEFOC']
        rv = self.writeStages(stagelist,'MOE','Off')
        rv = self.writeStages(stagelist,'MOO','Off')
        rv = self.writeStages(stagelist,'MOD','Pos')
        stagelist = ['ADC','DECKER','DEWARFOC']
        rv = self.writeStages(stagelist,'MOE','On')
        rv = self.writeStages(stagelist,'MOO','On')
        rv = self.writeStages(stagelist,'MOD','Pos')

        return rv

    def enableCalInst(self):

        rv = True
        stagelist = ['ADC','CALMIRROR','CALSOURCE','IODINE','GUIDEFOC']
        rv = self.writeStages(stagelist,'MOE','Off')
        rv = self.writeStages(stagelist,'MOO','Off')
        rv = self.writeStages(stagelist,'MOD','Pos')
        stagelist = ['DECKER','DEWARFOC']
        rv = self.writeStages(stagelist,'MOE','On')
        rv = self.writeStages(stagelist,'MOO','On')
        rv = self.writeStages(stagelist,'MOD','Pos')
        return rv

       
    def disableInst(self):

        stagelist = ['ADC','GUIDEFOC','CALMIRROR','CALSOURCE','IODINE']
        rv = self.writeStages(stagelist,'MOE','Off')
        rv = self.writeStages(stagelist,'MOO','Off')        
        rv = self.writeStages(['DECKER','DEWARFOC'],'MOE','On')
        return rv

    def turnOffInst(self):

        stagelist = ['ADC','GUIDEFOC','CALMIRROR','CALSOURCE','IODINE','DECKER','DEWARFOC']
        rv = self.writeStages(stagelist,'MOE','Off')
        rv = self.writeStages(stagelist,'MOO','Off')        
        return rv
    
    def focusinstr(self):
        self.instrPermit()
        rv = self.enableCalInst()
        if rv is False:
            try:
                ip = checkapf['INSTR_PERM'].read()
            except:
                ip = 'Unknown'
            apflog("Cannot enable instrument to move stages but instr_perm is %s" % (ip), level='alert',echo=True)
            return rv
        owner = self.apfschedule('OWNRHINT').read()        
        self.apfschedule('OWNRHINT').write('public')        
        
        lastfocus_dict = APFTask.get("focusinstr", ["lastfocus","nominal"])
        if float(lastfocus_dict["lastfocus"]) > DEWARMAX or float(lastfocus_dict["lastfocus"]) < DEWARMIN:
            lastfocus_dict["lastfocus"] =  lastfocus_dict["nominal"]
        result = self.focus()

        apfmot = ktl.Service('apfmot')
        dewarfocraw = apfmot['DEWARFOCRAW'].read(binary=True)
        
        if not result or (dewarfocraw > DEWARMAX or dewarfocraw < DEWARMIN):
            flags = "-b"
            focusdict = APFTask.get("focusinstr", ["PHASE"])
            instr_perm = ktl.read("checkapf", "INSTR_PERM", binary=True)
            if not instr_perm:
                self.instrPermit()
                if len(focusdict['PHASE']) > 0:
                    flags = " ".join(["-p", focusdict['phase']])
            else:
                apflog("Focusinstr has failed. Setting to %s and trying again." % (lastfocus_dict["lastfocus"]), level='error', echo=True)
                APFLib.write("apfmot.DEWARFOCRAW", lastfocus_dict["lastfocus"])
            result = self.focus(flags=flags)
            if not result:
                apflog("Focusinstr has failed. Setting to %s and exiting." % (lastfocus_dict["lastfocus"]), level='error', echo=True)
        self.apfschedule('OWNRHINT').write(owner)        

        return result
        
    def calibrate(self, script, time):
        s_calibrate = os.path.join(ScriptDir,"calibrate")
        if self.test: 
            print "Test Mode: calibrate %s %s." % (script, time)
            APFTask.waitFor(self.task, True, timeout=10)
            return True
        if time == 'pre' or 'post':

            rv = self.enableCalInst()
            if rv is False:
                try:
                    ip = checkapf['INSTR_PERM'].read()
                except:
                    ip = 'Unknown'
                apflog("Cannot enable instrument to move stages but instr_perm is %s" % (ip), level='alert',echo=True)
                return rv
            
            try:
                APFLib.write("apfmot.DEWARFOCRAW",ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True))
            except:
                apflog("Cannot read the last best fitting focus value or write the dewar focus value", level='error')
            if self.dewarfoc > DEWARMAX or self.dewarfoc < DEWARMIN:
                apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values." % (self.dewarfoc), level = "error", echo=True)
                return False
            apflog("Running calibrate %s %s" % (script, time), level = 'info')
            owner = self.apfschedule('OWNRHINT').read()        
            self.apfschedule('OWNRHINT').write('public')        
            
            cmd = '%s %s %s' % (s_calibrate,script, time)
            result, code = cmdexec(cmd,debug=True,cwd=os.getcwd())
            if not result:
                apflog("%s %s failed with return code %d" % (s_calibrate, script, code),echo=True)
            expression="($apftask.CALIBRATE_STATUS != 0) and ($apftask.CALIBRATE_STATUS != 1) "
            if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                apflog("%s %s failed to exit" % (s_calibrate,script),echo=True)
            self.apfschedule('OWNRHINT').write(owner)        
                
            return result
        else:
            print "Couldn't understand argument %s, nothing was done." % time

    def focus(self,flags="-b"):
        """Runs the focus routine appropriate for the user."""

        if self.test: 
            APFTask.waitFor(self.task, True, timeout=10)
            print "Test Mode: Would be running focusinstr."
            return True
        else:
            supplies = ('PS1_48V_ENA', 'PS2_48V_ENA')
            for keyword in supplies:
                value = motor[keyword].read(binary=True)
                if value != 1:
                    motor[keyword].write('Enabled', wait=False)
                    
            apflog("Running focusinstr routine.",echo=True)
            cmdpath = '/usr/local/lick/bin/robot/'
            execstr = " ".join(['focusinstr',flags])
            cmd = os.path.join(cmdpath,execstr)
            result, code = cmdexec(cmd,debug=True,cwd=os.getcwd())
            if not result:
                apflog("focusinstr failed with code %d" % code, echo=True)
                result = False
                
            expression="($apftask.FOCUSINSTR_STATUS == 3)"
            if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                apflog("focusinstr failed" ,echo=True, level="error")
                result = False
            expression="($apftask.FOCUSINSTR_LASTFOCUS > 0)"
            if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                apflog("focusinstr failed to find an adequate focus" ,echo=True, level="error")
                result = False
            return result


    def findStar(self):
        ra = self.tel['RA'].read()
        dec = self.tel['DEC'].read()
        rah,ram,ras = ra.split(":") 
        decd,decm,decs = dec.split(":")
        cmdpath = '/usr/local/lick/bin/robot/'
        cmd = os.path.join(cmdpath,"closest")
        cmdargs =  [cmd, rah,ram, ras, decd,decm,decs, "5","1","8"]
        sfncat = "/usr/local/lick/data/apf/StarCatalog.dat"
        try:
            starcat = open(sfncat)
        except:
            apflog("Cannot open file %s" % (sfncat), level="warn",echo=True)
            return False
        #$line[3] $line[4] $line[5] $line[6] $line[7] 5 1 "8"] < /usr/local/lick/data/apf/StarCatalog.dat"               
        p = subprocess.Popen(cmdargs, stdin=starcat, stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=os.path.curdir)
        out, err = p.communicate()
        ret_code = p.returncode
        if ret_code != 0:
            apflog(out,echo=True)            
            apflog(err, level="warn",echo=True)
            return False
        else:
            return out.split()
    
    def slew(self,star):
        cmdpath = '/usr/local/lick/bin/robot/'
        cmd = os.path.join(cmdpath,'slewlock')
        try:
            ra = float(star[1])
            ra *= 3.819718
            dec = float(star[2])
            dec *= 57.295779
        except:
            return False
        cmd +=  ' %s %s %f %f %s %s %d ' % ("reference",star[0],ra, dec,star[4],star[5],210)
        if self.test:
            apflog("Would slew by executing %s" %(cmd), echo=True)
        else:
            apflog("Slewing by executing %s" %(cmd), echo=True)
            result, code = cmdexec(cmd,cwd=os.path.curdir,debug=True)
            if not result:
                apflog("Failed at slewing: %s" %(code), level="error", echo=True)

        return result


    def runFocustel(self):
        """Runs the telescope focus routine."""
        el = self.tel['EL'].read(binary=True)
        cfspos = self.fspos.read(binary=True)
        crspos = self.rspos.read(binary=True)

        if abs(el - cfspos) < 2.5 or abs(el - crspos) < 2.5:
            apflog("Cannot focus, telescope too close to shutter", level="warn", echo=True)
            return False
           
        if self.test: 
            APFTask.waitFor(self.task, True, timeout=10)
            apflog("Test Mode: Would be running focus_telescope.",echo=True)
            return True
        else:
            apflog("Running focus_telescope routine.",echo=True)
            cmdpath = '/usr/local/lick/bin/robot/'
            cmd = os.path.join(cmdpath,'focus_telescope')
            result, code = cmdexec(cmd,cwd=os.path.curdir)
            try:
                self.guide['MODE'].write('Guide')
            except:
                apflog('Cannot modify apfguide.MODE to Guide.',level='error',echo=True)

            if not result:
                apflog("focustel failed with code %d" % code, echo=True)
                expression="($apftask.FOCUSINSTR_STATUS != 0) and ($apftask.FOCUSINSTR_STATUS != 1) "
                if not APFTask.waitFor(self.task,True,expression=expression,timeout=30):
                    apflog("focus_telescope failed to exit" ,echo=True)
                return result
            return True

    def runAutoexposure(self,ind=5):
        cmdpath = '/usr/local/lick/bin/robot/'
        cmd = os.path.join(cmdpath,'autoexposure')
        istr = "%d" % (ind)
        cmdargs = cmd
#       cmdargs = cmd + " -i " + istr
        result, code = cmdexec(cmdargs,cwd=os.path.curdir)

#        result, code = cmdexec([cmd,istr],cwd=os.path.curdir)
        if not result:
            apflog("autoexposure failed with code %d" % code, echo=True)
        return result

    def runCenterup(self):
        cmdpath = '/usr/local/lick/bin/robot/'
        cmd = os.path.join(cmdpath,'centerup')
        result, code = cmdexec(cmd,cwd=os.path.curdir)
        if not result:
            apflog("centerup failed with code %d" % code, echo=True)
        return result

    def focusTel(self):
        star = self.findStar()
        if not star:
            apflog("Cannot find star near current position!?",level='error',echo=True)
            return False
        apflog("Targeting telescope on %s" % star[0], echo=True)
        try:
            self.vmag.write(star[6])
        except Exception, e:
            apflog("Cannot write SCRIPTOBS_VMAG: %s" % (e), level='error',echo=True)
        try:
            sline = "%s %s %s pmra=%s pmdec=%s vmag=%s # end" % (star[0],star[1],star[2],star[4],star[5],star[6])
            self.line.write(sline)
        except Exception, e:
            apflog("Cannot write SCRIPTOBS_LINE: %s" % (e), level='error',echo=True)
        if self.slew(star):
            return self.runFocustel()
        return False
    
                
    def setTeqMode(self, mode):
        apflog("Setting TEQMode to %s" % mode)
        if self.test: 
            print "Would be setting TEQMode to %s" % mode
            return
        self.teqmode.write(mode,wait=False)
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
        result, code = cmdexec(cmd,debug=True,cwd=os.getcwd())
        if result:
            return True
        else:
            return False

    def statesSet(self):
        # there are three states - but we do not care about ESTOPST, that is will be cleared in openatsunset/openatnight
        if self.dome['ECLOSEST']:
            return True
        if self.dome['ESECURST']:
            return True
        return False


    def homeTelescope(self):
        rv, rc = cmdexec("/usr/local/lick/bin/robot/slew --home")
        try:
            homed = self.apfmon('ELHOMERIGHTSTA').read(binary=True)
            if rc == 0 and homed == 2:
                return True
            else:
                apflog("cannot home telescope" % (e),level='Alert',echo=True)
                return False
        except:
            apflog("cannot home telescope and/or cannot read apfmon keyword" % (e),level='Alert',echo=True)
            return False

    def checkHome(self,home=True):
        try:
            homed = self.apfmon('ELHOMERIGHTSTA').read(binary=True)
        except Exception, e:
            apflog("apfmon.ELHOMERIGHTSTA cannot be read: %s" % (e),level='Alert',echo=True)
            return False
        if homed == 2:
            return True
        else:
            if homed == 5 or homed == 6:
                if home:
                    self.homeTelescope()
                else:
                    apflog("Telescope needs to be homed",level='Alert',echo=True)
                    return False
            else:
                apflog("apfmon.ELHOMERIGHTSTA value is %d" % (homed),level='Alert',echo=True)
                return False
                
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

        if self.statesSet():
            apflog("An unusal emergency state is set.", level="error",echo=True)
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
            if not result:
                apflog("Second openup attempt also failed. Exit code %d. Giving up." % code,echo=True)
                return False
        rv = self.checkHome()
        if rv == False:
            return False
        try:
            APFLib.write("eostele.FOCUS",ktl.read("apftask","FOCUSTEL_LASTFOCUS",binary=True))
        except:
            apflog("Cannot move secondary focus.",level="error")
            return False
        return True

    def powerDownTelescope(self):
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
        # one last check
        
        apflog("Running power_down_telescope script")
        result, code = cmdexec(cmd)
        if not result:
            apflog("power_down_telescope has failed. Human intervention likely required.", level='error', echo=True)
        else:
            pass
        if result:    
            return True
        else:
            return False


    def servoFailure(self):
        """checks for amplifier faults"""
        servo_failed = False
        prefixs = ["AZ","EL","FA","FB","FC","TR" ]
        for pr in prefixs:
            nm = pr + "AMPFLT"
            val = tel[nm].read(binary=True)
            if val:
                servo_failed = True
                
        if servo_failed:
            return self.powerDownTelescope()
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
            result = APFTask.waitFor(self.task, False, chk_mv, timeout=300)
            if not result:
                apflog("Didn't have move permission after 5 minutes. ", echo=True) 
                break
            attempts += 1
            result, code = cmdexec(cmd)
            if not result:
                apflog("Closeup failed with exit code %d" % code, echo=True)
                if attempts == 2:
                    if self.servoFailure():
                        apflog("Servo amplifier failure, power cycled telescope",echo=True)
                if attempts == 3:
                    lstr = "Closeup has failed 3 times consecutively. Human intervention likely required."
                    areopen, whatsopen = self.isOpen()
                    if areopen == True:
                        # truly dire, the telescope is open
                        apflog(lstr, level='Alert', echo=True)
                    else:
                        # telescope powered on, and possibly in the wrong place, but not open
                        apflog(lstr, level='error', echo=True)
                APFTask.waitFor(self.task, True, timeout=30)
            else:
                break
        if result:
            try:
                APFTask.set(self.task, suffix='LAST_CLOSE', value=time.time())
            except:
                apflog("cannot write apftask.MASTER_LAST_CLOSE",level='warn',echo=True)
            return True
        else:
            apflog("Closeup could not successfully complete.")
            return False
        return False

    def updateLastObs(self,obsnum):
        """ If the last observation was a success, this function updates the file storing the last observation number and the hit_list which is required by the dynamic scheduler."""
        if obsnum['populated']:
            APFLib.write(self.robot["MASTER_LAST_OBS_UCSC"], obsnum)
                
        return


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
            if currState == 'enable' and self.wvel <= WINDSHIELD_LIMIT and float(self.temp) > TEMP_LIMIT:
                apflog("Setting scriptobs_windshield to Disable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Disable")
            if currState == 'disable' and (self.wvel > WINDSHIELD_LIMIT or float(self.temp) < TEMP_LIMIT):
                apflog("Setting scriptobs_windshield to Enable")
                APFLib.write(self.robot["SCRIPTOBS_WINDSHIELD"], "Enable")

    def eveningStar(self):
        """Aim the APF at the desired target. This calls prep-obs, slewlock, and focus-telescope. A workaround to relying on scriptobs."""
        if self.isOpen()[0] == False:
            apflog("APF is not open. Can't target a star while closed.",level='error',echo=True)
            return
        self.DMReset()
        # Call prep-obs
        apflog("Calling prep-obs.",echo=True)
        result, ret_code = cmdexec('prep-obs')
        if result == False:
            # try again
            self.DMReset()
            result, ret_code = cmdexec('prep-obs')
            if result is False:
                apflog("Prep-obs returned error code %d. Targeting object has failed." % (ret_code),level='error',echo=True)
                return
            
        self.DMReset()
        apflog("Slewing to lower el",echo=True)
        result, ret_code = cmdexec('slew -e 75')
        if result == False:
            apflog("Slew returned error code %d. Targeting object has failed." % (ret_code),level='error',echo=True)
            return
        # Slew to the specified RA and DEC, set guide camera settings, and centerup( Slewlock )
        # Focus the telescope - all of this, including finding the star, is done in focusTel
        self.DMReset()
        if self.focusTel():
            return True
        else:
            return False
    
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
        robot['MASTER_LAST_OBS_UCSC'].write(obs_num)
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
        return

    def DMReset(self):
        try:
            APFLib.write(self.checkapf['ROBOSTATE'], "master operating",timeout=10)
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
                APFLib.write(self.checkapf['DMTIME'], -1,timeout=10)
        except Exception, e:
            ostr = "Error: cannot touch DM Timer: %s " %( e)
            apflog(ostr,level='error',echo=True)

    def findRobot(self):
        """Trys to find a running instance of robot.csh. Returns the PID along with a boolean representing if the robot was succesfully found."""
        rpid = self.robot['SCRIPTOBS_PID'].read(binary=True)
        if rpid == '' or rpid == -1:
            return rpid, False
        else:
            return rpid, True

    def startRobot(self,observation=None,skip=False,raster=False):
        """Start an instance of scriptobs. Returns the result from subprocess.Popen()."""
        # For running in test mode
        if self.test:
            apflog("Would start robot",echo=True)
            if observation is not None:
                apflog("Would be taking observation in starlist %s" % observation,echo=True)
            APFTask.waitFor(self.task, True, timeout=10)
            return
        
        # Make sure the telescope autofocus is enabled 
        APFLib.write(self.autofoc, "robot_autofocus_enable")
        chk_foc = '$apftask.SCRIPTOBS_AUTOFOC == robot_autofocus_enable'
        result = APFTask.waitFor(self.task, False, chk_foc, timeout=60)
        if not result:
            apflog("Error setting scriptobs_autofoc", level='error',echo=True)
            return
        
        # Make sure APFTEQ is in night mode for observations
        if self.teqmode.read() != 'Night':
            self.setTeqMode('Night')

            # Check the instrument focus for a reasonable value
        if self.dewarfoc > DEWARMAX or self.dewarfoc < DEWARMIN:
            lastfit_dewarfoc = ktl.read("apftask","FOCUSINSTR_LASTFOCUS",binary=True)
            apflog("Warning: The dewar focus is currently %d. This is outside the typical range of acceptable values. Resetting to last derived value %d" % (self.dewarfoc,lastfit_dewarfoc), level = "error", echo=True)
            APFLib.write("apfmot.DEWARFOCRAW",lastfit_dewarfoc)

        robotdir = "/usr/local/lick/bin/robot/"

        telstate = tel['TELSTATE'].read()
        if telstate == 'Disabled':
            rv, retc = cmdexec(os.path.join(robotdir,"slew --hold"))
            if not rv:
                return rv
        rv, retc = cmdexec(os.path.join(robotdir,"prep-obs"))
        if not rv:
            return rv
        # Start scriptobs

        outfile = open("robot.log", 'a')
        if raster:
            args = ['/home/holden/src/raster_scan']
        else:
            if skip:
                args = ['/usr/local/lick/bin/robot/scriptobs', '-dir', os.getcwd(),'-skip']
            else:
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

        ripd, running = self.findRobot()
        if running:
            apflog("Killing Robot %s" % (str(ripd)))
            try:
                APFLib.write(self.robot['SCRIPTOBS_CONTROL'], "abort")
            except Exception, e:
                errstr = "Cannot abort scriptobs: %s" % (e)
                apflog(errstr,level="Warn",echo=True)

    def ucam_powercycle(self, fake=False):

        if fake:
            apflog("would have executed @LROOT/bin/robot/robot_power_cycle_ucam")
            return True
        else:
            val = subprocess.call("/usr/local/lick/bin/robot/robot_power_cycle_ucam")
            if val > 0:
                apflog("power cycle of UCAM failed",level='alert')
                return False

            return True

        return True



    def ucam_reboot(self,fake=False):
        if fake:
            apflog("Would have rebooted UCAM host ",echo=True)
            return True

        apftask = ktl.Service('apftask')
        command = apftask['UCAMLAUNCHER_UCAM_COMMAND']
        ucamstat = apftask['UCAMLAUNCHER_UCAM_STATUS']
        status = apftask['UCAMLAUNCHER_STATUS']
        
        try:
            command.write("Stop")
            apflog("Stopping UCAM software",echo=True)
            
            self.combo_ps.waitFor(" == MissingProcesses",timeout=30)
            command.write("Reboot")
            apflog("Rebooting UCAM host",echo=True)

        except:	      
            apflog("UCAM status bad, cannot restart",level='alert')
            return False
        ucamstat.waitFor(" != running",timeout=60)
        status.waitFor(" != Running",timeout=60)
        status.waitFor(" == Running",timeout=60)
        
        try:
            command.write("Run")
            ucamstat.waitFor(" == running",timeout=300)
            apflog("UCAM software running",echo=True)
            
        except:
            apflog("UCAM status bad, cannot restart",level='alert')
            return False
            
        nv = self.combo_ps.waitFor(" == Ok",timeout=30)
        apflog("UCAM software combo_ps keyword OK",echo=True)
        if nv:
            return nv
        else:
            apflog("UCAM host reboot failure, combo_ps still not ok" , level="alert", echo=True)

        self.obsnum.monitor()
        self.obsnum.callback(self.updateLastObs)

        self.event.monitor()
        self.event.callback(eventmon)


    def ucam_restart(self,fake=False):

        # modify -s apftask UCAMLAUNCHER_UCAM_COMMAND=Stop
        if fake:
            # would have restarted software
            apflog("Would have restarted UCAM software ")
            return True
        else:
            try:
                apflog("Stop and restarting UCAM software",echo=True)
                ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","stop")
                if self.combo_ps.waitFor(" == MissingProcesses",timeout=30):
                    ktl.write("apftask","UCAMLAUNCHER_UCAM_COMMAND","run")
                    nv = self.combo_ps.waitFor(" == Ok",timeout=30)
                    if nv:
                        return nv
                    else:
                        apflog("UCAM  restart failure, combo_ps still not ok" , level="error", echo=True)
            except:	      
                apflog("UCAM status bad, cannot restart",level='alert')
                return False

            self.ucam_reboot()

            
        return False



    def ucam_status(self,fake=False):

        comb = ucam['combo_ps']
        ctalk = ucam['ctalkto']
        if ctalk.read(binary=True) > 0:
            rv = self.ucam_powercycle(fake=fake)
            return rv

        if comb.read(binary=True) > 0:
            # brains!
            rv = self.ucam_restart(fake=fake)
            return rv

        ucamsta = self.apfmon['UCAMSTA'].read(binary=True)
        if ucamsta > 2:
            rv = self.ucam_restart(comb,fake=fake)
            return rv

        return True
    
if __name__ == '__main__':
    print "Testing telescope monitors, grabbing and printing out current state."

    task = 'example'

    APFTask.establish(task, os.getpid())
    apf = APF(task=task,test=False)

    # Give the monitors some time to start up
    APFTask.waitFor(task, True,timeout=2)

    
    print str(apf)

    while True:
        print str(apf)
        APFTask.wait(task,True,timeout=10)


        






