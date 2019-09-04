import os
import datetime
import shutil
import threading

import numpy
import ephem

import ktl
import APFTask
import APF


import ParseGoogledex
from SchedulerConsts import * 

IMMEDIATE = 30
NEXT_EXP = 20

class TOO(threading.Thread):
    def __init__(self,sheetns=""):

        threading.Thread.__init__(self)
        self.setDaemon(True)
        
        self.sheetns = sheetns
        
        self.apftask = ktl.Service('apftask')
        self.scriptobs_status = self.apftask['SCRIPTOBS_STATUS']
        self.scriptobs_line = self.apftask['SCRIPTOBS_LINE']
        self.scriptobs_step = self.apftask['SCRIPTOBS_STEP']

        self.apfucam = ktl.Service('apfucam')
        self.event = self.apfucam['EVENT']
        self.stop = self.apfucam['STOP'] # hopefully we will never need this but it will stop an exposure
        

    def end_exposure(self):
        try:
            self.stop.write(1)
            rv = True
        except:
            APF.log(str("Cannot communicate with APFUCAM"), level='error', echo=True)
            rv = False
        return rv

    def end_scriptobs(self, now=False):
        """ In case during an exposure there is a need to stop the robot and close up."""
        APF.log("Terminating scriptobs",echo=True)
        if now:U
            rv = self.end_exposure()
            APF.log("Abort exposure, terminating scriptobs now.",echo=True)
        else:
            if not self.event.read(binary) == 7:
                apflog("Waiting for current exposure to finish.")
                self.event.waitfor(" = ReadoutBegin", timeout=3600)
        APF.log("Killing scriptobs.",echo=True)
        ripd = self.apftask['SCRIPTOBS_PID'].read(binary=True)
        if ripd == -1:
            rv = True
        else:
            try:
                APFLib.write(self.apftask['scriptobs_control'], "abort")
                rv = True
            except Exception, e:
                errstr = "Cannot abort scriptobs: %s" % (e)
                APF.log(errstr,level="Error",echo=True)
                rv = False
        return rv

    
    def run(self):

        while True:
            rv = self.scriptobs_status.waitfor("== Running",timeout=86400)
            if rv:
                rv = self.event.waitfor("== ExposureBegin")
                # during an exposure get ToO sheet
                if rv: 
                    names, star_table, flags, stars = ParseGoogledex.parseGoogledex(sheetns=self.sheetns,outfn='too.dat',force_download=True)
                    dt = datetime.datetime.utcnow()
                    good_cadence = (ephem.julian_date(dt) - star_table[:, DS_LAST]) > star_table[:, DS_CAD]
                    immediate = (star_table[:,DS_APFoPRI] > IMMEDIATE) & good_cadence
                    next_exp = (star_table[:,DS_APFPRI] < IMMEDIATE) & (star_table[:,DS_APFPRI] >NEXT_EXP) & good_cadence

                    if len(star_table[:,DS_APFPRI][immediate]) > 0:
                        rv = self.end_exposure()
                        rv = self.end_scriptobs(now=True)
                        if rv is False:
                            APF.log("Cannot abort scriptobs which means all kinds of bad things are happening",level='error',echo=True)
                    if len(star_table:,DS_APFPRI][next_exp]) > 0:
                        rv self.end_scriptobs()
                        if rv is False:
                            APF.log("Cannot abort scriptobs which means all kinds of bad things are happening",level='error',echo=True)                    
            

