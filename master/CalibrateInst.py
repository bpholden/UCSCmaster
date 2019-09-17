class CalibrateInst(threading.Thread):
    def __init__(self, apf, opt, inst='Levy', task='master'):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.apf = apf
        self.task = task
        self.user = opt.name
        self.owner = 'public'
        self.signal = True
        self.debug = opt.test
        self.inst = inst


    def focus(self):

        result = self.apf.focusinstr()
        if result is False:
            apflog("Focus has failed.")
        else:
            apflog("Focus has finished. Setting phase to Cal-Pre")
        if not self.debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

        return result

    def calibrate(self):
        self.apf.instr_permit()

        if self.inst == 'Levy':
            result = self.apf.ucam_status()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
                return result
            
            result = self.apf.calibrate(script = opt.calibrate, time = 'pre')
            if not debug:
                APFTask.set(parent, suffix="LAST_OBS_UCSC", value=self.apf.ucam["OBSNUM"].read())

        if result == False:
            apflog("Calibrate has failed. Trying again",level='warn',echo=True)
            self.apf.instr_permit()
            if self.inst == 'Levy':
                result = self.apf.ucam_status()
                if result is False:
                    apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
                    return result
                result = self.apf.calibrate(script = opt.calibrate, time = 'pre')
                if not result:
                    apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                    self.apf.turn_off_lamps()
        return result
    
    def run(self):
        """ CalibrateInst.run() - runs the calibration
        """

        apflog("Starting Instrument Focus script.", level='Info', echo=True)
        result = self.focus()

        if result:
            apflog("Starting calibrate script.", level='Info', echo=True)
            result = self.calibrate()

        return result
            
