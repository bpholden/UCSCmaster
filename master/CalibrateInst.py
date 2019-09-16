class CalibrateInst(threading.Thread):
    def __init__(self, apf, opt, inst='Levy', task='master'):
        threading.Thread.__init__(self)
        self.setDaemon(True)
        self.APF = apf
        self.task = task
        self.user = opt.name
        self.owner = 'public'
        self.signal = True
        self.debug = opt.test
        self.inst = inst
        
    def run(self):
        """ CalibrateInst.run() - runs the observing
        """

        apflog("Starting focusinstr script.", level='Info', echo=True)
        apf.instr_permit()
        if self.inst == 'Levy':
            result = apf.ucam_status()
            if result is False:
                apflog("Failure in UCAM status and restart!", level='Alert', echo=True)

        result = apf.focusinstr()
        
        apflog("Focus has finished. Setting phase to Cal-Pre")
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())
            
        apflog("Starting calibrate pre script.", level='Info', echo=True)
        apf.instr_permit()
        result = apf.ucam_status()
        if result is False:
            apflog("Failure in UCAM status and restart!", level='Alert', echo=True)
        result = apf.calibrate(script = opt.calibrate, time = 'pre')
        if not debug:
            APFTask.set(parent, suffix="LAST_OBS_UCSC", value=apf.ucam["OBSNUM"].read())

        if result == False:
            apflog("Calibrate Pre has failed. Trying again",level='warn',echo=True)
            apf.instr_permit()
            result = apf.calibrate(script = opt.calibrate, time = 'pre')
            if not result:
                apflog("Error: Calibrate Pre has failed twice. Observer is exiting.",level='error',echo=True)
                apf.turn_off_lamps()

            
