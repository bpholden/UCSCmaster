import re
try:
    from apflog import *
except:
    from fake_apflog import *


class ObservedLog():
    def __init__(self):
        self.names = []
        self.times = []
        self.tempobs = []
        self.filename = ""
    def read_observed_log(self,filename):
        """ read_observed_log parses a file to find the object names and times
        ObservedLog.read_observed_log(filename)
        names - list of names, must be first column of file called filename
        times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line
        temps - a list of template observations

        """

        try:
            f = open(filename, 'r')
        except IOError:
            apflog( "Couldn't open %s" % filename,level="warn",echo=True)
            return obs, times
        else: 
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == '#' or line == "":
                        pass
                    else:
                        ls = line.split()
                        self.obs.append(ls[0])
                        if len(ls) > 15:
                            self.times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
                        else:
                            self.times.append(float(ls[1]))
                        length = len(ls)
                        mtch = re.search("temp\=Y",ls[length-1])
                        if mtch:
                            self.temps.append("Y")
                        else:
                            self.temps.append("N")

            
        self.obs.reverse()
        self.times.reverse()
        self.temps.reverse()
        return 
        

def getObserved(filename):
    """ getObserved parses a file to find the object names and times
    names, times, temps = getObserved(filename)
    names - list of names, must be first column of file called filename
    times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line
    temps - a list of template observations

    """
    obs = []
    times = []
    temps = []
    nobs = dict()
    try:
        f = open(filename, 'r')
    except IOError:
        apflog( "Couldn't open %s" % filename,level="warn",echo=True)
        return obs, times
    else: 
        for line in f:
            line = line.strip()
            if len(line) > 0:
                if line[0] == '#' or line == "":
                    pass
                else:
                    ls = line.split()
                    obs.append(ls[0])
                    if len(ls) > 15:
                        times.append( (int(ls[14].split('=')[1]), int(ls[15].split('=')[1])) )
                    else:
                        times.append(float(ls[1]))
                    length = len(ls)
                    mtch = re.search("temp\=Y",ls[length-1])
                    if mtch:
                        temps.append("Y")
                    else:
                        temps.append("N")

            
    obs.reverse()
    times.reverse()
    temps.reverse()
    return obs, times, temps
	
