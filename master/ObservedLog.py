import re
try:
    from apflog import *
except:
    from fake_apflog import *


class ObservedLog():
    def __init__(self,filename=''):
        self.names = []
        self.times = []
        self.tempobs = []
        self.owners = []        
        self.filename = filename


    def parse_key_vals(line_list):
        keyvals = dict()
        ovals = []
        for ls in line_list:
            split_vals = ls.split('=')
            if len(split_vals) == 2:
                keyvals[split_vals[0]] = split_vals[1]
            else:
                ovals.append(split_vals[0])

        return ovals, keyvals
        
        
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
            return (), (), (), ()
        else: 
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == '#' or line == "":
                        pass
                    else:
                        ovals, keyvals = self.parse_key_vals(line.split())
                        self.obs.append(ovals[0])
                        if 'uth' in keyvals.keys():
                            self.times.append( ( int(keyvals['uth']), int(keyvals['utm']) ) )
                        else:
                            self.times.append(float(ovals[1]))

                        
                        if 'temp' in keyvals.keys():
                            self.temps.append("Y")
                        else:
                            self.temps.append("N")
                            
                        if 'owner' in keyvals.keys():
                            self.owners.append(keyvals['owner'])
            
        self.obs.reverse()
        self.times.reverse()
        self.temps.reverse()
        self.owners.reverse()
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
	
