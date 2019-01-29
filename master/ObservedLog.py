try:
    from apflog import *
except:
    from fake_apflog import *


def getObserved(filename):
    """ getObserved parses a file to find the object names and times
    names, times = getObserved(filename)
    names - list of names, must be first column of file called filename
    times - times either as a timestamp in second column or a (hour,minute) tuple from a scriptobs line

    """
    obs = []
    times = []
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
            
    obs.reverse()
    times.reverse()
    return obs, times
	
