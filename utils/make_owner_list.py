import sys
sys.path.append("../master")
import UCSCScheduler_V2 as ds
import numpy as np
from optparse import OptionParser
import datetime
import time
import re

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-o","--owner",dest="owner",default="S.Vogt")
    parser.add_option("-s","--sheet",dest="sheet",default="2017B")
    (options, args) = parser.parse_args()    


    (names, star_table, flags, stars) = ds.parseGoogledex(sheetn=options.sheet)
#    print options.owner 
     
    for i in range(0,len(names)):
#        print flags['owner'][i]
        if options.owner in flags['owner'][i]:
            print names[i]
