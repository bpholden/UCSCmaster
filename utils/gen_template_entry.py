#!/usr/bin/env  /opt/kroot/bin/kpython
import sys
sys.path.append("../master")
#from ExposureCalc import *
import UCSCScheduler_V2 as ds
import numpy as np
from optparse import OptionParser
from datetime import datetime

if __name__ == "__main__":

    parser = OptionParser()
    (options, args) = parser.parse_args()    
    if len(args) < 1:
        print "needs a name"
        sys.exit()

     
    allnames, star_table, do_flag, stars  = ds.parseGoogledex()

   
    dt = datetime.utcnow()
    for arg in args:
        if arg in allnames:
            i = allnames.index(arg)
            row = []

            row.append(star_table[i, ds.DS_RA])
            row.append( star_table[i, ds.DS_DEC])
            row.append(star_table[i, ds.DS_PMRA])
            row.append(star_table[i, ds.DS_PMDEC])
            row.append(star_table[i, ds.DS_VMAG])
            row.append(1200)
            row.append(1e9)
            row.append( star_table[i, ds.DS_APFPRI])
            row.append(0)
            if star_table[i, ds.DS_VMAG] > 10:
                row.append(9)
            else:
                row.append(5)                
            line = ds.makeScriptobsLine(allnames[i],row,do_flag['do'][i],dt,decker="N",I2="N")
            print "%s #  %s" % (line,"pri = %s" % (star_table[i, ds.DS_APFPRI]))

