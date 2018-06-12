from __future__ import print_function
import sys
sys.path.append("../master")
import UCSCScheduler_V2 as ds
from  fake_apflog import *
import numpy as np
from optparse import OptionParser
from datetime import datetime
import time
import os
import pickle
import re
import ephem

def parsetemplateGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat"):

    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    try:
        f = open(os.path.join(os.getcwd(),outfn),'r')
    except IOError:
        print ("need to download googledex first")
    else:
        full_codex = pickle.load(f)
        f.close()
        
    col_names = full_codex[0]
    codex = full_codex[1:]

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "lastobs", "owner", \
                "B-V", "APF Desired Precision", "Close Companion", "Template",
                ]

    try:
        idx = [col_names.index(v) for v in req_cols]
    except ValueError:
        print("%s Not found in list" % (v) )
    didx = dict()
    for i,n in enumerate(req_cols):
        didx[n] = idx[i]
    names = []
    star_table = []
    do_flag = []
    stars = []
    template = []
    owner = []
    # Build the star table to return to 
    for ls in codex:
        if ls[0] == '':
            continue
        if float(ls[didx["APFpri"]]) < 0.5: continue
        row = []
        # Get the star name
        names.append(ls[didx["Star Name"]])
        # Get the RA
        row.append(ds.getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]))
        # Get the DEC
        row.append(ds.getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]))
        for i in ("pmRA", "pmDEC", "Vmag"):
            try:
                row.append(float(ls[didx[i]]))
            except ValueError:
                row.append(0.0)
        # For now use the old 1e9 count value
        row.append(900)
        row.append(1.e9)
        for i in ("APFpri", "APFcad"):
            try:
                row.append(float(ls[didx[i]]))
            except ValueError:
                row.append(0.0)
        row.append(7)
        for i in ("lastobs", "B-V", "APF Desired Precision" ):
            try:
                row.append(float(ls[didx[i]]))
            except ValueError:
                row.append(0.0)
        match = re.search("\A(y|Y)",ls[didx["Close Companion"]])
        if match:
            do_flag.append("y")
        else:
            do_flag.append("")
        match = re.search("\A(y|Y)",ls[didx["Template"]])
        if match:
            template.append(True)
        else:
            template.append(False)

        owner.append(ls[didx['owner']])
            
        star_table.append(row)
        star = ephem.FixedBody()
        star._ra = ephem.hours(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]]))
        star._dec = ephem.degrees(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]]))
        stars.append(star)

    return (names, np.array(star_table), do_flag, stars, template, owner)



if __name__ == "__main__":

    dt = datetime.utcfromtimestamp(int(time.time()))

    parser = OptionParser()
    (options, args) = parser.parse_args()    
    allnames, star_table, do_flag, stars, template  = parsetemplateGoogledex()
#    allnames, star_table, do_flag, stars = ds.parseGoogledex()
    for i in range(len(stars)):
        if star_table[i, ds.DS_APFPRI] > 5:
            row = star_table[i,:]
            print (ds.makeScriptobsLine(allnames[i],row,do_flag['do'][i],dt,decker="N",I2="N"),)
            print ("# pri=%.1f" % (star_table[i, ds.DS_APFPRI]))
