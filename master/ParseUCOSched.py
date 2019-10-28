from __future__ import print_function
from datetime import datetime, timedelta
import os
import re
import pickle
import sys
import time
import ephem
import numpy as np

import gspread
import json
from oauth2client.service_account import ServiceAccountCredentials

import ObservedLog
import Coords
from SchedulerConsts import MIN_TOTOBS, DS_BV, DS_ERR
import ExposureCalculations as ec

try:
    from apflog import *
except:
    from fake_apflog import *

def checkFlag(key,didx,line,regexp,default):
    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default


def parseStarname(starname):

    ostarname = starname.strip()
    m= re.search("HD\s+\d+",starname)
    if m:
        ostarname = re.sub("HD\s+","HD",starname)
    m = re.search("\s+",ostarname)
    while m:
        ostarname = re.sub("\s+","_",ostarname)
        m = re.search("\s+",ostarname)
    m = re.search("\+",ostarname)
    while m:
        ostarname = re.sub("\+","p",ostarname)

        
    return ostarname

    


def parseUCOSched(sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="sched.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '' },force_download=False):
    """ parseUCOSched parses google sheets and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseUCOSched(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="sched.dat")
    names - a list of stars in the starlist
    star_table - a numpy array
    flags - a dictionary of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "texp", "I2", "expcount", "decker","Close Companion", "APFnshots", \
                "owner", \
                "APFpri", "APFcad", "lastobs", "APFmax", "B-V", \
                "uth","utm","duration", "Template", "Nobs", "Total Obs"
                ]

    
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    apflog( "Starting Googledex parse",echo=True)    
    if not outdir :
        outdir = os.getcwd()
    if os.path.exists(os.path.join(outdir,outfn)) and force_download is False:
        try:
            f = open(os.path.join(outdir,outfn),'rb')
            full_codex = pickle.load(f)
            f.close()
        except:
            full_codex = make_local_copy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))
    else:
        full_codex = make_local_copy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))

    col_names = full_codex[0]
    codex = full_codex[1:]

    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    flags = { "do" : [], "decker" : [], "I2" : [], "owner" : [], "template" : [] }
    stars = []
    # Build the star table to return to 
    for ls in codex:
        row = []
        if ls[0] == '':
            continue
        apfpri = float_or_default(ls[didx["APFpri"]])
        nobs = int_or_default(ls[didx["Nobs"]])
        totobs = int_or_default(ls[didx["Total Obs"]],default=-1)

        if totobs > 0 and nobs >= totobs: continue
        if apfpri < 0.5: continue
        # Get the star name
        names.append(parseStarname(ls[didx["Star Name"]]))
        
        # Get the RA
        raval = Coords.getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval:
            row.append(raval)
        else:
            row.append(-1.)
        # Get the DEC
        decval = Coords.getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
        if decval:
            row.append(decval)
        else:
            row.append(-3.14)

        for coln in ("pmRA", "pmDEC"):
            row.append(float_or_default(ls[didx[coln]]))

        # Vmag
        row.append(float_or_default(ls[didx["Vmag"]],default=15.0))
        row.append(float_or_default(ls[didx["texp"]],default=1200))
        row.append(float_or_default(ls[didx["expcount"]],default=1e9))
        row.append(int_or_default(ls[didx["APFnshots"]],default=1))


        # scheduler specific
        row.append(apfpri)
        row.append(float_or_default(ls[didx["APFcad"]],default=0.7))
        row.append(float_or_default(ls[didx["lastobs"]],default=0))
        # APFmax
        row.append(float_or_default(ls[didx["APFmax"]]))

        inval = float_or_default(ls[didx["B-V"]],default=0.7)
        if inval < 0:
            inval = 1.
        if coln is 'B-V' and inval > 2:
            inval = 1
        row.append(inval)
                    
        for coln in ["uth", "utm"]:
            row.append(int_or_default(ls[didx[coln]]))
                
        # duration:
        row.append(float_or_default(ls[didx["duration"]]))
                
        # Nobs
        row.append(nobs)
                
        # Total Obs
        if totobs >= 0:
            row.append(totobs)
        else:
            row.append(0)

        check = checkflag("Close Companion",didx,ls,"\A(y|Y)","")
        if check == "Y" or check == "y" :
            flags['do'].append(check)
        else:
            flags['do'].append("")
            
        flags['decker'].append(checkflag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        i2select = checkflag("I2",didx,ls,"\A(n|N)",config["I2"])
        flags['I2'].append(i2select.upper())
        tempselect = checkflag("Template",didx,ls,"\A(n|N)",'Y')
        flags['template'].append(tempselect.upper())

        flags['owner'].append(checkflag("owner",didx,ls,"\A(\w?\.?\w+)",config["owner"]))

            
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(str(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]])))
        star._dec = ephem.degrees(str(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]])))
        stars.append(star)

    return (names, np.array(star_table), flags, stars)
