from __future__ import print_function
from datetime import datetime, timedelta
import os
import re
import pickle
import sys
import time
import ephem
import numpy as np
import astropy
import astropy.table
import astropy.io.ascii

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


def int_or_default(value,default=0):
    try:
        attr = int(value)
    except:
        attr = default
    return attr

def float_or_default(value,default=0.0):
    try:
        rv = float(value)
    except:
        rv = default
    return rv

def makeLocalCopy(req_cols,sheetns=["The Googledex"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="./googledex.dat"):
    full_codex = []
    # These are the columns we need for scheduling
    full_codex.append(req_cols)
        
    for sheetn in sheetns:
        worksheet = getSpreadsheet(sheetn=sheetn,certificate=certificate)
        if worksheet:
            cur_codex = worksheet.get_all_values()
            if len(cur_codex) <= 0:
                apflog("Worksheet %s exists but is empty, skipping" % (sheetn), level='error', echo=True)

                continue
            didx = findColumns(cur_codex[0],req_cols)
            
            for row in cur_codex[1:]:
                nrow = []
                for c in req_cols:
                    if c in didx.keys():
                        nrow.append(row[didx[c]])
                    else:
                        nrow.append(None)
                full_codex.append(nrow)
        
    f = open(outfn,'wb')
    pickle.dump(full_codex, f)
    f.close()
                
    return full_codex
    

def getSpreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json'):
    """ Get the spreadsheet from google

    worksheet = getSpreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')
    worksheet - the worksheet object returned by the gspread module

    sheetn - name of the google sheet, defaults to "The Googledex"
    certificate - certificate used to control access to the google sheet
    
    """
    # this downloads the googledex from the Google Drive
    # the certificate must be available
    # these certificates are generated through the Google Developer Interface
    # the developer must select the correct API for access

    # the certificate has an email associated with it, that email must
    # have the document shared with it to allow access 

    certificate_path = os.path.dirname(__file__)
    
    json_key = json.load(open(os.path.join(certificate_path, certificate)))
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']

    credentials = ServiceAccountCredentials.from_json_keyfile_name(os.path.join(certificate_path, certificate), scope)
    try:
        gs = gspread.authorize(credentials)
        apflog("Successfully logged in.", echo=True)
    except:
        apflog("Cannot log into Google API.", echo=True,level='error')
        return None
    apflog("Attempting to Open %s" % (sheetn),echo=True)
    try:
        spreadsheet = gs.open(sheetn)
        apflog("Loaded Main %s" % (sheetn),echo=True)
        worksheet = spreadsheet.sheet1
        apflog("Got spreadsheet", echo=True)
    except Exception as e:
        apflog("Cannot Read %s: %s"  % (sheetn, e), echo=True, level='error')
        worksheet = None
    return worksheet


def findColumns(col_names,req_cols,opt_cols=[]):
    """ findColumns finds the indices for the column names in the list of required columns
    indices = findColumns(col_names, req_cols)
    
    indices - a list of indices, each index maps to where in col_names the column is found and in the order of req_cols
    col_names - list of column names to be searched
    req_cols - list of names that should be in the first list
    """
    idx = []
    didx = dict()

    for r in req_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
        else:
            apflog("%s Not found in column names from google spreadsheet" % (r) , level="Warn",echo=True)

    for r in opt_cols:
        if r in col_names:
            didx[r] = col_names.index(r)
            
    # hack to handle an error
    if req_cols[0] == "Star Name" and req_cols[0] not in didx.keys():
        didx[req_cols[0]] = 0
        apflog("Pasting 'Star Name' into column 0 of google spreadsheet" , level="Warn",echo=True)

    return didx


def parseFracTable(sheet_table_name='2019B_frac',certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn=None):
    
    apflog( "Starting parse of %s" % (sheet_table_name),echo=True)    
    if not outdir :
        outdir = os.getcwd()
    if outfn is not None and os.path.exists(os.path.join(outdir,outfn)):
        sheetns=[]
        frac=[]
        with open(os.path.join(outdir,outfn)) as fp:
            lines = fp.readlines()
            for ln in lines:
                row = ln.strip().split()
                sheetns.append(row[0])
                try:
                    frac.append(float(row[0]))
                except:
                    frac.append(0)
        return sheetns,frac

    sheetns = []
    frac = []
    twod = []
    worksheet = getSpreadsheet(sheetn=sheet_table_name,certificate=certificate)
    if worksheet:
        cur_codex = worksheet.get_all_values()
        if len(cur_codex) <= 0:
            apflog("Worksheet %s exists but is empty, skipping" % (sheetn), level='error', echo=True)
            return None, None
        for row in cur_codex:
            sheetns.append(row[0])
            frac.append(float_or_default(row[1]))
            twod.append([row[0],row[1]])

        if outfn is not None:
            np.savetxt(os.path.join(outdir,outfn),twod,fmt="%s",delimiter=" ")
            
    return sheetns,frac

def parseUCOSched(sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="sched.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '', 'mode' : '', 'obsblock' : '', 'Bstar' : 'N' , 'raoff' : None, 'decoff' : None },force_download=False):
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
                    "texp", "I2", "expcount","decker","Close Companion", "APFnshots", \
                    "owner", "mode", "raoff", "decoff", \
                    "Bstar", "obsblock",\
                    "APFpri", "APFcad", "lastobs", "B-V", \
                    "uth","utm","duration", \
                    "Template", "Nobs", "Total Obs"
                    ]
    
    
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    if not outdir :
        outdir = os.getcwd()
    if os.path.exists(os.path.join(outdir,outfn)) and force_download is False:
        try:
            f = open(os.path.join(outdir,outfn),'rb')
            full_codex = pickle.load(f)
            f.close()
        except:
            full_codex = makeLocalCopy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))
    else:
        full_codex = makeLocalCopy(req_cols,sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))

    col_names = full_codex[0]
    codex = full_codex[1:]

    didx = findColumns(col_names,req_cols)
    
    star_table = { "name" : [], "ra" : [], 'dec' : [], 'pmRA' : [], 'pmDec' : [], 'Vmag' : [], 'texp' : [], 'expcount' : [], 'APFnshots' : [], 'APFpri' : [], 'APFcad' : [], 'lastobs' : [], 'BmV' : [], 'uth' : [], 'utm' : [], 'duration' : [], 'nobs' : [], 'totobs' : [])
    
    flags = { "do" : [], "decker" : [], "I2" : [], "owner" : [], "template" : [], "obsblock" : [], "mode" : [], "Bstar" : [], "inst" : [], "raoff" : [], "decoff" : [] }
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
        #        names.append(parseStarname(ls[didx["Star Name"]]))
        star_table['name'].append(parseStarname(ls[didx["Star Name"]]))
        # Get the RA
        raval = Coords.getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval:
            star_table['ra'].append(raval)
        else:
            star_table['ra'].append(-1.)
        # Get the DEC
        decval = Coords.getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
        if decval:
            star_table['dec'].append(decval)
        else:
            star_table['dec'].append(-3.14)

        for coln in ("pmRA", "pmDEC"):
            star_table[coln].append(float_or_default(ls[didx[coln]]))


        star_table['Vmag'].append(float_or_default(ls[didx["Vmag"]],default=15.0))
        star_table['texp'].append(float_or_default(ls[didx["texp"]],default=1200))
        expcount = float_or_default(ls[didx["expcount"]],default=1e9)
        if expcount > 3e9:
            expcount = 3e9
        star_table['expcount'].append(expcount)
        star_table['APFnshots'].append(int_or_default(ls[didx["APFnshots"]],default=1))


        # scheduler specific
        star_table['APFpri'].append(apfpri)
        star_table['APFcad'].append(float_or_default(ls[didx["APFcad"]],default=0.7))
        star_table["lastobs"]].append(float_or_default(ls[didx["lastobs"]],default=0))

        inval = float_or_default(ls[didx["B-V"]],default=0.7)
        if inval < 0:
            inval = 1.
        if coln is 'B-V' and inval > 2:
            inval = 1
        star_table['BmV'].append(inval)
                    
        for coln in ["uth", "utm"]:
            star_table[coln].append(int_or_default(ls[didx[coln]]))
                
        # duration:
        star_table['duration'].append(float_or_default(ls[didx["duration"]]))
                
        # Nobs
        star_table['nobs'].append(nobs)
                
        # Total Obs
        if totobs >= 0:
            star_table['totobs'].append(totobs)
        else:
            star_table['totobs'].append(0)

        check = checkFlag("Close Companion",didx,ls,"\A(y|Y)","")
        if check == "Y" or check == "y" :
            flags['do'].append(check)
        else:
            flags['do'].append("")
            
        flags['decker'].append(checkFlag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        i2select = checkFlag("I2",didx,ls,"\A(n|N)",config["I2"])
        flags['I2'].append(i2select.upper())
        tempselect = checkFlag("Template",didx,ls,"\A(n|N)",'Y')
        flags['template'].append(tempselect.upper())

        flags['owner'].append(checkFlag("owner",didx,ls,"\A(\w?\.?\w+)",config["owner"]))
        flags['mode'].append(checkFlag("mode",didx,ls,"\A(b|B|o|O)",config["mode"]).upper())
        flags['obsblock'].append(checkFlag("obsblock",didx,ls,"\A(\w+)",config["obsblock"]))
        flags['Bstar'].append(checkFlag("Bstar",didx,ls,"\A(Y|y)",config["Bstar"]).upper())
        flags['inst'].append(checkFlag("inst",didx,ls,"(levy|darts)",config['inst']).lower())

        flags['raoff'].append(checkFlag("raoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["raoff"]))
        flags['decoff'].append(checkFlag("decoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["decoff"]))
        
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(str(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]])))
        star._dec = ephem.degrees(str(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]])))
        stars.append(star)

    star_table =  astropy.table.Table(star_table,names=['name','ra','dec','pmRA','pmDec','Vmag','texp','expcount','APFnshots','APFpri','APFcad','lastobs','BmV','uth','utm','duration','nobs','totobs'])
    star_table = astropy.table.hstack([star_table,astropy.table.Table(flags)])
    return (star_table, stars)


def updateLocalGoogledex(intime,googledex_file="googledex.dat", observed_file="observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time.
        update_local_googledex(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """
    # name, times, temps = ObservedLog.getObserved(observed_file)
    obslog = ObservedLog.ObservedLog(observed_file)
    try:
        g = open(googledex_file, 'rb')
        full_codex = pickle.load(g)
        g.close()
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
        return obslog.names
    except EOFError:
        apflog("googledex file corrupt, so can't be updated",echo=True)
        return obslog.names


    codex_cols = full_codex[0]

    starNameIdx = codex_cols.index("Star Name")
    lastObsIdx = codex_cols.index("lastobs")
    try:
        nObsIdx = codex_cols.index("nObsIdx")
    except:
        nObsIdx = -1
    
    for i in range(1, len(full_codex)):
        row = full_codex[i]
        if row[starNameIdx] in obslog.names:
            # We have observed this star, so lets update the last obs field
            obstime = obslog.times[obslog.names.index(row[starNameIdx])]
            if isinstance(obstime,float):
                t = datetime.utcfromtimestamp(obstime)
            else:
                hr, min = obstime
                if type(intime) != datetime:
                    ctime = datetime.now()
                    td = timedelta(0,3600.*7)
                    intime = ctime + td
                t = datetime(intime.year, intime.month, intime.day, hr, min)


            jd = round(float(ephem.julian_date(t)), 4) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)),echo=True)
            row[lastObsIdx] = str(jd)
            if nObsIdx > 0:
                row[nObsIdx] = row[nObsIdx] + 1
            full_codex[i] = row

    with open(googledex_file, 'wb') as f:
        pickle.dump(full_codex, f)
    f.close()
    
    return obslog.names

