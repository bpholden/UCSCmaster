from __future__ import print_function
from datetime import datetime, timedelta
import os
import re

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
from SchedulerConsts import MIN_TOTOBS, DS_BV, DS_ERR, EXP_LIM
import ExposureCalculations as ec

try:
    from apflog import *
except:
    from fake_apflog import *

def checkFlag(key,didx,line,regexp,default):
    """ checkFlag(key, dict_ind, line, regexp, default)
    
    - key : the key in the dictionary indices, which key from the input list you want to check
    - didx : dictionary of indices so that keys can be used instead of indices for the line
    - line : a list of entries, map to keys by didx
    - regexp : if the regexp is matched, that value is returned as a string
    - default : value if the regexp fails

    """
    
    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default


def parseStarname(starname):
    """parseStarname(starname)

    starname - input value which should be the name of the a star, duh
    returns the starname value, doing some clean up to meet requirements (no spaces, trim trialing spaces, etc.) 
    """
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


def intDefault(value,default=0):
    """
    intDefault(value,default=0)
    returns the input value as an integer, and a failure to cast returns default
    """
    try:
        attr = int(value)
    except:
        attr = default
    return attr

def floatDefault(value,default=0.0):
    """
    floatDefault(value,default=0.0)
    returns the input value as an float, and a failure to cast returns default
    """
    try:
        rv = float(value)
    except:
        rv = default
    return rv


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
    worksheet = None
    tries = 0
    while worksheet is None and tries < 3:
        tries = tries + 1
        try:
            spreadsheet = gs.open(sheetn)
            apflog("Loaded Main %s" % (sheetn),echo=True)
            worksheet = spreadsheet.sheet1
            apflog("Got spreadsheet", echo=True)
        except Exception as e:
            apflog("Cannot Read %s: %s"  % (sheetn, e), echo=True, level='error')
            time.sleep(1)
    return worksheet

def retrieveCodex(req_cols,sheetns=["The Googledex"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json'):
    """retrieveCodex(req_cols,sheetns=["The Googledex"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')

    returns the "codex", a list of lists containing all of the columns
    in the req_cols list, source of the data are the Google sheets named
    in sheetns, needs a certificate to authenticate.

    - req_cols : a list of column names in the sheets that are required
    for the final list of lists
    - sheetns : sheets to download from
    - certificate : the thing that allows authentication
    """
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
                        if c is 'sheetn':
                            nrow.append(sheetn)
                        else:
                            nrow.append(None)

                full_codex.append(nrow)
        
                
    return full_codex
    

def findColumns(col_names,req_cols,opt_cols=[]):
    """findColumns finds the indices for the column names in the list of

    required columns indices = findColumns(col_names, req_cols)

    indices - a dictionary of indices, each index maps to where in
    col_names the column is found 

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


def parseFracTable(sheet_table_name='2019B_frac',certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn=None,outdir=None):
    
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
            frac.append(floatDefault(row[1]))
            twod.append([row[0],row[1]])

            
    return sheetns,frac

def parseRankTable(sheet_table_name='2019B_ranks',certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json'):
    
    apflog( "Starting parse of %s" % (sheet_table_name),echo=True)    

    sheetns = []
    rank = []

    worksheet = getSpreadsheet(sheetn=sheet_table_name,certificate=certificate)
    if worksheet:
        cur_codex = worksheet.get_all_values()
        if len(cur_codex) <= 0:
            apflog("Worksheet %s exists but is empty, skipping" % (sheetn), level='error', echo=True)
            return None, None
        for row in cur_codex[1:]:
            sheetns.append(row[0])
            rank.append(floatDefault(row[1]))
            
    return sheetns,rank


def initStarTable(col_list):

    """
    star_table = initStarTable(column_list)
    star_table - a Astropy Table object that has the columns needed, most are in column_list
    forces certain columns to be added


    """
    
# star_table = { "name" : [], "ra" : [], 'dec' : [], 'pmRA' : [], 'pmDEC' : [], 'Vmag' : [], 'texp' : [], 'expcount' : [], 'APFnshots' : [], 'APFpri' : [], 'APFcad' : [], 'lastobs' : [], 'BmV' : [], 'uth' : [], 'utm' : [], 'duration' : [], 'nobs' : [], 'totobs' : [], "do" : [], "decker" : [], "I2" : [], "owner" : [], "template" : [], "obsblock" : [], "mode" : [], "Bstar" : [],  "raoff" : [], "decoff" : [], 'sheetn' : [] }

    star_table = dict()
    for col in col_list:
        star_table[col] = []
    star_table['name'] = []
    star_table['do'] = []
    star_table['nobs'] = []
    star_table['totobs'] = []
    star_table['ra'] = []
    star_table['dec'] = []


    return star_table


def parseCodex(config,sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',prilim=0.5):
    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                    "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                    "texp", "I2", "expcount","decker","Close Companion", "APFnshots", \
                    "owner", "APFpri", "APFcad", "lastobs", "B-V", \
                    "uth","utm","duration", \
                    "Template", "Nobs", "Total Obs", \
                    "mode", "raoff", "decoff", "Bstar", "obsblock",\
                    'sheetn' \
                    ]

    full_codex = retrieveCodex(req_cols,sheetns=sheetns,certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')

    col_names = full_codex[0]
    codex = full_codex[1:]

    didx = findColumns(col_names,req_cols)
    star_table = initStarTable(req_cols)
    
    stars = []
    # Build the star table to return to 
    for ls in codex:
        row = []
        if ls[0] == '':
            continue
        apfpri = floatDefault(ls[didx["APFpri"]])
        nobs = intDefault(ls[didx["Nobs"]])
        totobs = intDefault(ls[didx["Total Obs"]],default=-1)

        if totobs > 0 and nobs >= totobs: continue
        if apfpri < prilim: continue
        # Get the star name
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

        for coln in ("RA hr","RA min","RA sec","Dec deg","Dec min","Dec sec"):
            star_table[coln].append(ls[didx[coln]])

        for coln in ("pmRA", "pmDEC"):
            star_table[coln].append(floatDefault(ls[didx[coln]]))


        star_table['Vmag'].append(floatDefault(ls[didx["Vmag"]],default=15.0))
        star_table['texp'].append(floatDefault(ls[didx["texp"]],default=1200))
        expcount = floatDefault(ls[didx["expcount"]],default=1e9)
        if expcount > EXP_LIM:
            expcount = EXP_LIM
        star_table['expcount'].append(expcount)
        star_table['APFnshots'].append(intDefault(ls[didx["APFnshots"]],default=1))


        # scheduler specific
        star_table['APFpri'].append(apfpri)
        star_table['APFcad'].append(floatDefault(ls[didx["APFcad"]],default=0.7))
        star_table["lastobs"].append(floatDefault(ls[didx["lastobs"]],default=0))

        inval = floatDefault(ls[didx["B-V"]],default=0.7)
        if inval < 0:
            inval = 1.
        if coln is 'B-V' and inval > 2:
            inval = 1
        star_table['B-V'].append(inval)
                    
        for coln in ["uth", "utm"]:
            star_table[coln].append(intDefault(ls[didx[coln]]))
                
        # duration:
        star_table['duration'].append(floatDefault(ls[didx["duration"]]))
                
        # Nobs
        star_table['nobs'].append(nobs)
                
        # Total Obs
        if totobs >= 0:
            star_table['totobs'].append(totobs)
        else:
            star_table['totobs'].append(0)

        check = checkFlag("Close Companion",didx,ls,"\A(y|Y)","")
        if check == "Y" or check == "y" :
            star_table['do'].append(check)
        else:
            star_table['do'].append("")
            
        star_table['decker'].append(checkFlag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        i2select = checkFlag("I2",didx,ls,"\A(n|N)",config["I2"])
        star_table['I2'].append(i2select.upper())
        tempselect = checkFlag("Template",didx,ls,"\A(n|N)",'Y')
        star_table['Template'].append(tempselect.upper())

        star_table['owner'].append(checkFlag("owner",didx,ls,"\A(\w?\.?\w+)",config["owner"]))
        star_table['mode'].append(checkFlag("mode",didx,ls,"\A(b|B|o|O)",config["mode"]).upper())
        star_table['obsblock'].append(checkFlag("obsblock",didx,ls,"\A(\w+)",config["obsblock"]))
#        star_table['inst'].append(checkFlag("inst",didx,ls,"(levy|darts)",config['inst']).lower())

        star_table['raoff'].append(checkFlag("raoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["raoff"]))
        star_table['decoff'].append(checkFlag("decoff",didx,ls,"\A((\+|\-)?\d+\.?\d*)",config["decoff"]))

        csheetn = checkFlag("sheetn",didx,ls,"\A(.*)",'public')

        if 'Bstar' in csheetn :
            star_table['Bstar'].append("Y")
            star_table['sheetn'].append('RECUR_A100')
        else:
            star_table['Bstar'].append("N")
            star_table['sheetn'].append(csheetn)

    badkeylist = []
    for k in star_table.keys():
        if len(star_table[k]) == 0:
            badkeylist.append(k)
    for k in badkeylist:
        del star_table[k]  

    return star_table

def genStars(star_table):
    """pyephem_objs = genStars(star_table)

    given a star_table returned by parseCodex (or initStarTable) returns
    a list of pyephem objects for every object in the table

    Inputs star_table - astropy Table that must have the RA and Dec in
    sexigresimal format with each column for each part of the
    coordinates separate
    """
    stars = []
    for i in range(0,len(star_table['name'])):
        star = ephem.FixedBody()
        star.name = star_table['name'][i]
        star._ra = ephem.hours(str(":".join([star_table["RA hr"][i], star_table["RA min"][i], star_table["RA sec"][i]])))
        star._dec = ephem.degrees(str(":".join([star_table["Dec deg"][i], star_table["Dec min"][i], star_table["Dec sec"][i]])))
        stars.append(star)
        
    return stars



def parseUCOSched(sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn="sched.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '', 'mode' : '', 'obsblock' : '', 'Bstar' : 'N' , 'raoff' : None, 'decoff' : None },force_download=False,prilim=0.5):
    """ parseUCOSched parses google sheets and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    star_table, stars = parseUCOSched(sheetns=["Bstars"],certificate='cert.json',outfn="sched.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '', 'mode' : '', 'obsblock' : '', 'Bstar' : 'N' , 'raoff' : None, 'decoff' : None },force_download=False,prilim=0.5)

    star_table - an astropy table
    stars - a list of pyEphem objects 

    Inputs:
    sheetns - list of google sheet names
    certificate - json file for authenticating to access google sheets
    outfn - output file name, will read this in if it already exists instead of downloading sheets if force_download is False
    outdir - output directory for outfn, defaults to ./
    config - default values for a number of flags
    force_download - force the google sheets to be downloaded even if outfn already exists
    prilim - limit on priority values, values below this are tossed

    """

    
    
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    if not outdir :
        outdir = os.getcwd()
    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn) and force_download is False:
        try:
            star_table = astropy.io.ascii.read(outfn)
        except:
            star_table  = parseCodex(config,sheetns=sheetns,certificate=certificate,prilim=prilim)

    else:
        star_table = parseCodex(config,sheetns=sheetns,certificate=certificate,prilim=prilim)

    stars = genStars(star_table)
    
    star_table = astropy.table.Table(star_table)
    astropy.io.ascii.write(star_table,outfn, format='ecsv', overwrite=True)
    
    return (star_table, stars)


def updateLocalStarlist(intime, observed_file="observed_targets",outfn='parsesched.dat',outdir=None):
    """
        Update the local copy of the googledex with the last observed star time.
        updateLocalStarlist(time,googledex_file="googledex.dat", observed_file="observed_targets")

    opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """

    if not outdir :
        outdir = os.getcwd()
    
    obslog = ObservedLog.ObservedLog(os.path.join(outdir,observed_file))

    outfn = os.path.join(outdir,outfn)
    if os.path.exists(outfn):
        star_table = astropy.io.ascii.read(outfn)
    else:
        return obslog.names, None

    toofn = os.path.join(outdir,toofn)
    if os.path.exists(toofn):
        too_table = astropy.io.ascii.read(toofn)
    else:
        too_table = None
        
    for name in obslog.names:
        index = obslog.names.index(name)
        obstime = obslog.times[index]
        owner = obslog.owners[index]        
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

        selection = (star_table['name'] == name) & (star_table['sheetn'] == owner)
        if any(selection):
            if jd > star_table['lastobs'][selection]:
                apflog( "Updating local googledex star %s from time %.4f to %.4f" % (name, star_table['lastobs'][selection], jd),echo=True)
                star_table['lastobs'][selection] = jd
                star_table['nobs'][selection] += 1
        elif too_table is not None:
            selection = (too_table['name'] == name) & (too_table['sheetn'] == owner)
            if any(selection) and jd > too_table['lastobs'][selection]:
                apflog( "Updating ToO target %s from time %.4f to %.4f" % (name, too_table['lastobs'][selection], jd),echo=True)
                too_table['lastobs'][selection] = jd
                too_table['nobs'][selection] += 1

    astropy.io.ascii.write(star_table,outfn, format='ecsv', overwrite=True)
    astropy.io.ascii.write(too_table,toofn, format='ecsv', overwrite=True)
    star_table = vstack(too_table,star_table)
    
    return obslog.names, star_table

def updateSheetLastobs(filename, sheetns=["Bstar"],ctime=None,certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json',outfn='parsesched.dat',outdir=None):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
        updateSheetLastobs(filename, sheetn="The Googledex",time=None,certificate='UCSC Dynamic Scheduler-4f4f8d64827e.json')

        filename - where the observations are logged
        sheetns - list of sheets that will be updated
        ctime - current time as a time stamp 
        certificate - required for authentication
        outfn - the local copy of the star list
        outdir - directory of data files
        returns the number of cells updated
    """

    if not outdir :
        outdir = os.getcwd()
    
    obslog = ObservedLog.ObservedLog(os.path.join(outdir,observed_file))
    if len(obslog.names) == 0:
        return
    if ctime is None:
        ctime = datetime.utcfromtimestamp(int(time.time()))

    outfn = os.path.join(outdir,outfn)
    star_table = astropy.io.ascii.read(outfn)

    # OK the following code is going to rely on the fact that owner
    # is the sheet name
    # this will need to be updated when we start using coverid
    needed_sheetns = list( set(obslog.owners))
        
    nupdates = 0
    for sheetn in needed_sheetns:
        ws = getSpreadsheet(sheetn=sheetn,certificate=certificate)
        
        if ws:
            vals = ws.get_all_values()
        else:
            continue
        col = vals[0].index("lastobs") 
        nobscol = vals[0].index("Nobs")
        tempcol = vals[0].index("Template")

        wait_time = len(vals)
        time.sleep(wait_time)
    
        for i, v in enumerate(vals):
            # Did we observe this target tonight?
            local_name = parseStarname(v[0])
            if local_name in obslog.names:
                # We observed this target, so update the cell in the worksheet
                # update_cell(row, col, val) - col and row are 1 indexed
                nameidx = obslog.names.index(local_name)
                otime = obslog.times[nameidx]
                taketemp = obslog.temps[nameidx]
                curowner = obslog.owners[nameidx]
                try:
                    star_table_row = star_table[(star_table['name'] == local_name)&(star_table['sheetn'] == sheetn)]
                except:
                    star_table_row = None

                if  star_table_row is not None:
                    jd = float(star_table_row['lastobs'][0])
                else:
                    if isinstance(otime,float):
                        t = datetime.utcfromtimestamp(otime)
                    else:
                        hr, mn = otime
                        t = datetime(ctime.year, ctime.month, ctime.day, hr, mn)
                    jd = float(ephem.julian_date(t))
                try:
                    pastdate = float(v[col])
                    try:
                        n = int(v[nobscol])
                    except:
                        n = 0
                    if jd > pastdate and curowner == sheetn:
                        ws.update_cell(i+1, col+1, round(jd, 4) )
                        ws.update_cell(i+1, nobscol+1, n + 1 )
                        nupdates += 2
                        apflog( "Updated %s to %.4f and %d in %s" % (v[0],round(jd, 4),n+1,sheetn),echo=True)
                except:
                    print (v[0], v[col])
                    ws.update_cell(i+1, col+1, round(jd,4) )
                    ws.update_cell(i+1, nobscol+1, 1 )
                    nupdates += 2
                try:
                   have_temp = v[tempcol]
                   if taketemp == "Y" and have_temp == "N" and curowner == sheetn:
                       ws.update_cell(i+1, tempcol+1, "Y")
                       nupdates += 1
                       apflog( "Updated %s to having a template in %s" % (v[0],sheetn),echo=True)                       
                except:
                    apflog( "Error logging template obs for %s" % (v[0]),echo=True,level='error')

    return nupdates

