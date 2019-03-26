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
from oauth2client.client import SignedJwtAssertionCredentials

import ObservedLog
import Coords
from SchedulerConsts import MIN_TOTOBS


try:
    from apflog import *
except:
    from fake_apflog import *

def checkflag(key,didx,line,regexp,default):
    try:
        match = re.search(regexp,line[didx[key]])
        if match:
            return match.group(1)
        else:
            return default
    except:
        return default


def parse_starname(starname):

    ostarname = starname
    m= re.search("HD\s+\d+",starname)
    if m:
        ostarname = re.sub("HD\s+","HD",starname)
    m = re.search("\s+",ostarname)
    while m:
        ostarname = re.sub("\s+","_",ostarname)
        m = re.search("\s+",ostarname)
        
    return ostarname

def parseGoogledex(sheetns=["Bstars"],certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat",outdir=None,config={'I2': 'Y', 'decker': 'W', 'owner' : '' }):
    """ parseGoogledex parses google sheets and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat")
    names - a list of stars in the starlist
    star_table - a numpy array
    flags - a dictionary of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    apflog( "Starting Googledex parse",echo=True)    
    if not outdir :
        outdir = os.getcwd()
    if os.path.exists(os.path.join(outdir,outfn)):
        try:
            f = open(os.path.join(outdir,outfn),'rb')
            full_codex = pickle.load(f)
            f.close()
        except:
            full_codex = make_local_copy(sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))
    else:
        full_codex = make_local_copy(sheetns=sheetns,certificate=certificate,outfn=os.path.join(outdir,outfn))

    col_names = full_codex[0]
    codex = full_codex[1:]

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", "APFmin", "APFmax", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2", "owner", "uth","utm","duration", "Template",
                "Nobs", "Total Obs"
                ]
    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    flags = { "do" : [], "decker" : [], "I2" : [], "owner" : [], "template" : [] }
    stars = []
    # Build the star table to return to 
    for ls in codex:
        if ls[0] == '':
            continue
        try:
            apfpri = float(ls[didx["APFpri"]])
        except:
            apfpri = 0.0
        try:
            nobs = int(ls[didx["Nobs"]])
        except:
            nobs = 0
        try:
            totobs = int(ls[didx["Total Obs"]])
        except :
            totobs = -1
        if totobs > 0 and nobs >= totobs: continue
        if apfpri < 0.5: continue
        row = []
        # Get the star name
        names.append(parse_starname(ls[didx["Star Name"]]))
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

        for coln in ("pmRA", "pmDEC", "Vmag"):
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln == "Vmag":
                    row.append(15.0)
                else:
                    row.append(0.0)
        # For now use the old 1e9 count value - these get recalculated 
        row.append(1200.0)
        row.append(1.e9)
        for coln in ["APFpri", "APFcad","APFnshots"] :
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln in ("APFpri","APFnshots"):
                    row.append(0)
                else:
                    row.append(1000.0)
                    
        
        for coln in ["lastobs", "B-V", "APF Desired Precision" ]:
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                if coln in ("lastobs", "B-V"):
                    row.append(0.0)
                else:
                    row.append(1000.0)

        for coln in ["uth", "utm"]:
            try:
                row.append(int(ls[didx[coln]]))
            except ValueError:
                row.append(0)
            except KeyError:
                row.append(0)
                
        for coln in ["duration"]:
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                row.append(0)
            except KeyError:
                row.append(0)
                
        for coln in ["APFmin"]:
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                row.append(MIN_TOTOBS)
            except KeyError:
                row.append(MIN_TOTOBS)
                
        for coln in ["APFmax"]:
            try:
                row.append(float(ls[didx[coln]]))
            except ValueError:
                row.append(0)
            except KeyError:
                row.append(0)

        for coln in ["Nobs"]:
            try:
                row.append(int(ls[didx[coln]]))
            except ValueError:
                row.append(0)
            except KeyError:
                row.append(0)
                
        for coln in ["Total Obs"]:
            try:
                row.append(int(ls[didx[coln]]))
            except ValueError:
                row.append(0)
            except KeyError:
                row.append(0)
                
                    
        check = checkflag("Close Companion",didx,ls,"\A(n|N)","Y")
        if check == "N" or check == "n":
            flags['do'].append("")
        else:
            flags['do'].append(check)
            
        flags['decker'].append(checkflag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        flags['I2'].append(checkflag("I2",didx,ls,"\A(n|N)",config["I2"]))
        flags['template'].append(checkflag("Template",didx,ls,"\A(n|N)",'Y'))

        flags['owner'].append(checkflag("owner",didx,ls,"\A(\w?\.?\w+)",config["owner"]))

            
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]]))
        star._dec = ephem.degrees(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]]))
        stars.append(star)

    return (names, np.array(star_table), flags, stars)

def update_googledex_lastobs(filename, sheetns=["2018B"],ctime=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
        update_googledex_lastobs(filename, sheetn="The Googledex",time=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')

        filename - where the observations are logged
    """
    names, times = ObservedLog.getObserved(filename)
    if len(names) == 0:
        return
    if ctime is None:
        ctime = datetime.utcfromtimestamp(int(time.time()))
    

    for sheetn in sheetns:
        ws = get_spreadsheet(sheetn=sheetn,certificate=certificate)
        if ws:
            vals = ws.get_all_values()
        else:
            next
        col = vals[0].index("lastobs") 
        nobscol = vals[0].index("Nobs")
    
        for i, v in enumerate(vals):
            # Did we observe this target tonight?
            if v[0] in names:
                # We observed this target, so update the cell in the worksheet
                # update_cell(row, col, val) - col and row are 1 indexed
                otime = times[names.index(v[0])]
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
                    if jd > pastdate:
                        ws.update_cell(i+1, col+1, round(jd, 2) )
                        ws.update_cell(i+1, nobscol+1, n + 1 )

                except:
                    print (v[0], v[col])
                    ws.update_cell(i+1, col+1, round(jd,2) )
                
            apflog( "Updated %s" % (sheetn),echo=True)

    return

def update_local_googledex(intime,googledex_file="googledex.dat", observed_file="observed_targets"):
    """
        Update the local copy of the googledex with the last observed star time.
        update_local_googledex(time,googledex_file="googledex.dat", observed_file="observed_targets")

        opens googledex_file and inputs date of last observation from observed_file
        in principle can use timestamps as well as scriptobs uth and utm values
    """
    names, times = ObservedLog.getObserved(observed_file)

    try:
        g = open(googledex_file, 'rb')
        full_codex = pickle.load(g)
        g.close()
    except IOError:
        apflog("googledex file did not exist, so can't be updated",echo=True)
        return names,times
    except EOFError:
        apflog("googledex file corrupt, so can't be updated",echo=True)
        return names,times


    codex_cols = full_codex[0]

    starNameIdx = codex_cols.index("Star Name")
    lastObsIdx = codex_cols.index("lastobs")
    try:
        nObsIdx = codex_cols.index("nObsIdx")
    except:
        nObsIdx = -1
    
    for i in range(1, len(full_codex)):
        row = full_codex[i]
        if row[starNameIdx] in names:
            # We have observed this star, so lets update the last obs field
            obstime = times[names.index(row[starNameIdx])]
            if isinstance(obstime,float):
                t = datetime.utcfromtimestamp(obstime)
            else:
                hr, min = obstime
                if type(intime) != datetime:
                    ctime = datetime.now()
                    td = timedelta(0,3600.*7)
                    intime = ctime + td
                t = datetime(intime.year, intime.month, intime.day, hr, min)

            # This keeps the JD precision to one decimal point. There is no real reason for this other than
            # the googledex currently only stores the lastObs field to one decimal precision. Legacy styles FTW.
            jd = round(float(ephem.julian_date(t)), 2) 
            apflog( "Updating local googledex star %s from time %s to %s" % (row[starNameIdx], row[lastObsIdx], str(jd)),echo=True)
            row[lastObsIdx] = str(jd)
            if nObsIdx > 0:
                row[nObsIdx] = row[nObsIdx] + 1
            full_codex[i] = row

    with open(googledex_file, 'wb') as f:
        pickle.dump(full_codex, f)
    f.close()
    
    return names, times

def make_local_copy(sheetns=["The Googledex"],certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="./googledex.dat"):
    full_codex = []
    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", "APFmin", "APFmax", \
                "B-V", "APF Desired Precision", "Close Companion", \
                "APF decker","I2", "owner", "uth","utm","duration", "Template",
                "Nobs", "Total Obs"
                ]
    full_codex.append(req_cols)
        
    for sheetn in sheetns:
        worksheet = get_spreadsheet(sheetn=sheetn,certificate=certificate)
        if worksheet:
            cur_codex = worksheet.get_all_values()
            didx = findColumns(cur_codex[0],req_cols)
            
            for row in cur_codex[1:]:
                nrow = []
                for c in req_cols:
                    nrow.append(row[didx[c]])
                full_codex.append(nrow)
        

    f = open(outfn,'wb')
    pickle.dump(full_codex, f)
    f.close()
    return full_codex
    
def get_spreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """ Get the spreadsheet from google

    worksheet = get_spreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')
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
    scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)

    try:
        gs = gspread.authorize(credentials)
        apflog("Successfully logged in.", echo=True)
    except:
        apflog("Cannot logg in.", echo=True,level='error')
        return None
    apflog("Attempting to Open %s" % (sheetn),echo=True)
    try:
        spreadsheet = gs.open(sheetn)
        apflog("Loaded Main %s" % (sheetn),echo=True)
        worksheet = spreadsheet.sheet1
        apflog("Got spreadsheet", echo=True)
    except:
        apflog("Cannot Read %s"  % (sheetn), echo=True, level='error')
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

