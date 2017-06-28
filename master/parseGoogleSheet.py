import gspread
import json
from oauth2client.client import SignedJwtAssertionCredentials
import numpy as np
import os
import pickle
import sys
import time

def get_spreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """ Get the spreadsheet from google

    worksheet = get_spreadsheet(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')
    worksheet - the worksheet object returned by the gspread module

    sheetn - name of the google sheet, defaults to "The Googledex"
    certificate - certificate used to control access to the google sheet
    
    this downloads the googledex from the Google Drive
    the certificate must be available
    these certificates are generated through the Google Developer Interface
    the developer must select the correct API for access

    the certificate has an email associated with it, that email must
    have the document shared with it to allow access 
    """
    certificate_path = os.path.dirname(__file__)
    
    json_key = json.load(open(os.path.join(certificate_path, certificate)))
    scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gs = gspread.authorize(credentials)

    apflog("Successfully logged in.", echo=True)
    spreadsheet = gs.open(sheetn)
    apflog("Loaded Main %s" % (sheetn),echo=True)
    worksheet = spreadsheet.sheet1
    apflog("Got spreadsheet", echo=True)

    return worksheet

def parseGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat",outdir=None,config={'I2': 'Y', 'decker': 'W' }):
    """ parseGoogledex parses the google sheet and returns the output as a tuple
    This routine downloads the data if needed and saves the output to a file. If the file exists, it just reads in the file.
    
    names, star_table, do_flag, stars = parseGoogledex(sheetn="The Googledex",certificate='UCSC Dynamic Scheduler-5b98d1283a95.json',outfn="googledex.dat")
    names - a list of stars in the starlist
    star_table - a numpy array
    flags - a dictionary of items on whether or not do="y" needs to be set for scriptobs 
    stars - a list of pyEphem objects 

    """
    # Downloading all the values is going slowly.
    # Try to only have to load this once a day
    if not outdir :
        outdir = os.getcwd()
    try:
        f = open(os.path.join(outdir,outfn),'r')
    except IOError:
        apflog( "Starting Googledex parse",echo=True)
        worksheet = get_spreadsheet(sheetn=sheetn,certificate=certificate)
        full_codex = worksheet.get_all_values()
        #time = (datetime.now() - start).total_seconds()
        #print "Loaded Values. Took {0:f} seconds.".format(time)
        f = open(os.path.join(outdir,"googledex.dat"),'w')
        pickle.dump(full_codex, f)
        f.close()
    else:
        full_codex = pickle.load(f)
        f.close()
        
    col_names = full_codex[0]
    codex = full_codex[1:]

    # These are the columns we need for scheduling
    req_cols = ["Star Name", "RA hr", "RA min", "RA sec", \
                "Dec deg", "Dec min", "Dec sec", "pmRA", "pmDEC", "Vmag", \
                "APFpri", "APFcad", "APFnshots", "lastobs", \
                "B-V", "SpecCounts", "ExpMeter", "Close Companion", \
                "APF decker","I2","owner"
                ]
    didx = findColumns(col_names,req_cols)
    
    names = []
    star_table = []
    flags = { "do" : [], "decker" : [], "I2" : [] }
    stars = []
    # Build the star table to return to 
    for ls in codex:
        if ls[0] == '':
            continue
        try:
            apfpri = float(ls[didx["APFpri"]])
        except:
            apfpri = 0.0
        if apfpri < 0.5: continue
        row = []
        # Get the star name
        names.append(ls[didx["Star Name"]])
        # Get the RA
        raval = getRARad(ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]])
        if raval:
            row.append(raval)
        else:
            row.append(-1.)
        # Get the DEC
        decval = getDECRad(ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]])
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
        # For now use the old 1e9 count value
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

        check = checkflag("Close Companion",didx,ls,"\A(n|N)","Y")
        if check == "N" or check == "n":
            flags['do'].append("")
        else:
            flags['do'].append(check)
            
        flags['decker'].append(checkflag("APF decker",didx,ls,"\A(W|N|T|S|O|K|L|M|B)",config["decker"]))
        flags['I2'].append(checkflag("I2",didx,ls,"\A(n|N)",config["I2"]))
                                         
            
        star_table.append(row)
        star = ephem.FixedBody()
        star.name = ls[0]
        star._ra = ephem.hours(":".join([ls[didx["RA hr"]], ls[didx["RA min"]], ls[didx["RA sec"]]]))
        star._dec = ephem.degrees(":".join([ls[didx["Dec deg"]], ls[didx["Dec min"]], ls[didx["Dec sec"]]]))
        stars.append(star)

    return (names, np.array(star_table), flags, stars)

def update_googledex_lastobs(filename, sheetn="The Googledex",ctime=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json'):
    """
        Update the online googledex lastobs column assuming things in filename have been observed.
        update_googledex_lastobs(filename, sheetn="The Googledex",time=None,certificate='UCSC Dynamic Scheduler-5b98d1283a95.json')

        filename - where the observations are logged
    """
    names, times = getObserved(filename)
    if len(names) == 0:
        return
    if ctime is None:
        ctime = datetime.utcfromtimestamp(int(time.time()))
    

    ws = get_spreadsheet(sheetn=sheetn,certificate=certificate)
    vals = ws.get_all_values()

    col = vals[0].index("lastobs") 
    
    for i, v in enumerate(vals):
        # Did we observe this target tonight?
        if v[0] in names:
            # We observed this target, so update the cell in the worksheet
            # update_cell(row, col, val) - col and row are 1 indexed
            otime = times[names.index(v[0])]
            if isinstance(otime,float):
                t = datetime.utcfromtimestamp(otime)
            else:
                hr, min = otime
                t = datetime(ctime.year, ctime.month, ctime.day, hr, min)
            jd = float(ephem.julian_date(t))
            try:
                pastdate = float(v[col])
                if jd > pastdate:
                    ws.update_cell(i+1, col+1, round(jd, 2) )
            except:
                print v[0], v[col]
                
    apflog( "Updated Googledex",echo=True)
    return
