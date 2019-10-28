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

    
