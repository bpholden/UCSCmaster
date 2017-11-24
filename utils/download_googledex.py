#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
import json
import gspread
import pickle
from oauth2client.client import SignedJwtAssertionCredentials

json_key = json.load(open('../master/UCSC Dynamic Scheduler-5b98d1283a95.json'))
scope = ['https://www.googleapis.com/auth/plus.login https://www.googleapis.com/auth/plus.me https://spreadsheets.google.com/feeds']

credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
gc = gspread.authorize(credentials)

sheetn ="The Googledex"

if len(sys.argv) > 1:
    sheetn = sys.argv[1]

spreadsheet = gc.open(sheetn)
#spreadsheet = gc.open_by_key("1VoNDlmtVSnqJqzWbbCzYGc7XAvR9ADNjNIDSIBZE-jE")
# 1VoNDlmtVSnqJqzWbbCzYGc7XAvR9ADNjNIDSIBZE-jE
#spreadsheet = gc.open_by_url("https://docs.google.com/spreadsheets/d/1VoNDlmtVSnqJqzWbbCzYGc7XAvR9ADNjNIDSIBZE-jE/")
#https://docs.google.com/spreadsheets/d/1VoNDlmtVSnqJqzWbbCzYGc7XAvR9ADNjNIDSIBZE-jE/edit?usp=sharing
#allspreads = gc.openall()
#print allspreads
print "got spreadsheet"
worksheet = spreadsheet.sheet1
full_codex = worksheet.get_all_values()
print "got all values from worksheet"
f = open("./googledex.dat",'w')
pickle.dump(full_codex, f)
print "dumped a pickled file"
f.close()
