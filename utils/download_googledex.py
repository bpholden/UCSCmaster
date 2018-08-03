#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
sys.path.append("../master")
import UCSCScheduler_V2 as ds
import os

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print "needs a sheet list"
        sys.exit()
    sheetns = sys.argv[1].split(",")
    if os.path.exists("googledex.dat"):
        os.unlink("googledex.dat")
    ds.parseGoogledex(sheetns=sheetns)

