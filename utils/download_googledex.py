#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
sys.path.append("../master")
import UCSCScheduler_V2 as ds

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print "needs a sheet list"
        sys.exit()
    sheetns = sys.argv[1]
    ds.parseGoogledex(sheetns=sheetns)

