#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
import time

sys.path.append("../master")
import ParseGoogledex 

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print "needs a filename and a list of sheet names"
        sys.exit()
    fn = sys.argv[1]
    if len(sys.argv) >= 3:
        sheetnl = sys.argv[2]
    else:
        sheetnl = "Bstars"

    sheetns = sheetnl.split(",")
    for sheetn in sheetns:
        n= ParseGoogledex.updateGoogledexLastobs(fn,sheetns=[sheetn])
        if n > 0:
            time.sleep(n)
