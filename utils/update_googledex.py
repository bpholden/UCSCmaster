#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
sys.path.append("../master")
import ParseGoogledex 

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print "needs a filename"
        sys.exit()
    fn = sys.argv[1]
    if len(sys.argv) >= 3:
        sheetnl = sys.argv[2]
    else:
        sheetnl = "The Googledex"

    sheetns = sheetnl.split(",")
    ParseGoogledex.update_googledex_lastobs(fn,sheetns=sheetns)
