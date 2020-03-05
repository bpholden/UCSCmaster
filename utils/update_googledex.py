#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
import time

sys.path.append("../master")
import ParseUCOSched
import ObservedLog

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print "needs a observed log and a local copy of the star table"
        sys.exit()
    fn = sys.argv[1]
    outfn = sys.argv[2]

    obslog = Observedlog.Observedlog(fn)

    if len(obslog.names) > 0:
        if obslog.sheetns[0] is None:
            sheetns = set(ol.owners)
        else:
            sheetns = set(ol.sheetns)
    for sheetn in sheetns:
        n= ParseUCOSched.updateSheetLastobs(fn,sheetns=[sheetn],outfn=outfn)
        if n > 0:
            time.sleep(n)
