#!/usr/bin/env  /opt/kroot/bin/kpython

import sys
import time
import datetime

sys.path.append("../master")
import ParseGoogledex 

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("needs a log filename and googledex filename")
        sys.exit()
        obfn = sys.argv[1]
        gdfn = sys.argv[2]

    dt = datetime.datetime.now()
    ParseGoogledex.update_local_googledex(dt,googledex_file=gdfn,observed_file=obfn)
