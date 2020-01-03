from __future__ import print_function
import sys
import os
sys.path.append("../master")

import APFTask
import ParseGoogledex

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print ("needs a sheet list")
        sys.exit()
    sheetns = sys.argv[1].split(",")
    if os.path.exists("googledex.dat"):
        os.unlink("googledex.dat")

    ParseGoogledex.parseGoogledex(sheetns=sheetns)

