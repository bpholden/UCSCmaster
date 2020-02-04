from __future__ import print_function
import sys
import os
sys.path.append("../master")

import ParseUCOSched

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print ("needs a sheet list")
        sys.exit()
    sheetns = sys.argv[1].split(",")
    outdir = "."
    outfn = "googledex.dat"
    if os.path.exists(os.path.join(outdir,outfn)):
        os.unlink(os.path.join(outdir,outfn))


    config = dict()
    config['I2'] = 'Y'
    config['decker']='W'
    config['mode']=''
    config['obsblock']=''
    config['Bstar']='N'
    config['owner']='public'
    config['inst']='levy'
    config['raoff'] = None
    config['decoff'] = None 

        
    ParseUCOSched.parseUCOSched(sheetns=sheetns,outfn=outfn,outdir=outdir,config=config)

