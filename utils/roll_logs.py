#!/usr/bin/env  /opt/kroot/bin/kpython

import os
import sys
sys.path.append("../master")
import apflog

if __name__ == "__main__":
    apflog.logpush(os.path.join(os.getcwd(),"observed_targets"))    
    apflog.logpush(os.path.join(os.getcwd(),"robot.log"))    
