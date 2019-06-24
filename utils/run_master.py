#!/opt/kroot/bin/kpython
from __future__ import print_function
import os
import sys
import re
from time import sleep

import subprocess
import shlex
import ConfigParser

import ktl
import apflog
import numpy as np

def readem_or_weep(service,keyword,binary=False):

    try:
        value = ktl.read(service,keyword,binary=binary)
    except:
        apflog.apflog("cannot read %s.%s" % (service,keyword),echo=True,level='error')
        sys.exit("cannot read %s.%s" % (service,keyword))

    return value

def writeem(service,keyword,value,binary=False):

    try:
        ktl.write(service,keyword,value,binary=binary)
    except:
        raise
        apflog.apflog("cannot write %s.%s" % (service,keyword),echo=True,level='error')
        sys.exit("cannot write %s.%s" % (service,keyword))

    return value


def primary_run(schedule):
    match = re.search("\A\s*(\d+)\s*",schedule)
    if match:
        return match.group(1)

def remap_config(configlist):
    config = dict()
    for it in configlist:
        config[it[0]] = it[1]
    return config

def read_config(configfile,runstr):
    srun = primary_run(runstr)
    config = ConfigParser.ConfigParser()
    if not os.path.exists(configfile):
        apflog.apflog("configuration file %s does not exist"  % configfile)
        sys.exit("configuration file %s does not exist"  % configfile)
    config.read(configfile)
    programs = remap_config(config.items('programs'))
    program = programs[srun]
    program_config = remap_config(config.items(program))

    if program_config['obsnum'] == "ucsc":
        program_config['obsnum'] = findUCSCObsNum()
    if program_config['name'] == "ucb":
        program_config['name'] = findUCBObsNum()

    return program_config

def findUCSCObsNum():
    last = int(ktl.read('apftask','MASTER_LAST_OBS_UCSC',binary=True))
        
    last += 100 - (last % 100)

    if last % 10000 > 9700:
        last += 10000 - (last % 10000)

    return last


def findUCBObsNum(lastcode=None):
    apftask = ktl.Service('apftask')
    if lastcode == None:
        lastcode = apftask['MASTER_LAST_OBS_UCB'].read()
        if os.path.isfile('/data/apf/ucb-%s100.fits' % lastcode):
            apflog.apflog( "Existing files detected for run %s. Not incrementing night code." % lastcode,echo=True)
            return lastcode

    zloc = list(np.where(np.array(list(lastcode)) == 'z')[0])

    if 2 not in zloc:
        ncode = lastcode[0:2] + chr(ord(lastcode[2])+1)   # increment last
    if 2 in zloc and 1 not in zloc:
        ncode = lastcode[0] + chr(ord(lastcode[1])+1) + 'a'   # increment middle
    if 1 in zloc and 2 in zloc:
        ncode = chr(ord(lastcode[0])+1) + 'aa'   # increment first

    apftask['MASTER_LAST_OBS_UCB'].write(ncode)
    ncode = 'ucb-' + ncode
    return ncode


def ok_config(config):

    return True

def build_exec_str(config):
    stuff_to_run = None
    if os.path.exists(config['executable']):

        arglist = config['args'].split()
        actual_args = [ config[a] for a in arglist] 
        args = config['argstr'] % tuple(actual_args)
        stuff_to_run = []
        stuff_to_run.append(config['executable'])
        stuff_to_run = stuff_to_run + shlex.split(args)
    else:
        sys.exit("file %s does not exist"  % (config['executable']))
    return stuff_to_run    

def define_observer(config):

    obs = "Observer: %s \n" % config['observer'] 
    obs = obs + "Phone: %s \n" % (config['phone'])
    obs = obs + "Email: %s \n" % (config['email'])
    obs = obs + "Locat: %s \n" % (config['locat'])

    return obs

def config_kwds(config,test=False):
    obs = define_observer(config)
    ownr = readem_or_weep('apfschedule','OWNRNAME')
    
    if test:
        print(obs)
        print(config['observer'])
        print(primary_run(schedule))
        print(ownr)
        print(config['owner'])
        print(config['name'])
        print(config['obsnum'])
    else:
        writeem('checkapf','OBSLOCAT',obs)
        writeem('apfucam','OBSERVER',config['observer'])
        writeem('apfschedule','ACTIVE_RUN',primary_run(schedule))


        writeem('apfschedule','OWNRHINT',ownr)

        if config['name'] == 'ucsc':
            writeem('apfucam','outfile',config['name'])
            writeem('apfucam','obsnum',config['obsnum'])
            
        if config['owner']:
            writeem('apfschedule','ownrhint',config['owner'])

    return True

def modify_env(config):
    cenv = os.environ
    cenv['PATH'] = config['pathvar']
    cenv['PYTHONPATH'] = config['pythonpathvar']
    return cenv

def gen_int_files(cmd_str,cpath):
    cln_str = re.sub("-o \d+ ","",cmd_str)
    
    watch_str = re.sub("Focus","Watching",cln_str)
    cal_str = re.sub("Focus","Cal-Pre",cln_str)
    
    fp = open(os.path.join(cpath,"watch.out"),"w+")
    fp.write(watch_str)
    fp.close()
    fp = open(os.path.join(cpath,"cal.out"),"w+")
    fp.write(cal_str)
    fp.close()

    return
           
if __name__ == "__main__":

    if len(sys.argv) > 1:
        test = True
    else:
        test = False
    
    schedule = readem_or_weep('apfschedule','SCHEDULED_RUNS')

    userkind = readem_or_weep('checkapf','USERKIND',binary=True)
    if userkind != 3 and test is False:
        sys.exit("checkapf not in robotic mode")


    master_status = readem_or_weep('apftask','master_status',binary=True)
    if master_status < 3 and test is False:
        sys.exit("master has been started")

    cpath = os.path.dirname(os.path.abspath(__file__))
    configfile = "master.config"
    config = read_config(os.path.join(cpath,configfile),schedule)

    if ok_config(config) and config_kwds(config):

        stuff_to_run = build_exec_str(config)
        env = modify_env(config)
        print(" ".join(stuff_to_run))
        gen_int_files(" ".join(stuff_to_run),cpath)
        if os.getuid() == int(config['user']) and test is False:
            p=subprocess.Popen(stuff_to_run,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
            sleep(10)
            if p.poll():
                print("Master failed for some reason or another.")
                (out,err) = p.communicate()
                print(err)
                print(out)
        elif test:
            pass
        else:
            print("Master cannot be run by this account")

        
