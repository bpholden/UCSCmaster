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
        program_config['obsnum'] = finducscObsNum()

    return program_config

def finducscObsNum():
    last = int(ktl.read('apftask','MASTER_LAST_OBS_UCSC',binary=True))
        
    last += 100 - (last % 100)

    if last % 10000 > 9700:
        last += 10000 - (last % 10000)

    

    return last

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

def config_kwds(config):
    obs = define_observer(config)

    writeem('checkapf','OBSLOCAT',obs)
    writeem('apfucam','OBSERVER',config['observer'])
    writeem('apfschedule','ACTIVE_RUN',primary_run(schedule))

    ownr = readem_or_weep('apfschedule','OWNRNAME')
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

        
if __name__ == "__main__":

    schedule = readem_or_weep('apfschedule','SCHEDULED_RUNS')

    userkind = readem_or_weep('checkapf','USERKIND',binary=True)
    if userkind != 3:
        sys.exit("checkapf not in robotic mode")


    master_status = readem_or_weep('apftask','master_status',binary=True)
    if master_status < 3:
        sys.exit("master has been started")

    cpath = os.path.dirname(os.path.abspath(__file__))
    configfile = "master.config"
    config = read_config(os.path.join(cpath,configfile),schedule)

    if ok_config(config) and config_kwds(config):

        stuff_to_run = build_exec_str(config)
        env = modify_env(config)
        print(" ".join(stuff_to_run))
        if os.getuid() == int(config['user']):
            p=subprocess.Popen(stuff_to_run,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
            sleep(10)
            if p.poll():
                print("Master failed for some reason or another.")
                (out,err) = p.communicate()
                print(err)
                print(out)
        else:
            print("Master cannot be run by this account")

        
