#!/opt/kroot/bin/kpython

import ktl
import APF
import subprocess
import shlex
from time import sleep
import ConfigParser
import os.path
import os
import sys
import re

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
    run = primary_run(runstr)
    config = ConfigParser.ConfigParser()
    if not os.path.exists(configfile):
        sys.exit("configuration file %s does not exist"  % configfile)
    config.read(configfile)
    programs = remap_config(config.items('programs'))
    program = programs[run]
    program_config = remap_config(config.items(program))

    if program == "ucsc":
        program_config['obsnum'] = finducscObsNum()

    return program_config

def finducscObsNum():

    last = int(ktl.read('apftask','MASTER_LAST_OBS_UCSC'))
        
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
    try:
        ktl.write('checkapf','OBSLOCAT',obs)
    except:
        raise
        sys.exit('Cannot communicate with checkapf service')
    try:
        ktl.write('apfucam','OBSERVER',config['observer'])
    except:
        raise
        sys.exit('Cannot communicate with checkapf service')
    try:
        ktl.write('apfschedule','ACTIVE_RUN',primary_run(schedule))
    except:
        raise
        sys.exit('Cannot communicate with apfschedule service')

    if config['name'] == 'ucsc':
        try: 
            ktl.write('apfucam','outfile',config['name'])
            ktl.write('apfucam','obsnum',config['obsnum'])
        except:
            raise
            sys.exit('Cannot communicate with apfucam service')

    if config['owner']:
        try: 
            ktl.write('apfschedule','ownrhint',config['owner'])
        except:
            raise
            sys.exit('Cannot communicate with apfschedule service')

    return True

def modify_env(config):
    cenv = os.environ
    cenv['PATH'] = config['pathvar']
    cenv['PYTHONPATH'] = config['pythonpathvar']
    return cenv
            
if __name__ == "__main__":

    try:
        schedule = ktl.read('apfschedule','SCHEDULED_RUNS')
    except:
        sys.exit("cannot read apfschedule.SCHEDULED_RUNS")
        
    cpath = os.path.dirname(os.path.abspath(__file__))
    configfile = "master.config"
    masterstatus=ktl.read('apftask','MASTER_PID',binary=True)
    config = read_config(os.path.join(cpath,configfile),schedule)

    if masterstatus < 0 and ok_config(config) and config_kwds(config):

        stuff_to_run = build_exec_str(config)
        env = modify_env(config)
        print " ".join(stuff_to_run)
        if os.getuid() == int(config['user']):
            p=subprocess.Popen(stuff_to_run,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
            sleep(10)
            if p.poll():
                print "Master failed for some reason or another."
                (out,err) = p.communicate()
                print err
                print out
        else:
            print "Master cannot be run by this account"
    else:
        print "Master script is already running"
