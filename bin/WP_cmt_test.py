#!/usr/bin/env python

from sys        import argv, stdout, stderr, exit
from os         import environ, mkdir, unlink
from shutil     import rmtree
from os.path    import exists, isdir
from subprocess import call

import utils

def RM(name):
    if exists(name):
        if isdir(name):
            rmtree(name)
        else:
            unlink(name)
    return

def mkdir_new(name):
    if exists(name):
        RM(name)
    mkdir(name)
    return

def main():
    BIN     = environ['WPHASE_HOME']+'/bin'
    EXTRACT = BIN + '/extract.csh'
    PREPARE = BIN + '/prepare_wp.csh'
    WPINVER = BIN + '/wpinversion'
    TRACES  = BIN + '/traces.py'

    if len(argv) > 1:
        CMTFILE = argv[1]
    elif not exists('i_master'):
        stderr.write('Error: file  i_master not available\n')
        exit(1)
    else: CMTFILE = utils.parseConfig('i_master')['CMTFILE']
    stdout.write('Using cmtfile = %s\n'%CMTFILE)
    mkdir_new('SYNTH')

    call(EXTRACT, shell=True, stdout=stdout)
    cmd = PREPARE + ' -a'
    call(cmd, shell=True, stdout=stdout)
    cmd = WPINVER + ' -osyndir SYNTH -nops -ocmtf /dev/null'
    oo  = open('/dev/null', 'r')
    call(cmd, shell=True, stdout=oo)
    oo.close()

    pdf = 'wp_pages.pdf'
    if exists(pdf): RM(pdf)
    cmd = TRACES + ' --icmtf ' + CMTFILE + ' --osydir SYNTH'
    call(cmd, shell=True, stdout=stdout)
    stdout.write('---> ' + pdf + '\n')

main()
