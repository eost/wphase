'''
My Utils

Written by Z. Duputel, September 2013
'''

# Externals
from os.path import exists,isdir
from numpy   import ndarray
import shutil
import os
import re


def grep(chaine, file):
    """
    Returns lines in file matching chaine
    """
    out = [];
    template = re.compile(chaine)
    fid  = open(file, 'r')
    for line in fid:
        if template.match(line):
            out.append(line)
    fid.close()
    
    # All done
    return(out)


def rm(ifiles):
    '''
    Remove file(s) or directory(ies)
    Args:
         ifiles:  file, directory or list of files, directories
    '''

    if type(ifiles)!=list and type(ifiles)!=ndarray:
        ifiles = [ifiles]
    for ifile in ifiles:
        if exists(ifile):
            if isdir(ifile):
                shutil.rmtree(ifile)
            else:
                os.remove(ifile)
                
    # All done
    return


def mkdir(idir,pflag=False):
    '''
    Create directory
    Args:
         idir:   directory to be created
         pflag: (optional, default=False)
               False: will     remove directory if exists
               True:  will NOT remove directory if exists
    '''        

    if pflag and not exists(idir):
        os.mkdir(idir)
    else:
        rm(idir)
        os.mkdir(idir)
        
    # All done
    return


class parseConfigError(Exception):
    """
    Raised if the config file is incorrect
    """
    pass


def parseConfig(cfg_file):
    '''
    Parse my config files and returns a config dictionary
    Args:
         cfg_file: configuration filename
    '''

    # Check if cfg_file exists
    assert exists(cfg_file), 'Cannot read %s: No such file'%(cfg_file)

    # Fill the config dictionary
    config = {}
    try:
        config_lines = open(cfg_file, 'r').readlines()
        for line in config_lines:
            if line.find('#')==0:
                continue
            if line.rstrip():
                key,value = line.strip().split(':')
                key   = key.strip()
                value = value.strip()
                if key in config:
                    config[key].append(value)
                else:
                    config[key]=value
    except:
        raise parseConfigError('Incorrect format in %s!\n'%cfg_file)

    # All done
    return config
