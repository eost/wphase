'''
My Utils

Written by Z. Duputel, September 2013
'''

############################################################################
#
#                  W phase source inversion package                 
#                               -------------
#
#        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#                      
# (c) California Institute of Technology and Universite de Strasbourg / CNRS 
#                                  April 2013
#
#    Neither the name of the California Institute of Technology (Caltech) 
#    nor the names of its contributors may be used to endorse or promote 
#    products derived from this software without specific prior written 
#    permission
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

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
                    config[key] = [config[key]]
                    config[key].append(value)
                else:
                    config[key]=value
    except:
        raise parseConfigError('Incorrect format in %s!\n'%cfg_file)

    # All done
    return config
