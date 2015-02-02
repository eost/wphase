#!/usr/bin/env python
# -- coding: iso-8859-1 --

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

from sys        import stdin, stdout, stderr
from os.path    import exists, isdir
from shutil     import rmtree
from os         import mkdir, unlink, environ, chdir, rename
from glob       import glob
from subprocess import call


def RM(name):
    if exists(name):
        if isdir(name):
            rmtree(name)
        else:
            unlink(name)
    # All done
    return;


def mkdir_new(name):
    if exists(name):
        RM(name)
    mkdir(name)
    # All done
    return;


def get_sac_header(file, var):
    tmpfile = '_tmp_'
    oo=open(tmpfile,'w')
    cmd = 'saclst %s f %s'%(var,file)
    call(cmd, shell=True, stdin=stdin, stdout=oo, stderr=oo)
    oo.close()
    val = open(tmpfile,'r').readlines()[0].strip().split()[1]
    unlink(tmpfile) 
    # All done
    return val
    

def ch_sac_header(file, var, value):
    cmd = 'sac<<FIN\nrh %s\nch %s %s\nwh\nq\nFIN\n'%(file,var,value)
    call(cmd, shell=True, stdin=stdin, stdout=stdout, stderr=stderr)


def main():
    TOL       = 5.   # degrees (cmpaz tolerance)
    RDSEED    = environ['RDSEED']
    DATA      = 'DATA_org'
    LOG       = 'LOG'
    i_master  = 'i_master'

    if not exists(i_master):
        stderr.write('Error: i_master file not available.\n')
        exit(1)

    mkdir_new(DATA)
    mkdir_new(LOG)

    seeds = []
    lines = open(i_master).readlines()
    for line in lines:
        if line[0:5] == 'SEED:':
            seed = line.strip().split()[1]
            if not exists(seed):
                stderr.write('Seed file: %s not available.\n'%seed)
            else:
                seeds.append(seed)
    if len(seeds) < 1:
        stderr.write('No data extracted.\n')
        exit(1)

    oo = open(LOG+'/_log_rdseed','w')
    oo.write('\n')
    for seed in seeds:
        stdout.write('Extracting data from %s\n'%seed)
        cmd = RDSEED + ' -dp -z 3 -q ' + DATA + ' -f ' + seed
        call(cmd, shell=True, stdin=stdin, stdout=oo, stderr=oo)
    oo.close()

    # Patch
    chdir(DATA)
    Hsacs = glob('*.SAC')
    for sac in Hsacs:
        doit = False
        kcmpnm = get_sac_header(sac, 'kcmpnm')
        k1,k2,k3 = kcmpnm
        if (k1 == 'H' or k1 == 'B' or k1 == 'L') and (k3 == '1' or k3 == '2'):    
            cmpaz = float(get_sac_header(sac, 'cmpaz'))
            if abs(cmpaz - 90.) < TOL:
                newcmp = kcmpnm.replace(k3,'E'); doit = True
            if abs(cmpaz) < TOL or abs(360.-cmpaz) < TOL:
                newcmp = kcmpnm.replace(k3,'N'); doit = True
        if doit:
            kstnm   = get_sac_header(sac, 'kstnm')
            knetwk  = get_sac_header(sac, 'knetwk')
            locid   = get_sac_header(sac, 'khole')
            ch_sac_header(sac, 'kcmpnm', newcmp)
            newname = sac.replace('.'+kcmpnm+'.','.'+newcmp+'.')
            stdout.write('Renaming: %s\n -------> %s\n'%(sac, newname))
            rename(sac, newname)

            # PZ (only renaming)
            PZs = glob('SAC_PZs_%s_%s_%s_%s_*'%(knetwk,kstnm,kcmpnm,locid))
            for pz in PZs:
                newname = pz.replace('_'+kcmpnm+'_','_'+newcmp+'_')
                stdout.write('Renaming: %s\n -------> %s\n'%(pz, newname))
                rename(pz, newname)
    # All done
    return 0

if __name__=='__main__':
    main()
        
    

    


