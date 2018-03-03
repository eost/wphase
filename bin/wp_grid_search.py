#!/usr/bin/env python
# *-* coding: iso-8859-1 *-*

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

# Insert WPHASE BIN in sys.path (Useful to get local versions)
import os,sys
WPHOME = os.path.expandvars('$WPHASE_HOME')
WPBIN = os.path.join(WPHOME,'bin')
sys.path.insert(0,'./')
sys.path.insert(1,WPBIN)

# Avoid writing pyc files
sys.dont_write_bytecode = True

# GRID SEARCH FOR WPHASE INVERSION
from Arguments import *

# Import external modules
import shutil,time,getopt

# Import internal modules
from EQ import *
import utils


def addrefsol(cmtref,cmtfile):
    '''
    Adding reference moment tensor included in file cmtref to cmtfile
    '''
    cmtf = open(cmtref,'r')
    L=cmtf.readlines()
    cmtf.close()
    cmtf = open(cmtfile,'a')
    if len(L) < 13:
        raise EOFError('incomplete cmtfile: %s'%(cmtref))
    for l in L[7:]:
        cmtf.write(l)
    cmtf.close()
    # All done
    return;


def grid_search(eq,cmtref,ts_Nit,ts_dt,tsb,xy_Nit,xy_dx,xy_Nx,xy_Nopt,fastflag,flagts,flagxy,sdrM0={},dz=0.,
        minz=3.5,ts_ofile='grid_search_ts_out',xy_ofile='grid_search_xy_out',stdoutput='stdout',
        logfile='LOG/gs_o_wpinversion.log', comments = []):
    '''
    Grid search
    Args:
        * eq: eq object
        * cmtref: reference CMTSOLUTION file
        * ts_Nit: number of iteration for time-shift grid-search
        * ts_dt: initial sampling step for time-shift grid-search
        * tsb: bounds of time-shift grid search
        * xy_Nit: number of iteration of lat/lon grid-search
        * xy_dx: sampling step for lat/lon grid-search 
        * xy_Nopt: Number of optimum neighbor regions to re-sample
        * fastflag: perform time-shift grid-search ? (True or False)
        * flagxy: perform lat/lon grid-search? 
        * sdrM0: input dictionary for double-couple inversions
        * dz: sampling step for depth grid-search
        * minz: minimum depth for depth grid-search
        * ts_ofile: time-shift grid search output file
        * xy_ofile: lat/lon grid search output file
        * stdoutput: standard output
        * logfile: log filename
        * comments: comments to be added to the output ps file
    '''

    # Standard output
    if stdoutput == 'stdout':
        fid = sys.stdout
        flag = False
    else:
        fid = open(stdoutput,'a+')
        flag = True
    fid.write('CENTROID GRID SEARCH\n')

    # Setting parameters
    cmttmp = cmtref
    optpar = ' -log %s -osyndir gs_SYNTH -icmtf %s '%(logfile,cmtref)
    for o,a in sdrM0.items():
        if len(a):
            optpar += ' %s %s '%(o,a)
        else:
            optpar += ' %s '%(o)
    if not os.access('gs_SYNTH',os.F_OK):
        os.mkdir('gs_SYNTH')    

    # Prepare command line
    EXE = WPINV_XY        
    if flagts: # time-shift
        if len(tsb) == 2:
            ts1 = tsb[0]
            ts2 = tsb[1]
        else:        
            if eq.mag < 5.5:
                ts1 = 1.
                ts2 = eq.ts*3.    
                if ts2 > 100.:
                    ts2 = 100.    
            else:
                ts1 =  1. 
                if eq.mag <= 7.0:
                    ts2 = 30. 
                elif eq.mag <= 8.0:
                    ts2 = 100. 
                else: 
                    ts2 = 168.
        optpar += ' -ts %10.4f %10.4f %10.4f -ts_Nit %d -otsgsf %s'%(ts1,ts_dt,ts2,ts_Nit,ts_ofile)
        if not fastflag:
            optpar += ' -hdsafe '
    else:
        optpar += ' -nots '
    if flagxy:
        optpar += ' -xy_Nit %d -dx %.2f -Nx %d -Nopt %d -oxygsf %s'%(xy_Nit,xy_dx,xy_Nx,xy_Nopt,xy_ofile)
        if dz>0.:
            optpar += ' -dz %.2f -minz %.2f '%(dz,minz)
        wcmtfile = 'xy_WCMTSOLUTION'        
    else:
        optpar += ' -noxy '
        if flagts:
            wcmtfile = 'ts_WCMTSOLUTION'
        else:
            wcmtfile = 'gs_WCMTSOLUTION'
    if flag:
        fid.close()
        optpar += ' > %s '%stdoutput
    for c in comments:
        optpar += ' -comments "'+c+'"'

    # Run grid-search
    fid.write('Command_line:'+EXE+optpar+'\n')
    fid.flush()
    os.system(EXE+optpar)
    fid.flush()

    # Update eq after grid-search
    eq.rcmtfile(wcmtfile)
    out = utils.grep(r'^Wmag:',logfile)
    eq.mag = float(out[-1].split()[1]) ;

    # All done
    return;


def usage():
    print('usage: wp_grid_search [-s] [-t] [-p] [-i] ... [--help]')
    # All done
    return;


def disphelp():
    print('Centroid time-shift and centroid position grid search\n')
    usage()
    print('\nAll parameters are optional:')
    print('   -s, --hdsafe         Use a  time grid-search considering ts=fd')
    print('   -t, --onlyts         Centroid time-shift grid search only')
    print('   -p, --onlyxy         Centroid position grid search only')
    print('   -S, --npar           Do not use the parallelized grid-search and use ')
    print('                          the sequential version instead (parallelized version is used)')
    print('   -i, --imas \'file\'    Set i_master file (i_master)')
    print('   -n, --noref          Do not use the reference solution in cmtfile (ref. sol. used)')
    print('   --nont               Full moment tensor inversion (no null trace)')    
    print('   --dc                 Double-couple inversion')
    print('   --strike \'strike\'    Double-couple inversion with fixed strike')
    print('   --dip \'dip\'          Double-couple inversion with fixed dip')
    print('   --rake \'rake\'        Double-couple inversion with fixed rake')
    print('   --mom \'mom\'          Double-couple inversion with fixed scalar moment')
    print('\n   -h, --help           Display this help and exit')
    print('\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>')
    # All done
    return;


def main(argv):
    # Extract command line options
    try:
        opts, args = getopt.gnu_getopt(argv[1:],'stpSdi:nhz',["hdsafe","onlyts","onlyxy","npar",
                                      "imas=","strike=","dc","nont","dip=",
                                      "rake=","mom=","noref","xyz","old",
                                      "help"])
    except getopt.GetoptError as err:
        usage()
        raise

    # Parse command line options
    i_master = IMASTER
    fastflag = True
    flagts   = True
    flagxy   = True
    flagxyz  = False
    flagref  = True
    sdrM0    = {}
    for o, a in opts:
        if o == '-h' or o == '--help':
            disphelp()
            sys.exit(0)
        if o == '-s' or o == '--hdsafe':
            fastflag = False
        if o == '-t' or o == '--onlyts':
            if not flagts:
                usage()
                raise getopt.GetoptError('options -t and -p cannot be used simultaneously')
            flagxy = False
            flagts = True
        if o == '-p' or o == '--onlyxy':
            if not flagxy:
                usage()
                raise getopt.GetoptError('options -t and -p cannot be used simultaneously')                
            flagts   = False
            fastflag = False
            flagxy = True
        if o == '--dc':
            sdrM0['-dc']=''
        if o == '--nont':
            sdrM0['-nont']=''            
        if o == '--strike':
            sdrM0['-strike']=a
        if o == '--dip':
            sdrM0['-dip']=a
        if o == '--rake':
            sdrM0['-rake']=a
        if o == '--mom':
            sdrM0['-mom']=a
        if o == '-i' or o == '--imas':
            i_master = a
        if o == '-n' or o == '--noref':
            flagref = False
        if o == '-z' or o == '--xyz':
            flagxyz = True
        if o == '--old':
            WPINV_XY += ' -old'

    # Read i_master
    iconfig = utils.parseConfig(i_master)
    cmtref  = iconfig['CMTFILE']
    evname  = iconfig['EVNAME'].replace(' ','_').replace(',','')

    # Set comments in output ps file
    Median    = '-med '
    if 'P2P_SCREENING' in iconfig:
        p2p_items = iconfig['P2P_SCREENING'].split()
        if p2p_items[0] != 'YES':
            Median = ' '
        elif len(p2p_items)>1:
            for p in p2p_items[1:]:
                Median += p+' '
    ths = '-th 5.0 3.0 0.9'
    if 'RMS_SCREENING' in iconfig:
        ths = '-th '+iconfig['RMS_SCREENING']
    comments = [VERSION,'GF_PATH: '+GF_PATH,'Screening: '+Median+ths]

    # Read reference CMTFILE
    eq   = EarthQuake()
    eq.rcmtfile(cmtref)
    eq.title = evname.strip().replace(' ','_').replace(',','')
    cmtf = open(cmtref,'r')
    L=cmtf.readlines()
    cmtf.close()
    if len(L) < 13:
        print('*** WARNING : no reference solution in %s'%(cmtref))
        flagref = False

    # TS and/or LAT/LON Grid-search
    if (flagts or flagxy) and not flagxyz:
        grid_search(eq,cmtref,TS_NIT,TS_DT,TSBOUNDS,XY_NIT,XY_DX,XY_NX,XY_NOPT,fastflag,
                    flagts,flagxy,sdrM0,ts_ofile=TS_OFILE,xy_ofile=XY_OFILE,comments=comments)

    # TS and LAT/LON/DEP Grid-search
    if flagxyz:
        grid_search(eq,cmtref,TS_NIT,TS_DT,TSBOUNDS,XYZ_NIT,XYZ_DX,XYZ_NX,XYZ_NOPT,fastflag,
                    flagts,flagxyz,sdrM0,dz=DDEP,minz=MINDEP,ts_ofile=TS_OFILE,xy_ofile=XYZ_OFILE,
                    comments=comments)
        if flagxy:
            eq.wcmtfile('_tmp_CMTSOLUTION.xyz')
            if flagref:
                addrefsol(cmtref,'_tmp_CMTSOLUTION.xyz')
            grid_search(eq,'_tmp_CMTSOLUTION.xyz',TS_NIT,TS_DT,TSBOUNDS,XY_NIT,XY_DX,XY_NX,XY_NOPT,
                    0,0,1,sdrM0,ts_ofile=TS_OFILE,xy_ofile=XY_OFILE,comments=comments)
            utils.rm('_tmp_CMTSOLUTION.xyz')
    
    # Cleaning up
    if os.path.exists('_tmp_ts_table'):        
        utils.rm('_tmp_ts_table')
    if os.path.exists('_tmp_xy_table'):        
        utils.rm('_tmp_xy_table')


if __name__ == "__main__":
    main(sys.argv)
