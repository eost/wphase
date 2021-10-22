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

# Append WPHASE BIN to pythonpath
import os,sys
WPHOME = os.path.expandvars('$WPHASE_HOME')
WPBIN = os.path.join(WPHOME,'bin')
sys.path.insert(0,'./')
sys.path.insert(1,WPBIN)

# Avoid writing pyc files
sys.dont_write_bytecode = True

# CONCATENATED W PHASE TRACES
from Arguments import *

#Customizing matplotlib
import matplotlib as mpl
mpl.use('PDF')


# Import external modules
import re
import getopt as go
import numpy as np
import matplotlib.pyplot as plt


# Import internal modules
import sacpy


# Internal functions
def usage(cmd):
    print('usage: %s [chan1 chan2 (default: LHZ LHN LHE LH1 LH2)] [option] (for help see %s -h)'%(cmd,cmd))
    # All done
    return;


def disphelp(cmd):
    print('Make CWP traces\n')
    usage(cmd)
    print('\nAll parameters are optional:')
    print('\n   -i, --ifort15             input fort.15 file (e.g. fort.15, ts_fort.15, xy_fort.15)')
    print('\n   -n, --noref               no reference solution')
    print('\n   -h, --help                display this help and exit')
    print('\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>')
    # All done
    return;


def main(argv):
    # Input parameters
    try:
        opts, args = go.gnu_getopt(argv[1:],'i:nh',["ifort15=","noref","help"])
    except getopt.GetoptError as err:
        raise
    o_wpfile = O_WPINVERSION
    predfile = ''
    isref    = 1
    CHAN     = ['LHZ', 'LHN', 'LHE', 'LH1', 'LH2']
    for o, a in opts:
        if o == '-h' or o == '--help':
            disphelp(sys.argv[0])
            sys.exit(0)
        if o == '-i' or o == '--ifort15':
            predfile = a
        if o == '-n' or o == '--noref':
            isref = 0    
    if len(args):
        CHAN  = args
    if not os.path.exists(o_wpfile):
        sys.stderr.write('Error: file %s not available\n'%(o_wpfile))
    if len(predfile) and not os.path.exists(predfile):
        raise IOError('No fort.15 file named %s\n'%(predfile))
    if not len(predfile):
        predfile = 'xy_fort.15'
        if not os.path.exists(predfile):
            predfile = 'ts_fort.15'
            if not os.path.exists(predfile):
                predfile = 'fort.15'
                if not os.path.exists(predfile):
                    raise IOError('No fort.15 file found\n')
    sys.stdout.write('Input fort.15 file: %s\n'%predfile)
    count = 0
    sys.stdout.write('Input channels are: ')
    for chan in CHAN:
        if not os.path.exists('%s_%s'%(predfile,chan)):
            continue 
        else:
            count += 1
        sys.stdout.write('%5s'%chan)
    if not count:
        ErrMsg = 'Error: No fort.15_ file for'
        for chan in CHAN:
            ErrMsg += '%5s'%(chan)
        raise IOError(ErrMsg)
    sys.stdout.write('\nRead %s ...\n%s pages:\n'%(o_wpfile,count))

    # Main loop
    sac = sacpy.Sac()
    L = open(o_wpfile).readlines()
    ppW = mpl.backends.backend_pdf.PdfPages('CWP_W.pdf')
    sys.stdout.write('CWP_W.pdf\n')
    if isref:
        ppR = mpl.backends.backend_pdf.PdfPages('CWP_R.pdf')
        sys.stdout.write('CWP_R.pdf')
    for chan in CHAN:
        cb = 0.0
        stat_label = []
        stat_posit = []
        # Read o_wpinversion
        for l in L:
            items = l.strip().split()
            sac.read(items[0])
            if sac.kcmpnm[2] != chan[2]:
                continue
            stat_label.append(sac.kstnm)
            i1 = int(items[3])
            i2 = int(items[4])
            npts = float(i2-i1)    
            stat_posit.append(cb+npts/2.0)
            cb += npts
        if not len(stat_label):
            sys.stderr.write('WARNING: No channel %s in %s\n'%(chan,o_wpfile))
            continue
        # Read predfile
        ifile = predfile+'_'+chan
        L2 = open(ifile).readlines()
        ncol = len(L2[0].strip().split())
        if ncol < 3 and isref:
            print('Warning No ref solution in %s'%ifile)
            isref = 0
        else:
            Wref = []
        Wdat = []
        Wsyn = []
        for l in L2:
            items = l.strip().split()
            Wdat.append(float(items[0])*1000.0)
            Wsyn.append(float(items[1])*1000.0)
            if isref:
                if len(items)<3:
                    raise IOError('ERROR: error reading %s\n'%(ifile))
                Wref.append(float(items[2])*1000.0)
        t = np.arange(0,len(Wdat),dtype='float')
        # Display
        fig=plt.figure(figsize=CWP_FIGSIZE) 
        fig.subplots_adjust(left=0.08,bottom=0.12,right=0.96,top=0.88,wspace=0.2,hspace=0.2)
        plt.plot(t,Wdat,'k')
        plt.plot(t,Wsyn,'r')
        ymin = 1.1*min(Wdat)
        ymax = 1.1*max(Wdat)
        for stnm,x in zip(stat_label,stat_posit):
            plt.text(x,ymax*0.6,stnm,rotation=90,fontsize=16,fontstyle='italic')
        plt.ylim([ymin,ymax])
        plt.xlim([0,t[-1]])
        plt.xlabel('time, sec')
        plt.ylabel('displacement, mm')
        plt.title('Data fit, W Phase solution, %s'%chan[2])
        ppW.savefig(papertype='a4',orientation='landscape')
        plt.close()
        if isref:
            fig = plt.figure(figsize=CWP_FIGSIZE) 
            fig.subplots_adjust(left=0.08,bottom=0.12,right=0.96,top=0.88,wspace=0.2,hspace=0.2)
            plt.plot(t,Wdat,'k')
            plt.plot(t,Wref,'r')
            ymin = 1.1*min(Wdat)
            ymax = 1.1*max(Wdat)
            for stnm,x in zip(stat_label,stat_posit):
                plt.text(x,ymax*0.6,stnm,rotation=90,fontsize=16,fontstyle='italic')
            plt.ylim([ymin,ymax])
            plt.xlim([0,t[-1]])
            plt.xlabel('time, sec')
            plt.ylabel('displacement, mm')
            plt.title('Data fit, Reference solution, %s'%chan[2])
            ppR.savefig(papertype='a4',orientation='landscape')
            plt.close()
        sys.stdout.write('\n')
    ppW.close()
    if isref:
        ppR.close()
        

if __name__=='__main__':
    main(sys.argv)
            
