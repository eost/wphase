#!/usr/bin/python
# *-* coding: iso-8859-1 *-*

FIGSIZE   = [11.69,8.27]

import matplotlib
matplotlib.use('PDF')
import os,sys,re
import getopt as go
import pylab as pyl

# Environment variables
WPHOME = os.path.expandvars('$WPHASE_HOME')
print 'WPHASE_HOME is %s'%(WPHOME)

def unpack_c(chararray):
	S = ''
	for c in chararray:
		if c == ' ':
			break
		S+=c
	return S

class Sac:
	def __init__(self):
                self.kstnm  = '-12345'
		self.kcmpnm = '-12345'
	def rsac(self,FILE,np=-1):
		try:
			fid     = open(FILE,'rb')
			fid.seek(440,0);
			self.kstnm   = unpack_c(pyl.fromfile(fid,'c',   8))
			fid.seek(600,0);
			self.kcmpnm  = unpack_c(pyl.fromfile(fid,'c',   8))
			fid.close()
		except IOError:
			sys.stderr.write('ERROR: Reading file '+FILE+'\n')
			sys.exit(1)
        
def usage(cmd):
	print 'usage: %s [chan1 chan2 (default: LHZ LHN LHE)] [option] (for help see %s -h)'%(cmd,cmd)

def disphelp(cmd):
	print 'Make CWP traces\n'
	usage(cmd)
	print '\nAll parameters are optional:'
	print '\n   -i, --ifort15             input fort.15 file (e.g. fort.15, ts_fort.15, xy_fort.15)'
        print '\n   -n, --noref               no reference solution'
	print '\n   -h, --help                display this help and exit'
	print '\nReport bugs to: <zacharie.duputel@eost.u-strasbg.fr>'

if __name__=='__main__':
    # Input parameters
    try:
        opts, args = go.gnu_getopt(sys.argv[1:],'i:nh',["ifort15=","noref","help"])
    except:
        sys.stderr.write('*** ERROR ***\n')
        usage(sys.argv[0])
        sys.exit(1)
    o_wpfile = 'o_wpinversion'
    predfile = ''
    isref    = 1
    CHAN     = ['LHZ','LHN','LHE']
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
	    sys.stderr.write('Error: No fort.15 file named %s\n'%predfile)
	    sys.exit(1)
    if not len(predfile):
	    predfile = 'xy_fort.15'
	    if not os.path.exists(predfile):
		    predfile = 'ts_fort.15'
		    if not os.path.exists(predfile):
			    predfile = 'fort.15'
			    if not os.path.exists(predfile):
				    sys.stderr.write('Error: No fort.15 file found\n')
				    sys.exit(1)	    
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
	    sys.stdout.write('\n')
	    sys.stderr.write('\nError: No fort.15_ file for')
	    for chan in CHAN:
		    sys.stderr.write('%5s'%chan)
	    sys.stderr.write('\n')
	    sys.exit(1)
    sys.stdout.write('\nRead %s ...\n%s pages:\n'%(o_wpfile,count))

    # Main loop
    sac = Sac()
    L = open(o_wpfile).readlines()
    ppW = matplotlib.backends.backend_pdf.PdfPages('CWP_W.pdf')
    sys.stdout.write('CWP_W.pdf\n')
    if isref:
	    ppR = matplotlib.backends.backend_pdf.PdfPages('CWP_R.pdf')
	    sys.stdout.write('CWP_R.pdf')
    for chan in CHAN:
	    cb = 0.0
	    stat_label = []
	    stat_posit = []
	    # Read o_wpinversion
	    for l in L:
		    items = l.strip().split()
		    sac.rsac(items[0])
		    if sac.kcmpnm != chan:
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
		    print 'Warning No ref solution in %s'%ifile
		    isref = 0
	    else:
		    Wref = []
	    Wdat = []
	    Wsyn = []
	    for l in L2:
		    items = map(float,l.strip().split())
		    Wdat.append(items[0]*1000.0)
		    Wsyn.append(items[1]*1000.0)
		    if isref:
			    if len(items)<3:
				    sys.stderr.write('ERROR: error reading %s\n'%(ifile))
				    sys.exit(1)
			    Wref.append(items[2]*1000.0)
	    t = pyl.arange(0,len(Wdat),dtype='float')
	    # Display
	    fig=pyl.figure(figsize=FIGSIZE) 
	    fig.subplots_adjust(left=0.08,bottom=0.12,right=0.96,top=0.88,wspace=0.2,hspace=0.2)
	    pyl.plot(t,Wdat,'k')
	    pyl.plot(t,Wsyn,'r')
	    ymin = 1.1*min(Wdat)
	    ymax = 1.1*max(Wdat)
	    for stnm,x in zip(stat_label,stat_posit):
		    pyl.text(x,ymax*0.6,stnm,rotation=90,fontsize=16,fontstyle='italic')
	    pyl.ylim([ymin,ymax])
	    pyl.xlim([0,t[-1]])
	    pyl.xlabel('time, sec')
	    pyl.ylabel('displacement, mm')
	    pyl.title('Data fit, W Phase solution, %s'%chan)
	    ppW.savefig(papertype='a4',orientation='landscape')
	    pyl.close()
	    if isref:
		    fig = pyl.figure(figsize=FIGSIZE) 
		    fig.subplots_adjust(left=0.08,bottom=0.12,right=0.96,top=0.88,wspace=0.2,hspace=0.2)
		    pyl.plot(t,Wdat,'k')
		    pyl.plot(t,Wref,'r')
		    ymin = 1.1*min(Wdat)
		    ymax = 1.1*max(Wdat)
		    for stnm,x in zip(stat_label,stat_posit):
			    pyl.text(x,ymax*0.6,stnm,rotation=90,fontsize=16,fontstyle='italic')
		    pyl.ylim([ymin,ymax])
		    pyl.xlim([0,t[-1]])
		    pyl.xlabel('time, sec')
		    pyl.ylabel('displacement, mm')
		    pyl.title('Data fit, Reference solution, %s'%chan)
		    ppR.savefig(papertype='a4',orientation='landscape')
		    pyl.close()
	    sys.stdout.write('\n')
    ppW.close()
    if isref:
	    ppR.close()
	    
