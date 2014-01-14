#!/usr/bin/env python
#*-* coding: iso-8859-1 *-*

############################################################################
#
#	              W phase source inversion package 	            
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

# W PHASE TRACES

IMASTER       = 'i_master'      # IMASTER FILENAME
O_WPINVERSION = 'o_wpinversion' # o_wpinversion filename
LOGDIR        = 'LOG'
LENGTH        = 1500;           # Traces lenght
DLAT, DLON    = 20., 20.        # Half-size of the map region
OPDFFILE      = 'wp_pages.pdf'

FIGSIZE   = [11.69,8.270]
#FIGSIZE   = [5.84,4.135]
YLIM_AUTO = True
YLIMFIXED = [-9,12] # Y lim if YLIM_AUTO = False
NC = 3 # Number of columns
NL = 5 # Number of lines

plotparams = {'backend': 'pdf',
			  'axes.labelsize': 10,
			  'text.fontsize': 10,
			  'maskedarray': 'obsolete',
			  'numerix': 'obsolete',
			  'xtick.labelsize': 10,
			  'ytick.labelsize': 10,
			  'legend.pad': 0.1,     # empty space around the legend box
			  'legend.fontsize': 10,
			  'lines.markersize': 6,
			  'font.size': 10,
			  'savefig.dpi': 200,
			  'keymap.all_axes': 'a',
			  'keymap.back': ['left', 'c', 'backspace'],
			  'keymap.forward': ['right', 'v'],
			  'keymap.fullscreen': 'f',
			  'keymap.grid': 'g',
			  'keymap.home': ['h', 'r', 'home'],
			  'keymap.pan': 'p',
			  'keymap.save': 's',
			  'keymap.xscale': ['k', 'L'],
			  'keymap.yscale': 'l',
			  'keymap.zoom': 'o',                  
			  'path.snap': True,
			  'savefig.extension': 'pdf',
			  'pdf.compression': 9,
			  'figure.figsize': FIGSIZE}


# Import modules
import matplotlib
matplotlib.use('PDF')
import os,sys,re
import getopt as go
import shutil as sh
import pylab as pyl
from subprocess import call

pyl.rcParams.update(plotparams)

# Environment variables
WPHOME = os.path.expandvars('$WPHASE_HOME')
print 'WPHASE_HOME is %s'%(WPHOME)
if WPHOME[-1] != '/':
	WPHOME += '/'
if LOGDIR[-1] != '/':
	LOGDIR += '/'
SYNTHS = WPHOME+'bin/synth_v6'

# Function - Class definitions

def parse_config(cfg_file):
	config = {}
	try:
		config_lines = open(cfg_file, 'r').readlines()
		for line in config_lines:
			if line.find('#')==0:
				continue
			if line.rstrip():
				key,value = line.strip().split(':')
				config[key.strip()]=value.strip()
	except:
		sys.stderr.write('Error: format  %s\n'%cfg_file)
		sys.exit(1)
	return config

def unpack_c(chararray):
	S = ''
	for c in chararray:
		if c == ' ':
			break
		S+=c
	return S

class Sac:
	def __init__(self):
		self.delta  =  -12345.
		self.b      =  -12345.
		self.o      =  -12345.
		self.stla   =  -12345.
		self.stlo   =  -12345.
		self.evla   =  -12345.
		self.evlo   =  -12345.
		self.cmpaz  =  -12345.
		self.cmpinc =  -12345.
		self.az     =  -12345.
		self.baz    =  -12345.
		self.gcarc  =  -12345.
		self.dist   =  -12345.
		self.nzyear =  -12345
		self.nzjday =  -12345
		self.nzhour =  -12345
		self.nzmin  =  -12345
		self.nzsec  =  -12345
		self.nzmsec =  -12345
		self.npts   =  -12345
		self.kstnm  = '-12345'
		self.kcmpnm = '-12345'
		self.knetwk = '-12345'
		self.khole  = '-12345'
		self.id     = self.knetwk+'_'+self.kstnm+'_'+self.khole+'_'+self.kcmpnm
		self.depvar =  []
	def rsac(self,FILE,np=-1,datflag=1):
		try:
			fid     = open(FILE,'rb')
			self.delta   = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(20,0)
			self.b       = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(28,0)
			self.o       = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(124,0)
			self.stla    = pyl.fromfile(fid,'float32',   1)[0]
			self.stlo    = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(140,0)
			self.evla    = pyl.fromfile(fid,'float32',   1)[0]
			self.evlo    = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(200,0)
			self.dist    = pyl.fromfile(fid,'float32',   1)[0]
			self.az      = pyl.fromfile(fid,'float32',   1)[0]
			self.baz     = pyl.fromfile(fid,'float32',   1)[0]
			self.gcarc   = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(228,0)
			self.cmpaz   = pyl.fromfile(fid,'float32',   1)[0]
			self.cmpinc  = pyl.fromfile(fid,'float32',   1)[0]
			fid.seek(280,0);
			self.nzyear  = pyl.fromfile(fid,  'int32',   1)[0]
			self.nzjday  = pyl.fromfile(fid,  'int32',   1)[0]
			self.nzhour  = pyl.fromfile(fid,  'int32',   1)[0]
			self.nzmin   = pyl.fromfile(fid,  'int32',   1)[0]
			self.nzsec   = pyl.fromfile(fid,  'int32',   1)[0]
			self.nzmsec  = pyl.fromfile(fid,  'int32',   1)[0]
			fid.seek(316,0)
			self.npts    = pyl.fromfile(fid,  'int32',   1)[0]
			fid.seek(440,0);
			self.kstnm   = unpack_c(pyl.fromfile(fid,'c',   8))
			fid.seek(464,0);
			self.khole   = unpack_c(pyl.fromfile(fid,'c',   8))
			fid.seek(600,0);
			self.kcmpnm  = unpack_c(pyl.fromfile(fid,'c',   8))
			self.knetwk  = unpack_c(pyl.fromfile(fid,'c',   8))
			fid.seek(632,0);
			if self.khole == '':
				self.khole = '--'
			self.id    = self.knetwk+'_'+self.kstnm+'_'+self.khole+'_'+self.kcmpnm
			if not datflag:
				fid.close()
				return
			np = int(np)
			if np < 0 or np > self.npts:
				np = self.npts
			if np > 0:
				self.depvar = pyl.array(pyl.fromfile(fid,'float32',np),dtype='d')
			fid.close()
			
		except IOError:
			sys.stderr.write('error reading file '+FILE+'!!!\n')
			sys.exit(1)

def rm(fd):
        if os.path.islink(fd) or os.path.isfile(fd):
                os.remove(fd)
                return 0
        elif os.path.isdir(fd):
                sh.rmtree(fd)
                return 0
        return 1

def change_label_size(ax,size=10.0):
	for l in ax.get_xticklabels() + ax.get_yticklabels():
		l.set_fontsize(size)

def show_basemap(ax,evla,evlo,stla,stlo,coords,m=None):	
	if not m:
		from mpl_toolkits.basemap import Basemap
		m = Basemap(llcrnrlon=evlo-DLON, llcrnrlat=evla-DLAT,
				    urcrnrlon=evlo+DLON, urcrnrlat=evla+DLAT,
				   projection='lcc',lat_1=evla-DLAT/2.,lat_2=evla+DLAT/2.,
				   lon_0=evlo, resolution ='h',area_thresh=50. )

		#m = Basemap(projection='ortho',lat_0=evla,lon_0=evlo,resolution='c')
	pos  = ax.get_position().get_points()
	W  = pos[1][0]-pos[0][0] ; H  = pos[1][1]-pos[0][1] ;		
	ax2 = pyl.axes([pos[1][0]-W*0.38,pos[0][1]+H*0.01,H*1.08,H*1.00])
	m.drawcoastlines(linewidth=0.5,zorder=900)
 	m.fillcontinents(color='0.75',lake_color=None)
 	m.drawparallels(pyl.arange(evla-DLAT,evla+DLAT,5.0),linewidth=0.2)
 	m.drawmeridians(pyl.arange(evlo-DLON,evlo+DLON,5.0),linewidth=0.2)
	m.drawmapboundary(fill_color='w')
 	xs,ys = m(coords[:,1],coords[:,0])
 	xr,yr = m(stlo,stla)
 	xc,yc = m(evlo,evla)
	m.plot(xs,ys,'o',color=(1.00000,  0.74706,  0.00000),ms=4.0,alpha=1.0,zorder=1000)
 	m.plot([xr],[yr],'o',color=(1,.27,0),ms=8,alpha=1.0,zorder=1001)
 	m.scatter([xc],[yc],c='b',marker=(5,1,0),s=120,zorder=1002)	
	return m

def show_polarmap(ax,az,dist,coords):
	distlabel  = 6371.0*pyl.arange(30.0,120.0,30.0)*pyl.pi/180.0
	pos  = ax.get_position().get_points()
	W  = pos[1][0]-pos[0][0] ; H  = pos[1][1]-pos[0][1] ;		
	ax2 = pyl.axes([pos[1][0]-W*0.3,pos[0][1]+H*0.38,H*0.73,H*0.73],polar=True) 
	ax2.plot(coords[:,2],coords[:,3],'yv',ms=4)
	ax2.plot(sacdata.az,sacdata.dist,'rv',ms=6)
	ax2.scatter([0.0],[0.0],c='r',marker=(5,1,0),s=70,zorder=1002)
	ax2.set_thetagrids([])
	ax2.set_rgrids(distlabel,labels='')

def usage(cmd):
	print 'usage: %s [option] (for help see %s -h)'%(cmd,cmd)

def disphelp(cmd,solfile,syndir):
	print 'Display W phase traces\n'
	usage(cmd)
	print '\nAll parameters are optional:'
	print '   -i, --icmtf          (w)cmt file name (default: %s)'%(solfile)
	print '   -d, --osyndir        output synthetic directory (default: %s)'%(syndir)
	print '\n   -h, --help         display this help and exit'
	print '\nReport bugs to: <zacharie.duputel@unistra.fr>'

if __name__ == '__main__':
	# Input parameters
	imaster = IMASTER
	length  = LENGTH;
	syndir  = 'SYNTH_traces'
	o_wpinversion = O_WPINVERSION
	nc = NC
	nl = NL
	solfile = None

	# Title
	conf  = parse_config(imaster)
	title = '_'.join(conf['EVNAME'].split())
	title += ',  filter = (%s, %s, %s, %s)'%(conf['filt_cf1'],conf['filt_cf2'],conf['filt_order'],conf['filt_pass']) 
	try:
		opts, args = go.gnu_getopt(sys.argv[1:],'i:d:h',["icmtf=","osydir=","help"])
	except go.GetoptError, err:
		sys.stderr.write('*** ERROR ***\n')
		usage(sys.argv[0])
		sys.exit(1)	
	for o, a in opts:
		if o == '-h' or o == '--help':
			disphelp(sys.argv[0],solfile,syndir)
			sys.exit(0)
		if o == '-i' or o=='--icmtf':
			solfile = a
			if not os.path.exists(solfile):
				print 'ERROR: no wcmtfile named %s'%(solfile)
				usage(sys.argv[0])
				sys.exit(1)
		if o == '-d' or o == '--osyndir':
			syndir = a
	if not solfile:
		for f in ['xy_WCMTSOLUTION','ts_WCMTSOLUTION','WCMTSOLUTION']:
			if os.path.exists(f):
				solfile = f
				break
		if not solfile:
			print 'ERROR: no available wcmtfile'
			usage(sys.argv[0])
			sys.exit(1)
	if os.path.exists(syndir) and syndir != '.' and syndir != './':
		rm(syndir)
	if syndir != '.' and syndir != './':
		os.mkdir(syndir)
	if not os.path.exists(LOGDIR):
		os.mkdir(LOGDIR)
	for l in os.listdir('.'):
		if l[:4]=='page' and l[-4:]=='.pdf':
			rm(l)
	# Compute synthetics
	cmd    = SYNTHS+' '+imaster+' '+solfile+' '+o_wpinversion+' '+syndir
	print cmd
	#status = call(cmd, shell=True, stdin=sys.stdin, stdout=sys.stdout);
	status = os.system(SYNTHS+' '+imaster+' '+solfile+' '+o_wpinversion+' '+syndir+' > %s_tmp_synths'%(LOGDIR))
	if status:
		print 'Error while running '+SYNTHS
		sys.exit(1)
	# Sac Objects
	sacdata = Sac()
	sacsynt = Sac()
	coords = []
	L = open(o_wpinversion).readlines()
	for l in L:
		sacf = l.strip().split()[0]
		sacdata.rsac(sacf,datflag=0)
		coords.append([sacdata.stla,sacdata.stlo,sacdata.az,sacdata.dist])
	coords = pyl.array(coords)
	# Display	
	print 'Input (W)CMTSOLUTION file is: %s'%(solfile)
	print 'Output synthetic directory is: %s'%(syndir)
	perpage = nl*nc
	ntot   = len(L)
	npages = pyl.ceil(float(ntot)/float(perpage))
	nchan = 1
	count = 1
	pages = 1
	fig = pyl.figure()
	fig.subplots_adjust(bottom=0.06,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.35)
	print '%d pages:'%(npages)
	pp = matplotlib.backends.backend_pdf.PdfPages(OPDFFILE)
	m = None
	for l in L:
		# Parse line
		items = l.strip().split()
		fic1 = items[0]
		sacdata.rsac(fic1)
		chan = sacdata.kcmpnm[0:3]
		loc  = sacdata.khole
		fic2 = syndir+'/%s.%s.%s.%s.complete_synth.bp.sac'\
		       %(sacdata.kstnm,sacdata.knetwk,chan,loc)
		sacsynt.rsac(fic2)		
		# pages
		if count > perpage:
			pyl.suptitle(title+ ',   p %d/%d'%(pages,npages), fontsize=16, y=0.95)
			ofic = 'page_W_%02d.pdf'%(pages)
			print ofic
			#fig.set_rasterized(True)
			pp.savefig(orientation='landscape',format='pdf')
			pyl.close()
			pages += 1
			count = 1
			fig = pyl.figure()
			fig.subplots_adjust(bottom=0.06,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.35)
		# Time - W phase window
		t1 = pyl.arange(sacdata.npts,dtype='double')*sacdata.delta + sacdata.b - sacdata.o
		t2 = pyl.arange(sacsynt.npts,dtype='double')*sacsynt.delta + sacsynt.b - sacsynt.o		
		wnb = float(items[5])
		wne = float(items[6])
		wtb = sacdata.b - sacdata.o + wnb * sacdata.delta
		wte = sacdata.b - sacdata.o + wne * sacdata.delta		
		# Plot trace
		ax = pyl.subplot(nl,nc,count)
		pyl.plot(t1,sacdata.depvar*1000.,'k')
		pyl.plot(t2,sacsynt.depvar*1000.,'r-')
		pyl.plot([wtb,wte],[0.,0.],'ro')
		# Axes limits
		B=wtb-150.0
		if B<0:
			B = 0.0
		pyl.xlim([B,B+length*sacsynt.delta])		
		if YLIM_AUTO:
			a    = pyl.absolute(sacsynt.depvar[:length]).max()*1000.
			ymin = -1.1*a
			ymax =  1.1*a
			ylims = [ymin,ymax]
		else:
			ylims = YLIMFIXED
		pyl.ylim(ylims)		
		# Annotations
		if sacdata.kcmpnm[2] == 'Z':
			label = r'%s %s %s %s $(\phi,\Delta) = %6.1f^{\circ}, %6.1f^{\circ}$'%(
					sacdata.knetwk,sacdata.kstnm, sacdata.kcmpnm, sacdata.khole,
			sacdata.az, sacdata.gcarc)
		else:
			label  = r'%s %s %s %s $(\phi,\Delta,\alpha) = %6.1f^{\circ},'
			label += '%6.1f^{\circ}, %6.1f^{\circ}$'
			label  = label%(sacdata.knetwk,sacdata.kstnm, sacdata.kcmpnm, sacdata.khole,
			sacdata.az, sacdata.gcarc, sacdata.cmpaz)	
		pyl.title(label,fontsize=10.0,va='center',ha='center')
		if not (count-1)%nc:
			pyl.ylabel('mm',fontsize=10)
		if (count-1)/nc == nl-1 or nchan+nc > ntot:
			pyl.xlabel('time, sec',fontsize=10) 
		pyl.grid()
 		try:
 			m = show_basemap(ax,sacdata.evla,sacdata.evlo,sacdata.stla,sacdata.stlo,coords,m)
			pass
 		except:
 			show_polarmap(ax,sacdata.az,sacdata.dist,coords)
 			print 'No basemap module'
		count += 1
		nchan += 1
	ofic = 'page_W_%02d.pdf'%(pages)
	print ofic
	#fig.set_rasterized(True)
	pyl.suptitle(title + ',    p %d/%d'%(pages,npages), fontsize=16, y=0.95)
	pp.savefig(orientation='landscape',format='pdf')
	pyl.close()
	pp.close()
