'''
Input parameters for W-phase inversion
'''

import os 

# Time-shift grid-search parameters
TS_NIT   = 3  # Nb of iterations
TS_DT    = 4. # Initial time step
TSBOUNDS = [] # Bounds (empty=automatically determined from mb or Ms in the PDE line)
TS_OFILE = 'grid_search_ts_out'

# Centroid Lat/Lon grid-search parameters
XY_NIT   = 3   # Nb of iterations
XY_DX    = 0.4 # Intial samp. period
XY_NX    = 3   # Half_width = XY_NX*XY_DX
XY_NOPT  = 5   # Nb of optimal-points
XY_OFILE = 'grid_search_xy_out'

# Centroid Depth grid-search parameters
XYZ_NIT   = 1    # Nb of iterations
XYZ_DX    = 0.6  # Intial samp. period
XYZ_NX    = 1    # Half_width = XYZ_NX*XYZ_DX (if XYZ_NX=0: no Lat/Lon grid-seach is performed)
XYZ_NOPT  = 4    # Nb of optimal-points
DDEP      = 50.  # Delta depth ( Z_SEARCH within Z_INITIAL +/- DDEP )
MINDEP    = 11.5 
XYZ_OFILE = 'grid_search_xyz_out'

# Filenames
IMASTER         = 'i_master'      # IMASTER FILENAME
O_WPINVERSION   = 'o_wpinversion' # o_wpinversion filename
LOGDIR          = 'LOG'

# traces plot parameters
LENGTH_GLOBAL   = 3000 ;        # Traces lenght (teleseismic data)
LENGTH_REGIONAL = 1500 ;        # Traces lenght (regional data)
DLAT,DLON       = 20.,20.
OPDFFILE        = 'wp_pages.pdf'
TRACES_FIGSIZE  = [11.69,8.270]
YLIM_AUTO = True
YLIMFIXED = [-9,12] # Y lim if YLIM_AUTO = False
NC = 3 # Number of columns
NL = 5 # Number of lines

# cwp plot parameters
CWP_FIGSIZE   = [11.69,8.27]

TRACES_PLOTPARAMS = {'backend': 'pdf', 'axes.labelsize': 10,
                     'font.size': 10,
                     'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
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
                     'savefig.format': 'pdf',
                     'pdf.compression': 9,
                     'figure.figsize': TRACES_FIGSIZE}


# W-phase home directory
WPHOME = os.path.expandvars('$WPHASE_HOME')
print('WPHASE_HOME is %s'%(WPHOME))
if WPHOME[-1] != '/':
    WPHOME += '/'

GF_PATH = os.path.expandvars('$GF_PATH')
print('GF_PATH is %s'%(GF_PATH))

# Code version
VERSION = 'Version: r250'

# Path to binaries
BIN = WPHOME+'bin/'

WPINV_XY = BIN+'wpinversion_gs -imas i_master -ifil '+O_WPINVERSION
SYNTHS   = WPHOME+'bin/synth_v6'

