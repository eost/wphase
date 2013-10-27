#!/bin/csh -f

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

if ( ! $?WPHASE_HOME ) then
    echo "*** ERROR : The WPHASE_HOME environment variable is not defined"
    exit(1)
endif
set SRC = `pwd`
set w_home1 = `dirname ${WPHASE_HOME}`/`basename $WPHASE_HOME`
if ( ! -e ${w_home1}/bin) then
    echo "*** ERROR : The WPHASE_HOME environment variable is incorrect" 
    echo "***         ($w_home1/bin does not exist)  "
    exit(1)
endif
cd $w_home1
set w_home2 = `pwd`
if ( $w_home1 != $w_home2 ) then
    echo "*** ERROR : The WPHASE_HOME environment variable must be an absolute path"
    exit(1)
endif
cd $SRC
cd ..
set w_home2 = `pwd`
if ( $w_home1 != $w_home2 ) then
    echo "*** WARNING ***"
    echo "The WPHASE_HOME environment variable is $WPHASE_HOME"
    echo "While the current wphase directory is $w_home2"
    echo "Continue will write ${w_home1}/bin/WP_HEADER.CSH"
    echo "Continue? [Y/(N)]"
    set rep = $<
    if ( $rep != "Y" && $rep != "y" ) then
	echo "Stop"
	exit(1)
    endif
endif

set BIN    = $WPHASE_HOME/bin
set o_file = $BIN/WP_HEADER.CSH

echo "Creating  $BIN/WP_HEADER.CSH"

rm -f ${o_file}

cat << EOF > ${o_file}
set      gf_path = \$GF_PATH
set      wp_home = \$WPHASE_HOME
set      gmt_bin = \$GMT_BIN
set      rdseed  = \$RDSEED
set      home    = \$HOME

unsetenv *

setenv   HOME      \$home
setenv   ARGV      "\$argv"
setenv   DISPLAY   :0.0
setenv   GF_PATH   \$gf_path
setenv   WPHASE_HOME   \$wp_home
setenv   GMT_BIN   \$gmt_bin
setenv   RDSEED    \$rdseed

unset    *

setenv   LC_ALL en_US

set BIN                  = \$WPHASE_HOME/bin
set TRIM_SAC_FILES	 = \$BIN/trim_sac_files
set TRIM_SAC_FILES_QRT	 = \$BIN/trim_sac_files_qrt
set MAKE_RESP_TABLE      = \$BIN/make_resp_lookup_table

set REC_DEC_FILT         = \$BIN/rec_dec_filt
set SYN_CONV_FILT        = \$BIN/syn_conv_filt
set ROT_HORIZ_CMP        = \$BIN/rot_horiz_cmp

set PREP_KERNELS         = \${BIN}/prep_kernels
set PREP_KERNELS_only_Z  = \${BIN}/prep_kernels_only_Z

set FAST_SYNTH_Z         = \$BIN/fast_synth_only_Z
set FAST_SYNTH           = \$BIN/fast_synth
set FAST_SYNTH_ROT       = \$BIN/fast_synth_rot

set DELTA_T              = \$BIN/delta_t
set READCMT              = \$BIN/readcmt
set SYNTHS               = \$BIN/synth_v6
set CONV_VBB             = \$BIN/conv_vbb

set SACLST               = \$BIN/saclst
set DECIMATE             = \$BIN/decim_one_sac_file_to_1sps
#GMT
set PSTEXT               = \$GMT_BIN/pstext
set PSCOAST              = \$GMT_BIN/pscoast
set PSXY                 = \$GMT_BIN/psxy
set PSBASEMAP            = \$GMT_BIN/psbasemap
set GMTSET               = \$GMT_BIN/gmtset
set PSHISTOGRAM          = \$GMT_BIN/pshistogram

# commands
EOF

set coms = "echo ls ln cp mv rm sed cat pwd grep expand date mkdir cut awk find head tail sort touch od tr wc seq xargs paste dirname basename ps2pdf gs"

set dirs = '/usr/bin /bin'
unset WHICH
foreach dir ($dirs)
	if -e $dir/which set WHICH = $dir/which
end
if ! $?WHICH then
	echo "ERROR: Command 'which' not found in $dirs"
        exit 1
endif

foreach com ($coms)
    unset compath
    set COM        = `echo $com | tr "[:lower:]" "[:upper:]"` 
    set compath    = `$WHICH $com`
    echo "set $COM = $compath" >> ${o_file}
    if ( ! -e $compath ) echo "WARNING: Command $com not found in $PATH"
end
