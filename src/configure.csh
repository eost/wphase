#!/bin/csh -f

if ( ! $?WPHASE_HOME ) then
    echo "*** ERROR : The WPHASE_HOME environment variable is not defined"
    exit(1)
else
    echo "WPHASE_HOME: $WPHASE_HOME"
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

set FAST_SYNTH_Z         = \$BIN/fast_synth_only_Z
set FAST_SYNTH           = \$BIN/fast_synth
set FAST_SYNTH_ROT       = \$BIN/fast_synth_rot

set DELTA_T              = \$BIN/delta_t
set READCMT              = \$BIN/readcmt
set SYNTHS               = \$BIN/synth_v5
set CONV_VBB             = \$BIN/conv_vbb

set MSEEDLST             = \$BIN/mseedlst
set SACLST               = \$BIN/saclst

#GMT
set PSTEXT               = \$GMT_BIN/pstext
set PSCOAST              = \$GMT_BIN/pscoast
set PSXY                 = \$GMT_BIN/psxy
set PSBASEMAP            = \$GMT_BIN/psbasemap
set GMTSET               = \$GMT_BIN/gmtset
set PSHISTOGRAM          = \$GMT_BIN/pshistogram

# commands
EOF

set coms = "echo ls ln cp mv rm sed cat pwd grep date mkdir cut awk find head tail sort touch od tr wc seq xargs paste dirname basename ps2pdf"

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
    if ! -e $compath echo "WARNING: Command $com not found in $PATH"
end
