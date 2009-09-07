set      gf_path = $GF_PATH
set      wp_home = $WPHASE_HOME
set      gmt_bin = $GMT_BIN
set      rdseed  = $RDSEED
set      home    = $HOME

unsetenv *

setenv   HOME      $home
setenv   ARGV      "$argv"
setenv   DISPLAY   :0.0
setenv   GF_PATH   $gf_path
setenv   WPHASE_HOME   $wp_home
setenv   GMT_BIN   $gmt_bin
setenv   RDSEED    $rdseed

unset    *

setenv   LC_ALL en_US

set BIN                  = $WPHASE_HOME/bin
set TRIM_SAC_FILES	 = $BIN/trim_sac_files
set TRIM_SAC_FILES_QRT	 = $BIN/trim_sac_files_qrt
set MAKE_RESP_TABLE      = $BIN/make_resp_lookup_table

set REC_DEC_FILT         = $BIN/rec_dec_filt_v1
set SYN_CONV_FILT        = $BIN/syn_conv_filt_v1
set ROT_HORIZ_CMP        = $BIN/rot_horiz_cmp

set FAST_SYNTH_Z         = $BIN/fast_synth_only_Z
set FAST_SYNTH           = $BIN/fast_synth

set DELTA_T              = $BIN/delta_t
set READCMT              = $BIN/readcmt
set SYNTHS               = $BIN/synth_v5
set CONV_VBB             = $BIN/conv_vbb

set SACLST               = $BIN/saclst

#GMT
set PSTEXT               = $GMT_BIN/pstext
set PSCOAST              = $GMT_BIN/pscoast
set PSXY                 = $GMT_BIN/psxy
set PSBASEMAP            = $GMT_BIN/psbasemap
set GMTSET               = $GMT_BIN/gmtset
set PSHISTOGRAM          = $GMT_BIN/pshistogram

# commands
set ECHO = /bin/echo
set LS = /bin/ls
set LN = /bin/ln
set CP = /bin/cp
set MV = /bin/mv
set RM = /bin/rm
set SED = /bin/sed
set CAT = /bin/cat
set PWD = /bin/pwd
set GREP = /bin/grep
set DATE = /bin/date
set MKDIR = /bin/mkdir
set CUT = /usr/bin/cut
set AWK = /usr/bin/awk
set FIND = /usr/bin/find
set HEAD = /usr/bin/head
set TAIL = /usr/bin/tail
set SORT = /usr/bin/sort
set TOUCH = /usr/bin/touch
set OD = /usr/bin/od
set WC = /usr/bin/wc
set SEQ = /usr/bin/seq
set XARGS = /usr/bin/xargs
set PASTE = /usr/bin/paste
set DIRNAME = /usr/bin/dirname
set BASENAME = /usr/bin/basename
set PS2PDF = /usr/bin/ps2pdf
