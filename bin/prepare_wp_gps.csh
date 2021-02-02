#!/bin/csh -f
# Command lines examples:
# prepare_wp.csh 
# prepare_wp.csh Z
# prepare_wp.csh -a
#

source $WPHASE_HOME/bin/WP_HEADER_GPS.CSH
set LOG        = LOG
set DATA       = DATA
set DATA_org   = DATA_org
set BEFORE     = 300          # number of points before origin
set AFTER      = 1200         # number of points after origin
set FILT_ORDER = `$GREP filt_order i_master | $CUT -d ':' -f2 | $AWK '{print $1}'`
set FILT_CF1   = `$GREP filt_cf1 i_master | $CUT -d ':' -f2 | $AWK '{print $1}'`
set FILT_CF2   = `$GREP filt_cf2 i_master | $CUT -d ':' -f2 | $AWK '{print $1}'`
set FILT_PASS  = `$GREP filt_pass i_master | $CUT -d ':' -f2 | $AWK '{print $1}'`
set CMPS       = "Z N E 1 2"

# command line arguments 
set trim_flag = '-u'
set my_argv = ($ARGV)
if ($#my_argv > 0) then
    set CMPS = $my_argv
    if ( $CMPS == '-a' ) then
	set CMPS = "Z N E 1 2"
	set trim_flag = "-a"
    endif
endif


if ( ! -e ${DATA_org} ) then
   ${ECHO} "Missing ${DATA_org}"
   exit 1
endif
 
if ( -e ${DATA} ) $RM -rf ${DATA}
${MKDIR} ${DATA}
${MKDIR} -p ${LOG}

###############################################
# Renaming files, cutting and filtering (sac) #
${ECHO} "Renaming files, cutting and filtering (sac)..."
cd ${DATA_org}
foreach file (*.sac)
    ${SAC} << EOF>/dev/null
    r $file
    SETBB TEMP0 $file
    SETBB TEMP1 &$file,NZYEAR
    SETBB TEMP2 %TEMP1%.
    SETBB TEMP3 %TEMP2%&$file,NZJDAY
    SETBB TEMP4 %TEMP3%.
    SETBB TEMP5 %TEMP4%&$file,NZHOUR
    SETBB TEMP6 %TEMP5%.
    SETBB TEMP7 %TEMP6%&$file,NZMIN
    SETBB TEMP8 %TEMP7%.
    SETBB TEMP9 %TEMP8%&$file,NZSEC
    SETBB TEMP10 %TEMP9%.0000.
    SETBB TEMP11 %TEMP10%&$file,KNETWK
    SETBB TEMP12 %TEMP11%.
    SETBB TEMP13 %TEMP12%&$file,KSTNM
    SETBB TEMP14 %TEMP13%..
    SETBB TEMP15 %TEMP14%&$file,KCMPNM
    SETBB TEMP16 %TEMP15%.
    SETBB TEMP17 %TEMP16%X.1sps.SAC
    cut o -${BEFORE} ${AFTER}
    r $file
    bp c ${FILT_CF1} ${FILT_CF2} n ${FILT_ORDER} p ${FILT_PASS}
    w %TEMP17
    q
EOF
end
cd ..
${MV} ${DATA_org}/*.SAC ${DATA}

###############################################
# Trimming files                              #
${LS} ${DATA}/*.SAC >! ${DATA}/_sac_files_list
$TRIM_SAC_FILES i_master ${DATA}/_sac_files_list scr_dat_fil_list $trim_flag # Use of "-u" will allow only one network per-channel

################################################
# Synthetics preparatioN:convolution and filter #
${ECHO} "Synthetics convolution and filter...  "
${PREP_KERNELS} scr_dat_fil_list l
if $status exit(1)

################################################
# Creating input file for inversion ...        #
${CUT} -d' ' -f1  scr_dat_fil_list >! i_wpinversion
