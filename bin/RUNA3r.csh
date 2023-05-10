#!/bin/csh -f
#
# W phase package - RUNA3r (mediane, rms and amplitude ratio screening)
#
# Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#

source $WPHASE_HOME/bin/WP_HEADER.CSH
##################################
if ( ! -e i_master ) then
   ${ECHO} "Error: file  i_master not available"
   exit
endif

set my_argv = ($ARGV)
set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract.csh
set PREPARE = ${BIN}/prepare_wp.csh
set WPINVER = ${BIN}/wpinversion

set median  = "-med"
set p2p_flag = `$GREP ^P2P_SCREENING i_master | $CUT -d':' -f2`
set p2p_yes = `$ECHO $p2p_flag | $CUT -d':' -f2 | $GREP YES`
set p2p_scr = `$ECHO $p2p_flag | $CUT -d':' -f2 | $SED 's/YES//'`
if ( $#p2p_flag != 0 && $#p2p_yes == 0 ) then
    set median = ""
else if ( $#p2p_scr != 0 ) then
    set median = "${median} $p2p_scr"
endif

set ths     = "5.0 3.0 0.9"
set rms_scr = `$GREP ^RMS_SCREENING i_master| $CUT -d':' -f2`
if ( $#rms_scr != 0 ) set ths = "$rms_scr"
 

if ( -e SYNTH ) ${RM} -rf SYNTH 
${MKDIR} SYNTH; 

$EXTRACT
$PREPARE
if $status exit 1

set version = "Version: r253"
set screening = "Screening: $median"
set comgfdir  = "GF_PATH: $GF_PATH"

$ECHO "COMMAND LINE: $WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth $median \
	${my_argv} -comment '$version' -comment '$screening' -comment '$comgfdir'"
$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth $median \
	${my_argv} -comment "$version" -comment "$screening" -comment "$comgfdir"
if $status exit 1

${CP} p_wpinversion.ps p_wpinversion.noth.ps
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set screening = "$screening -th"
set nr  = 2.0

foreach th ($ths)
	set screening = "$screening  $th"
    $ECHO -e "\nCOMMAND LINE: $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
	-log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
	-osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -old ${my_argv} \
	-comment '$version' -comment '$screening' -comment '$comgfdir'" 
    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th}.ps \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -old ${my_argv} \
	-comment "$version" -comment "$screening" -comment "$comgfdir" 
    if $status exit 1
    ${CP} -f o_wpinv.th_$th o_wpinversion
end


$ECHO -e "\nCOMMAND LINE: $WPINVER -nr $nr -ifil o_wpinversion -ofil o_wpinv.r_${nr} -log LOG/wpinversion.r_${nr}.log -ps p_wpinversion.r_${nr} -osyndir SYNTH -ocmtf  WCMTSOLUTION.r_${nr} -old ${my_argv}"
$WPINVER -nr $nr -ifil o_wpinversion -ofil o_wpinv.r_${nr} \
    -log LOG/wpinversion.r_${nr}.log -ps p_wpinversion.r_${nr} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.r_${nr} -old ${my_argv}
if $status exit 1

${CP} -f o_wpinv.r_${nr} o_wpinversion

${CP} WCMTSOLUTION.r_${nr} WCMTSOLUTION
${CP} p_wpinversion.r_${nr} p_wpinversion
${CP} LOG/wpinversion.r_${nr}.log LOG/wpinversion.log

${ECHO} -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
${ECHO} "              fort.15 fort.15_LHZ fort.15_LHN fort.15_LHE"
