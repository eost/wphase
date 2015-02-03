#!/bin/csh -f
#
# W phase package - RUNA3  (median and rms screening)
#
# Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#

source $WPHASE_HOME/bin/WP_HEADER.CSH
##################################

if ( ! -e i_master ) then
   ${ECHO} "Error: file  i_master not available"
   exit 1
endif

set my_argv = ($ARGV)
set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract.csh
set PREPARE = ${BIN}/prepare_wp.csh
set WPINVER = ${BIN}/wpinversion

set median  = "-med"
set p2p_scr = `$GREP ^P2P_SCREENING i_master| $CUT -d':' -f2`
if ( $#p2p_scr != 0 && $p2p_scr != "YES" ) set median = ""
set ths     = "5.0 3.0 0.9"
set rms_scr = `$GREP ^RMS_SCREENING i_master| $CUT -d':' -f2`
if ( $#rms_scr != 0 ) set ths = "$rms_scr"
 

if ( -e SYNTH ) ${RM} -rf SYNTH 
${MKDIR} SYNTH; 
$EXTRACT
if $status exit 1
$PREPARE
if $status exit 1

#if ( -e $WPHASE_HOME/.svn/entries ) set version = "${version}`$HEAD -4 $WPHASE_HOME/.svn/entries | $TAIL -1`"
set version = "Version: r251"
set screening = "Screening: $median"
set comgfdir  = "GF_PATH: $GF_PATH"

$ECHO "COMMAND LINE: $WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth $median \
	${my_argv} -comment '$version' -comment '$screening' -comment '$comgfdir'"
$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth $median \
	${my_argv} -comment "$version" -comment "$screening" -comment "$comgfdir"

${CP} p_wpinversion.ps p_wpinversion.noth.ps
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv


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

${CP} WCMTSOLUTION.th_${th} WCMTSOLUTION
${CP} p_wpinversion.th_${th}.ps p_wpinversion.ps
${CP} LOG/wpinversion.th_${th}.log LOG/wpinversion.log
#${RM} -rf GF

${ECHO} -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion.ps"
${ECHO} "              fort.15 fort.15_LHZ fort.15_LHN fort.15_LHE"
