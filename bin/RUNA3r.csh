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

${MKDIR} -p SYNTH; ${RM} -rf SYNTH/*
$EXTRACT
$PREPARE
if $status exit(1)

$ECHO "COMMAND LINE: $WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth -med ${my_argv}"
$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -pdata fort.15.noth -med ${my_argv}

${CP} p_wpinversion p_wpinversion.noth
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set ths = "5.0 3.0 0.9"
set nr  = 2.0

foreach th ($ths)
    $ECHO -e "\nCOMMAND LINE: $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -old ${my_argv}"
    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} \
    -old ${my_argv}
    if ($status == 1) exit(1)
    ${CP} -f o_wpinv.th_$th o_wpinversion
end

$ECHO -e "\nCOMMAND LINE: $WPINVER -nr $nr -ifil o_wpinversion -ofil o_wpinv.r_${nr} -log LOG/wpinversion.r_${nr}.log -ps p_wpinversion.r_${nr} -osyndir SYNTH -ocmtf  WCMTSOLUTION.r_${nr} -old ${my_argv}"
$WPINVER -nr $nr -ifil o_wpinversion -ofil o_wpinv.r_${nr} \
    -log LOG/wpinversion.r_${nr}.log -ps p_wpinversion.r_${nr} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.r_${nr} -old ${my_argv}

if ($status == 1) exit(1)
    ${CP} -f o_wpinv.r_${nr} o_wpinversion

${CP} WCMTSOLUTION.r_${nr} WCMTSOLUTION
${CP} p_wpinversion.r_${nr} p_wpinversion
${CP} LOG/wpinversion.r_${nr}.log LOG/wpinversion.log

${ECHO} -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
${ECHO} "              fort.15 fort.15_LHZ fort.15_LHN fort.15_LHE"
