#!/bin/csh -f
#
# W phase package - RUNA3 for L, T, Z components
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
set PREPARE = ${BIN}/prepare_wp_LTZ.csh
set WPINVER = ${BIN}/wpinversion_LTZ


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
set ths =  "5.0 3.0 0.9"


foreach th ($ths)
    $ECHO -e "\nCOMMAND LINE: $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -old ${my_argv}"
    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -old ${my_argv}
    if ($status == 1) exit(1)
    ${CP} -f o_wpinv.th_$th o_wpinversion
end

${CP} WCMTSOLUTION.th_${th} WCMTSOLUTION
${CP} p_wpinversion.th_${th} p_wpinversion
${CP} LOG/wpinversion.th_${th}.log LOG/wpinversion.log


$ECHO -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
$ECHO "                fort.15 fort.15_LHZ fort.15_LHL fort.15_LHT"
