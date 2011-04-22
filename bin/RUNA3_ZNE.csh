#!/bin/csh -f
#
# W phase package - RUNA3 for Z, N, E components
#
# Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#

source $WPHASE_HOME/bin/WP_HEADER.CSH
##################################

set my_argv = ($ARGV)
if ($#my_argv < 1) then
    set wZ = 1.
    set wN = 1.
    set wE = 1.
else if ($#my_argv == 3) then
    set wZ = $my_argv[1]
    set wN = $my_argv[2]
    set wE = $my_argv[3]
else
    $ECHO "*** ERROR ($0) ***"
    $ECHO "Syntax: $0 [wZ wN wE]"
    exit
endif

set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract.csh
set CALC    = ${BIN}/calc_fast_synths_ZNE.csh
set PREPARE = ${BIN}/prepare_wp_ZNE.csh
set WPINVER = ${BIN}/wpinversion_ZNE

${RM} -rf SYNTH
${MKDIR} SYNTH

$EXTRACT
$CALC
if $status exit(1)
$PREPARE

if (-e i_master) then
        ${GREP} -v "^#" i_master >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif
set gf_dir   = "./GF"
set tmp      = `${GREP} GFDIR   i_tmp`
if ! $status then
        set gf_dir   = `$ECHO $tmp | ${HEAD} -1 | ${CUT} -d: -f2`
endif
${RM} -f i_tmp

$ECHO $WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -gfdir ${gf_dir} \
	 -pdata fort.15.noth -wn ${wN} -we ${wE} -wz ${wZ} -med

$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -gfdir ${gf_dir} \
	 -pdata fort.15.noth -wn ${wN} -we ${wE} -wz ${wZ} -med

${CP} p_wpinversion p_wpinversion.noth
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set ths = "5.0 3.0 0.9"

foreach th ($ths)
    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -gfdir ${gf_dir} \
    -wn ${wN} -we ${wE} -wz ${wZ} -old
    if ($status == 1) exit(1)
    ${CP} -f o_wpinv.th_$th o_wpinversion
end

${CP} WCMTSOLUTION.th_${th} WCMTSOLUTION
${CP} p_wpinversion.th_${th} p_wpinversion
${CP} LOG/wpinversion.th_${th}.log LOG/wpinversion.log

${ECHO} -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
${ECHO} "              fort.15 fort.15_LHZ fort.15_LHL fort.15_LHT"
