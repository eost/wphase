#!/bin/csh -f
#
# W phase package - RUNA3r for L, T, Z components
#
# Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#

source $WPHASE_HOME/bin/WP_HEADER.CSH
##################################

set my_argv = ($ARGV)
if ($#my_argv < 1) then
    set wZ = 1.
    set wL = 1.
    set wT = 1.
else if ($#my_argv == 3) then
    set wZ = $my_argv[1]
    set wL = $my_argv[2]
    set wT = $my_argv[3]
else
    $ECHO "*** ERROR ($0) ***"
    $ECHO "Syntax: $0 [wZ wL wT]"
    exit
endif

set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract.csh
set CALC    = ${BIN}/calc_fast_synths_LTZ.csh
set PREPARE = ${BIN}/prepare_wp_LTZ.csh
set WPINVER = ${BIN}/wpinversion_LTZ

${RM} -rf SYNTH
${MKDIR} SYNTH

$EXTRACT
$CALC
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

$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -gfdir ${gf_dir} \
	 -pdata fort.15.noth -wl ${wL} -wt ${wT} -wz ${wZ} -med

${CP} p_wpinversion p_wpinversion.noth
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set ths = "5.0 3.0 0.9"
set nr  = 2.0

foreach th ($ths)
    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -gfdir ${gf_dir} \
    -wl ${wL} -wt ${wT} -wz ${wZ} -old
    if ($status == 1) exit(1)
    ${CP} -f o_wpinv.th_$th o_wpinversion
end

$WPINVER -nr $nr -ifil o_wpinversion -ofil o_wpinv.r_${nr} \
    -log LOG/wpinversion.r_${nr}.log -ps p_wpinversion.r_${nr} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.r_${nr} -gfdir ${gf_dir} \
    -wl ${wL} -wt ${wT} -wz ${wZ}  -old

if ($status == 1) exit(1)
    ${CP} -f o_wpinv.r_${nr} o_wpinversion

${CP} WCMTSOLUTION.r_${nr} WCMTSOLUTION
${CP} p_wpinversion.r_${nr} p_wpinversion
${CP} LOG/wpinversion.r_${nr}.log LOG/wpinversion.log


${ECHO} -e "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
${ECHO} "              fort.15 fort.15_LHZ fort.15_LHL fort.15_LHT"
