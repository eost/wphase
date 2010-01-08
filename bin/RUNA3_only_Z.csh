#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH
##################################

set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract_only_Z.csh
set CALC    = ${BIN}/calc_fast_synths_only_Z.csh
set PREPARE = ${BIN}/prepare_wp_ZNE.csh
set WPINVER = ${BIN}/wpinversion_ZNE


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
        set gf_dir   = `echo $tmp | ${HEAD} -1 | ${CUT} -d: -f2`
endif
${RM} -f i_tmp 

$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -gfdir ${gf_dir} \
  	 -pdata fort.15.noth 

${CP} p_wpinversion p_wpinversion.noth
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set ths =  "5.0 3.0 0.9"

foreach th ($ths)

    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th}\
    -gfdir ${gf_dir} -old
    if ($status == 1) exit(1)
    ${CP} -f o_wpinv.th_$th o_wpinversion        
end

${CP} WCMTSOLUTION.th_${th} WCMTSOLUTION
${CP} p_wpinversion.th_${th} p_wpinversion
${CP} LOG/wpinversion.th_${th}.log LOG/wpinversion.log
