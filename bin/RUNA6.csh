#!/bin/csh -f

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
    echo "*** ERROR ($0) ***"
    echo "Syntax: $0 [wZ wN wE]"
    exit
endif

set BIN     = $WPHASE_HOME/bin
set EXTRACT = ${BIN}/extract.csh
set CALC    = ${BIN}/calc_fast_synths.csh
set PREPARE = ${BIN}/prepare_wp.csh
set WPINVER = ${BIN}/wpinversion

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

# ${RM} -rf GF
# ${CP} -f /home/zac/WP6/run_test_martinique07/*.dec.bp.int DATA/
# ${CP} -rf /home/zac/WP6/run_test_martinique07/GF ./


$WPINVER -log LOG/wpinversion.noth.log -osyndir SYNTH -gfdir ${gf_dir} \
	 -pdata fort.15.noth -wn ${wN} -we ${wE} -wz ${wZ} -nt -med

${CP} p_wpinversion p_wpinversion.noth
${CP} o_wpinversion o_wpinversion.noth
${CP} WCMTSOLUTION WCMTSOLUTION.noth

${CP} o_wpinversion o_wpinv
set ths =  "5.0 4.0 3.0 2.0 1.0 0.9"


foreach th ($ths)

    $WPINVER -th ${th} -ifil o_wpinversion -ofil o_wpinv.th_${th} \
    -log LOG/wpinversion.th_${th}.log -ps p_wpinversion.th_${th} \
    -osyndir SYNTH -ocmtf  WCMTSOLUTION.th_${th} -gfdir ${gf_dir} \
    -wn ${wN} -we ${wE} -wz ${wZ} -nt -old

#     set NBSTA = `${CAT} o_wpinversion | ${WC}  -l`
#     if ( $NBSTA < 20 ) then
#  	${ECHO} "${NBSTA} stations : stop for th=${th}"
#  	break
#     endif
    ${CP} -f o_wpinv.th_$th o_wpinversion        
end

${CP} WCMTSOLUTION.th_${th} WCMTSOLUTION
${CP} p_wpinversion.th_${th} p_wpinversion
${CP} LOG/wpinversion.th_${th}.log LOG/wpinversion.log


echo "\nOutput files: o_wpinversion WCMTSOLUTION p_wpinversion"
echo "              fort.15 fort.15_LHZ fort.15_LHL fort.15_LHT"
