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
set CALC    = ${BIN}/calc_fast_synths_ZNE.csh
set PREPARE = ${BIN}/prepare_wp_ZNE.csh
set WPINVER = ${BIN}/wpinversion_ZNE
set TRACES  = ${BIN}/traces_6t_global.gmt 

${RM} -rf SYNTH
${MKDIR} SYNTH

if (-e i_master) then
        ${GREP} -v "^#" i_master >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif

set CMTFILE  = `${GREP} CMTFILE i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
${RM} -f i_tmp

$EXTRACT -a
$CALC
$PREPARE
$WPINVER -osyndir SYNTH 

${RM} -f page_6t*.ps page_6t*.pdf
${TRACES} ${CMTFILE} SYNTH
