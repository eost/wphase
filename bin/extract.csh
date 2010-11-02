#!/bin/csh -f

########################################

source $WPHASE_HOME/bin/WP_HEADER.CSH
set CHANS = "LHE LHN LHZ"

set my_argv = ($ARGV)
if ($#my_argv < 1) then
    set trim_flag = "-u"
else if ($#my_argv == 1) then

    set trim_flag = $my_argv[1]
endif
set DATA = DATA
set LOG  = LOG
########################################
if (-e i_master) then
        ${GREP} -v "^#" i_master | ${SED} -e 's/ *:/:/' >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif
set seeds  = `${GREP} SEED       i_tmp |            ${CUT} -d: -f2`
set rundir = `${PWD}`
${RM} -f i_tmp i_seedfiles
foreach seed ($seeds)
    set dname = `${DIRNAME}  $seed`
    set fname = `${BASENAME}  $seed`
    cd $dname
    set dname = `${PWD}`
    cd $rundir
    ${ECHO} "${dname}/${fname}" >> i_seedfiles
end
set seeds = `${CAT} i_seedfiles`
########################################

${RM} -rf ${DATA} ${LOG}

${MKDIR} ${DATA} ${LOG}

########################################
# Extract, select and ${SORT} data
${ECHO} "Extracting data from seed           ( >! ${LOG}/_log_rdseed )"

${ECHO} '' >! ${LOG}/_log_rdseed

cd ${DATA}
foreach seed ($seeds)
    foreach chan ($CHANS)
	${ECHO} -e "${seed}\n\n\nd\n\n\n${chan}\n\n\n\nn\ny\n3\n\n\n\n\n\nQuit" | $RDSEED >>& ../${LOG}/_log_rdseed
    end
end
${RM} -f _sac_files_list
cd ../

foreach CHAN ($CHANS)
    ${LS}  ${DATA}/*${CHAN}*.SAC    >> ${DATA}/_sac_files_list
end

$TRIM_SAC_FILES i_master ${DATA}/_sac_files_list scr_dat_fil_list ${trim_flag}

########################
# Responses Lookup table
${ECHO} "Creating the responses lookup table ( >! ${LOG}/_log_resps_lookup_table )"

${LS} -1 ${DATA}/SAC_PZs*   > ${DATA}/pz_fil_list
$MAKE_RESP_TABLE ${DATA}/pz_fil_list i_master coeffs_rec_lut >! ${LOG}/_log_resps_lookup_table
