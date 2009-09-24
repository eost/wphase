#!/bin/csh -f

########################################
source $WPHASE_HOME/bin/WP_HEADER.CSH

set DATA = DATA
set LOG  = LOG

########################################
if (-e i_master) then
        ${GREP} -v "^#" i_master | ${SED} -e 's/ *:/:/' >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif

${RM} -rf ${LOG}
${MKDIR} ${LOG}

########################
# Responses Lookup table
${ECHO} "Creating the responses lookup table ( >! ${LOG}/_log_resps_lookup_table )"

${LS} -1 ${DATA}/SAC_PZs*   > ${DATA}/pz_fil_list
$MAKE_RESP_TABLE ${DATA}/pz_fil_list i_master coeffs_rec_lut !> ${LOG}/_log_resps_lookup_table
