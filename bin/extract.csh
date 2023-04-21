#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

########################################
set DATA = DATA_org
set LOG  = LOG
########################################

if (! -e i_master) then
        ${ECHO} "Error: file  i_master not available"
        exit 1
endif

########################################
if ( -e ${DATA} ) ${RM} -rf ${DATA} 
if ( -e ${LOG} )  ${RM} -rf ${LOG}
${MKDIR} ${DATA} ${LOG}
########################################

# Extract

${ECHO} "Extracting data from seed           ( >! ${LOG}/_log_rdseed )"
${ECHO} '' >! ${LOG}/_log_rdseed

set seeds  = `${GREP} SEED i_master | ${GREP} -v "^#" | ${CUT} -d: -f2`
foreach seed ($seeds)
    if (! -e ${seed}) then
        ${ECHO} "Error: seed file $seed not available"
        exit 1
    endif
    ${RDSEED} -dp -q ${DATA} -f ${seed} >>& ${LOG}/_log_rdseed
end
