#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH
set LOG     = LOG
set DATA    = DATA

###############################################
if (-e i_master) then
        ${GREP} -v "^#" i_master >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif
set CMTFILE  = `${GREP} CMTFILE i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
${RM} -f i_tmp
if -e ${LOG} then
    ${RM} -rf ${LOG}/_log_dec_filt
else
    ${MKDIR} ${LOG}
endif

################################################
# deconvolution and filtering 
${ECHO} "Data deconvolution and filter...                   ( >! ${LOG}/_log_dec_filt )"
$REC_DEC_FILT coeffs_rec_lut i_master scr_dat_fil_list dec_bp_dat_fil_list >> ${LOG}/_log_dec_filt

################################################
# Rotation of data horizontal components (from E/N to L/T)
${ECHO} "Rotation of horizontal components...               ( >! ${LOG}/_log_rot_data )"
$ROT_HORIZ_CMP dec_bp_dat_fil_list rot_dec_bp_dat_fil_list ${DATA} -icmtf ${CMTFILE} >! ${LOG}/_log_rot_data

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Synthetics convolution and filter...         "
${PREP_KERNELS_LTZ} scr_dat_fil_list l
if $status exit(1)

###############################################
# Creating input file for inversion ... 
${CUT} -d' ' -f1  rot_dec_bp_dat_fil_list >! i_wpinversion
${RM} -f dat_fil_list
