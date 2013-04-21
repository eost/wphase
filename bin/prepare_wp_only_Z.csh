#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH
set LOG     = LOG

###############################################
if -e ${LOG} then
    ${RM} -rf LOG/_log_dec_filt
else
    ${MKDIR} ${LOG}
endif

################################################
# deconvolution and filtering 
${ECHO} "Data deconvolution and filter...                   ( >! ${LOG}/_log_dec_filt )"
$REC_DEC_FILT coeffs_rec_lut i_master scr_dat_fil_list dec_bp_dat_fil_list >> ${LOG}/_log_dec_filt

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Synthetics convolution and filter... "
${PREP_KERNELS_only_Z} scr_dat_fil_list l
if $status exit(1)

################################################
# Creating input file for inversion ... 
${CUT} -d' ' -f1  dec_bp_dat_fil_list  > i_wpinversion
