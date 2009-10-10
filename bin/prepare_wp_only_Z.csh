#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

################################################

set CHAN    = LHZ
set gf_dir  = GF
set LOG     = LOG

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
set CMTFILE  = `${GREP} CMTFILE i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
${RM} -f i_tmp

################################################
# deconvolution and filtering 
${ECHO} "Data deconvolution and filter...                   ( >! ${LOG}/_log_dec_filt )"

if -e ${LOG} then
    ${RM} -rf LOG/_log_dec_filt
else
    ${MKDIR} ${LOG}
endif

#${CUT} -d' ' -f1 screened_data_file_list | ${AWK} '{print $1}' >! dat_fil_list
$REC_DEC_FILT coeffs_rec_lut i_master scr_dat_fil_list dec_bp_dat_fil_list >> ${LOG}/_log_dec_filt

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Synthetics convolution and filter...               ( >! ${LOG}/_log_synths_conv_filt    )"

${FIND} ${gf_dir} -name "*sac*" -exec ${RM} \{\} \+
${RM} -rf ${LOG}/_log_synths_conv_filt
${AWK} '{printf "%s.%s\n", $1, $2}' ${gf_dir}/STAT_LIST | ${SED} -e "s/.*/&.${CHAN}.SAC/" > syn_fil_list

$SYN_CONV_FILT syn_fil_list l  -imas i_master -gfdir ${gf_dir} >! ${LOG}/_log_synths_conv_filt

################################################
# Creating input file for inversion ... 
${CUT} -d' ' -f1  dec_bp_dat_fil_list  > i_wpinversion

${RM} -f dat_fil_list
