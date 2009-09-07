#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

################################################

set CHANS   = "LHL LHT LHZ"

set LOG     = LOG
set DATA    = DATA

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

################################################
# deconvolution and filtering 
${ECHO} "Data deconvolution and filter...                   ( >! ${LOG}/_log_dec_filt )"

if -e ${LOG} then
    ${RM} -rf ${LOG}/_log_dec_filt
else
    ${MKDIR} ${LOG}
endif

$REC_DEC_FILT coeffs_rec_lut i_master scr_dat_fil_list dec_bp_dat_fil_list >> ${LOG}/_log_dec_filt

################################################
# Rotation of data horizontal components (from E/N to L/T)
${ECHO} "Rotation of horizontal components...               ( >! ${LOG}/_log_rot_data )"
$ROT_HORIZ_CMP dec_bp_dat_fil_list rot_dec_bp_dat_fil_list ${DATA} >! ${LOG}/_log_rot_data

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Synthetics convolution and filter...               ( >! ${LOG}/_log_synths_conv_filt    )"

${FIND} GF -name "*sac*" -exec ${RM} \{\} \+
${RM} -rf ${LOG}/_log_synths_conv_filt syn_fil_list
foreach CHAN ($CHANS)
    ${AWK} '{printf "%s.%s\n", $1, $2}' $gf_dir/STAT_LIST | ${SED} -e "s/.*/&.${CHAN}.SAC/" >> syn_fil_list
end

$SYN_CONV_FILT syn_fil_list l  -imas i_master -gfdir ${gf_dir} >! ${LOG}/_log_synths_conv_filt

###############################################
# Creating input file for inversion ... 
${CUT} -d' ' -f1  rot_dec_bp_dat_fil_list >! i_wpinversion

${RM} -f dat_fil_list
