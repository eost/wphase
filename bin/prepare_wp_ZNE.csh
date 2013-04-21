#!/bin/csh -f
# Command lines examples:
# prepare_wp.csh 
# prepare_wp.csh Z
# prepare_wp.csh -a
#

source $WPHASE_HOME/bin/WP_HEADER.CSH
set LOG      = LOG
set DATA     = DATA
set DATA_org = DATA_org
set CMPS     = "Z N E"

# command line arguments 
set trim_flag = '-u'
set my_argv = ($ARGV)
if ($#my_argv > 0) then
    set CMPS = $my_argv
    if ( $CMPS == '-a' ) then
	set CMPS = "Z N E"
	set trim_flag = "-a"
    endif
endif


if ( ! -e ${DATA_org} ) then
   ${ECHO} "Missing ${DATA_org}"
   exit 1
endif
 
if ( -e ${DATA} ) $RM -rf ${DATA}
${MKDIR} ${DATA}
${MKDIR} -p ${LOG}
###############################################
foreach CMP ($CMPS)
    ${CP} -f ${DATA_org}/SAC_PZs_*_*_??${CMP}_* ${DATA}
    ${SACLST} kcmpnm f ${DATA_org}/*.SAC | ${EXPAND} | ${GREP} " ..${CMP}" | ${CUT} -d' ' -f1 | \
	$XARGS -n 1 -I {} ${DECIMATE} {} ${DATA} >! ${LOG}/_log_decimate
end
${LS} ${DATA}/*.SAC >! ${DATA}/_sac_files_list

$TRIM_SAC_FILES i_master ${DATA}/_sac_files_list scr_dat_fil_list $trim_flag # Use of "-u" will allow only one network per-channel

########################
# Responses Lookup table
${ECHO} "Creating the responses lookup table ( >! ${LOG}/_log_resps_lookup_table )"
${LS} -1 ${DATA}/SAC_PZs* >! ${DATA}/pz_fil_list
$MAKE_RESP_TABLE ${DATA}/pz_fil_list i_master coeffs_rec_lut >! ${LOG}/_log_resps_lookup_table

################################################
# deconvolution and filtering                  #
${ECHO} "Data deconvolution and filter...                   ( >! ${LOG}/_log_dec_filt )"
$REC_DEC_FILT coeffs_rec_lut i_master scr_dat_fil_list dec_bp_dat_fil_list >! ${LOG}/_log_dec_filt

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Synthetics convolution and filter...  "
${PREP_KERNELS_ZNE} scr_dat_fil_list l
if $status exit(1)

################################################
# Creating input file for inversion ...        #
${CUT} -d' ' -f1  dec_bp_dat_fil_list >! i_wpinversion
