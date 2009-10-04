#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

set my_argv = ($ARGV)
if ($#my_argv < 1) then
    set pre="xy_"
else if ($#my_argv == 1) then
    set pre = $my_argv[1]
endif

################################################

set CHANS   = "LHN LHE LHZ"
set LOG     = LOG
set DATA    = ${pre}DATA

########################################

if (-e ${pre}i_master) then
        ${GREP} -v "^#" ${pre}i_master >! i_tmp
else
        ${ECHO} "Error: file  ${pre}i_master not available"
        exit
endif

set CMTFILE  = `${GREP} CMTFILE i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
set gf_dir   = "./${pre}GF"
set tmp      = `${GREP} GFDIR   i_tmp`
if ! $status then
        set gf_dir   = `echo $tmp | ${HEAD} -1 | ${CUT} -d: -f2`
endif

${RM} -f i_tmp
if -e ${LOG} then
    ${RM} -rf ${LOG}/_${pre}log_rot_data ${LOG}/_${pre}log_synths_conv_filt 
else
    ${MKDIR} ${LOG}
endif

################################################
# Rotation of data horizontal components (from E/N to L/T)
${ECHO} "Rotation of horizontal components...               ( >! ${LOG}/_log_rot_data )"
$CHANGE_SAC_HEADERS dec_bp_dat_fil_list upd_dec_bp_dat_fil_list ${DATA} -icmtf ${CMTFILE} >! ${LOG}/_log_rot_data

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "Centroid Grid Search : synthetics convolution and filter...( >! ${LOG}/_${pre}log_synths_conv_filt    )"

${FIND} ${gf_dir} -name "*sac*" -exec ${RM} \{\} \+

if -e syn_fil_list then
    $SYN_CONV_FILT syn_fil_list l -imas ${pre}i_master -gfdir ${gf_dir} >! ${LOG}/_${pre}log_synths_conv_filt
    #$SYN_CONV_FILT ${pre}i_master syn_fil_list l rr tt pp rt rp tp >! ${LOG}/_${pre}log_synths_conv_filt
else
    foreach CHAN ($CHANS)
	${AWK} '{printf "%s.%s\n", $1, $2}' $gf_dir/STAT_LIST | ${SED} -e "s/.*/&.${CHAN}.SAC/" >> syn_fil_list
    end
    $SYN_CONV_FILT syn_fil_list l -imas ${pre}i_master -gfdir ${gf_dir} >! ${LOG}/_${pre}log_synths_conv_filt
    #$SYN_CONV_FILT ${pre}i_master syn_fil_list l rr tt pp rt rp tp >! ${LOG}/_${pre}log_synths_conv_filt
endif

###############################################
# Creating input file for inversion ...
#${CAT} o_wpinversion | ${SED} s/DATA/${DATA}/ >! i_wpinversion
