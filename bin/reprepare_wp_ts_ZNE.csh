#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

################################################

set CHANS   = "LHN LHE LHZ"
set LOG     = LOG

################################################

if (-e ts_i_master) then
        ${GREP} -v "^#" ts_i_master >! i_tmp
else
        ${ECHO} "Error: file  ts_i_master not available"
        exit
endif
set gf_dir   = "./ts_GF"
set tmp      = `${GREP} GFDIR   i_tmp`
if ! $status then
        set gf_dir   = `echo $tmp | ${HEAD} -1 | ${CUT} -d: -f2`
endif

${RM} -f i_tmp
${FIND} ${gf_dir} -name "*sac*" -exec ${RM} \{\} \+

if -e ${LOG} then
    ${RM} -rf ${LOG}/_ts_log_synths_conv_filt
else
    ${MKDIR} ${LOG}
endif

################################################
# Synthetics preparatio:convolution and filter #
${ECHO} "ts grid search : synthetics convolution and filter...( >! ${LOG}/_ts_log_synths_conv_filt    )"

if -e syn_fil_list then
    $SYN_CONV_FILT syn_fil_list l  -imas ts_i_master -gfdir ${gf_dir} >! ${LOG}/_ts_log_synths_conv_filt
else
    foreach CHAN ($CHANS)
	${AWK} '{printf "%s.%s\n", $1, $2}' $gf_dir/STAT_LIST | ${SED} -e "s/.*/&.${CHAN}.SAC/" >> syn_fil_list
    end
    $SYN_CONV_FILT syn_fil_list l -imas ts_i_master -gfdir ${gf_dir} >! ${LOG}/_ts_log_synths_conv_filt
endif

###############################################
# Creating input file for inversion ... 
#${CAT} o_wpinversion >! i_wpinversion
