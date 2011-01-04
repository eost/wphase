#!/bin/csh -f
#
# W phase package - Extract data from miniseed volumes
#
# Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#

##############################################
# Extract for quasi-real time implementation
# using slarchive
########################################

source $WPHASE_HOME/bin/WP_HEADER.CSH

set DATA = DATA
set LOG  = LOG

set loclst = loc_list

########################################
if (-e i_master) then
        ${GREP} -v "^#" i_master | ${SED} -e 's/ *:/:/' >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit
endif
set cmtfile      = `${GREP} CMTFILE  i_tmp | ${CUT} -d: -f2`
set miniseed_dir = `${GREP} SEED     i_tmp | ${CUT} -d: -f2 | ${SED} 's/\/$//'`
set dataless_dir = `${GREP} DATALESS i_tmp | ${CUT} -d: -f2 | ${SED} 's/\/$//'`
set date   = `${HEAD} -1 ${cmtfile} | ${CUT} -c6-15 | ${AWK} '{print $2"/"$3"/"$1}'`
set year   = `${DATE} -d${date} +%Y`
set jday   = `${DATE} -d${date} +%j`
set epdate = `${DATE} -d${date} -u +%s`
set stalst = "miniseed_list"

echo $year $jday
${RM} -f i_tmp

########################################
${RM} -f ${stalst}
${RM} -f  ${loclst}

${RM} -rf ${DATA}
${RM} -rf ${LOG}
${MKDIR} ${DATA}
${MKDIR} ${LOG}
########################################
${ECHO} "Creating channel list"

${FIND} ${miniseed_dir}/ -name "*BH?.*${year}.*${jday}*" > ${stalst} ### MSEED FILE NAMES
#${LS} -1 ${miniseed_dir}*/*/*LH*.${jday} > ${stalst}

set mseeds = `${CAT} ${stalst}`    
set Nchan = `${CAT} ${stalst} | ${WC}  -l`
${ECHO} "====> total channel : ${Nchan}"

########################################

${RM} -f _tmp_o_ids
${TOUCH} _tmp_o_ids
if (-e o_ids) then
    ${GREP} '0$' o_ids | $AWK '{print $1}' >! _tmp_o_ids
endif

${RM} -f *PZs* rdseed.stations
${ECHO} "Extracting dataless information   (>> ${LOG}/_log_rdseed_pz   >>${LOG}/_log_rdseed_staloc)"
foreach mseed (${mseeds})
    set channel = `${MSEEDLST} $mseed`
    set sta = `${ECHO} $channel | ${AWK} '{print $1}'`
    set net = `${ECHO} $channel | ${AWK} '{print $2}'`    
    set lid = `${ECHO} $channel | ${AWK} '{print $3}'`
    set cha = `${ECHO} $channel | ${AWK} '{print $4}'`
    set id = ${net}.${sta}.${lid}.${cha}    
    set buf = `$GREP ${id} _tmp_o_ids`
    if ! $status then
	echo "${id} rejected (status 0 in o_ids)" 
	@ i++
	continue
    endif
    ${ECHO} "   ${net} ${sta} ${lid} ${cha}"
 
    # PZ files
    ${ECHO} -e "${dataless_dir}/${net}.dataless\n\n\np\n${sta}\n${cha}\n${net}\n${lid}\nQuit\n" | $RDSEED >>& ${LOG}/_log_rdseed_pz
    set flag = 1
    set npz    = `${LS} -1 *PZs* | ${WC}  -l`
    if ( ${npz} == 1 ) then
	${MV} *PZs* $DATA
	set flag = 0
    else
	set pzfiles = `${LS} -1 *PZs*`
	foreach pzfile ($pzfiles)
	    set pzyear1 = `$ECHO $pzfile | ${CUT} -d'.' -f1 | ${SED} 's/^.*_//'`
	    set pzjday1 = `$ECHO $pzfile | ${CUT} -d'.' -f2`
	    set pzyear2 = `$ECHO $pzfile | ${CUT} -d'.' -f6 | ${CUT} -d'_' -f2`
	    set pzjday2 = `$ECHO $pzfile | ${CUT} -d'.' -f7`
	    @ pzjday1--
	    @ pzjday2--
	    set pzepoch1 = `$DATE -d"${pzyear1}-01-01 +${pzjday1}days" -u +'%s'`
	    set pzepoch2 = `$DATE -d"${pzyear2}-01-01 +${pzjday2}days" -u +'%s'`
	    if ( $pzepoch1 <= ${epdate} && $pzepoch2 >= ${epdate} ) then
		set flag = 0
		${MV} $pzfile $DATA
		break
	    endif
	end
	${RM} -f *PZs*
    endif
    if ( ${flag} == 1) then 
	${ECHO} "WARNING : No matching PZFILE for ${net} ${sta} ${lid} ${cha}"
	continue
    endif

    # Coordinates
    ${ECHO} -e "${dataless_dir}/${net}.dataless\n\n\nS\n${sta}\n${cha}\n${net}\n${lid}\nQuit\n" | $RDSEED >>& ${LOG}/_log_rdseed_staloc
    set N = `${CAT} rdseed.stations | ${GREP} ^${sta} | ${WC}  -l`
    if ( 1 < ${N} ) then
	${CAT} rdseed.stations | ${GREP} 2500,365,23:59:59.9999 > ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV
    else
	${CAT} rdseed.stations > ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV
    endif
    set stanm = `${HEAD} -1 ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV | ${AWK} '{print $1}'`
    set netnm = `${HEAD} -1 ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV | ${AWK} '{print $2}'`
    set stla  = `${HEAD} -1 ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV | ${AWK} '{print $3}'`
    set stlo  = `${HEAD} -1 ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV | ${AWK} '{print $4}'`
    set stel  = `${HEAD} -1 ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV | ${AWK} '{print $5}'`
    if ( ${stanm} != ${sta} || ${netnm} != ${net}) then
	${ECHO} "WARNING ${net} ${sta} ${lid} ${cha} rejected (LOC file incorrect)"
	${ECHO} "WARNING ${net} ${sta} ${lid} ${cha} rejected (LOC file incorrect)" >> ${LOG}/_log_rdseed_staloc
	${RM} -f ${DATA}/SAC_LOC_${net}_${sta}_${cha}_${lid}_LAT_LON_ELEV
	${RM} -f ${DATA}/SAC_PZs_${net}_${sta}_${cha}_${lid}_*
	continue
    endif
    ${ECHO} ${net} ${sta} ${cha} ${lid} ${stla} ${stlo} ${stel} >> ${loclst}
end 

########################################

$TRIM_SAC_FILES_QRT i_master $loclst $stalst ${DATA} scr_dat_fil_list -u


########################
# Responses Lookup table
${ECHO} "Creating the responses lookup table ( >! ${LOG}/_log_resps_lookup_table )"

${LS} -1 ${DATA}/SAC_PZs*   > ${DATA}/pz_fil_list
$MAKE_RESP_TABLE ${DATA}/pz_fil_list i_master coeffs_rec_lut >! ${LOG}/_log_resps_lookup_table
