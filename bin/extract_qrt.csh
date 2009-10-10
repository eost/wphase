#!/bin/csh -f

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
set miniseed_dir = `${GREP} SEED     i_tmp | ${CUT} -d: -f2`
set dataless_dir = `${GREP} DATALESS i_tmp | ${CUT} -d: -f2`
set date = `${HEAD} -1 ${cmtfile} | ${AWK} '{print $3"/"$4"/"$2}'`
set year = `${DATE} -d${date} +%Y`
set jday = `${DATE} -d${date} +%j`
set stalst = "station_list.${jday}"

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

${FIND} ${miniseed_dir} -name "*LH?.${year}.${jday}" > ${stalst}
#${LS} -1 ${miniseed_dir}*/*/*LH*.${jday} > ${stalst}

set Nchan = `${CAT} ${stalst} | ${WC}  -l`
set NETS  = `${CAT} ${stalst} | ${SED} 's/^.*\///' | ${CUT} -d. -f1`
set STATS = `${CUT} -d. -f2 ${stalst}`
set LOCID = `${CAT} ${stalst} | ${SED} 's/\.\./\.--\./' | ${CUT} -d. -f3`
set CHANS = `${CUT} -d. -f4 ${stalst}`
${ECHO} "====> total channe${LS} : ${Nchan}"

########################################

${RM} -f *PZs* rdseed.stations
set i = 1
${ECHO} "Extracting dataless information   (>> ${LOG}/_log_rdseed_pz   >>${LOG}/_log_rdseed_staloc)"
while ( $i <= $Nchan)
    set net = $NETS[$i]
    set sta = $STATS[$i]
    set lid = $LOCID[$i]
    set cha = $CHANS[$i]
    ${ECHO} "   ${net} ${sta} ${lid} ${cha}"

    # PZ files
    ${ECHO} -e "${dataless_dir}${net}.dataless\n\n\np\n${sta}\n${cha}\n${net}\n${lid}\nQuit\n" | $RDSEED >>& ${LOG}/_log_rdseed_pz
    set npz    = `${LS} -1 *PZs* | ${WC}  -l`
    if ( ${npz} == 0 ) then 
	@ i++
	continue
    else if ( ${npz} == 1 ) then
	${MV} *PZs* $DATA
    else
	set pzfile = `${LS} -1 *PZs* | ${SORT} | ${TAIL} -1`
	${MV} $pzfile $DATA
	${RM} -f *PZs*
    endif

    # coordinates
    ${ECHO} -e "${dataless_dir}${net}.dataless\n\n\nS\n${sta}\n${cha}\n${net}\n${lid}\nQuit\n" | $RDSEED >>& ${LOG}/_log_rdseed_staloc
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
	@ i++
	continue
    endif
    ${ECHO} ${net} ${sta} ${cha} ${lid} ${stla} ${stlo} ${stel} >> ${loclst}
    @ i++
end 

########################################
$TRIM_SAC_FILES_QRT i_master $loclst ${DATA} scr_dat_fil_list -u


########################
# Responses Lookup table
${ECHO} "Creating the responses lookup table ( >! ${LOG}/_log_resps_lookup_table )"

${LS} -1 ${DATA}/SAC_PZs*   > ${DATA}/pz_fil_list
$MAKE_RESP_TABLE ${DATA}/pz_fil_list i_master coeffs_rec_lut !> ${LOG}/_log_resps_lookup_table
