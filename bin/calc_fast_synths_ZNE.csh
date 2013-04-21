#!/bin/csh -f

source $WPHASE_HOME/bin/WP_HEADER.CSH

###################################
if (-e i_master) then
	${GREP} -v "^#" i_master >! i_tmp
else
        ${ECHO} "Error: file  i_master not available"
        exit 1
endif
########################################

set MOM_POW  = 28
set EVNAME   = `${GREP} EVNAME  i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
set CMTFILE  = `${GREP} CMTFILE i_tmp | ${HEAD} -1 | ${CUT} -d: -f2`
set gf_dir   = "./GF"
set tmp      = `${GREP} GFDIR   i_tmp`
if ! $status then
	set gf_dir   = `echo $tmp | ${HEAD} -1 | ${CUT} -d: -f2`
endif

${RM} -f i_tmp
${RM} -rf $gf_dir
${MKDIR} $gf_dir

# setting up the list of stations for which to calculate synthetics.
${AWK} '{printf "%-4s %-2s %9.4f %9.4f\n",$2,$3,$5,$6}' scr_dat_fil_list | ${SORT} --unique >! $gf_dir/STAT_LIST

# Elementary moment tensors
${HEAD} -2 $CMTFILE                                      >! $gf_dir/CMT_tmp
${GREP} -i shift     $CMTFILE | ${SED} -e "s/[0-9]/0/g"     >> $gf_dir/CMT_tmp
${GREP} -i duration  $CMTFILE | ${SED} -e "s/[0-9]/0/g"     >> $gf_dir/CMT_tmp
${GREP} -i latitude  $CMTFILE                            >> $gf_dir/CMT_tmp
${GREP} -i longitude $CMTFILE                            >> $gf_dir/CMT_tmp
${GREP} -i depth     $CMTFILE                            >> $gf_dir/CMT_tmp

set cmps  = "rr tt pp rt rp tp"
set cmps2 = "rr tt pp rt rp tp"
foreach cmp ($cmps)
	${MKDIR} $gf_dir/gf_$cmp
	cd $gf_dir/gf_$cmp
        ${CP} ../CMT_tmp CMTSOLUTION_$cmp
        foreach cmp2 ($cmps2)
                if ($cmp == $cmp2) then
                        ${ECHO} "M${cmp2}:       1.000000e+$MOM_POW" >> CMTSOLUTION_$cmp
                else
                        ${ECHO} "M${cmp2}:       0.000000e+$MOM_POW" >> CMTSOLUTION_$cmp
                endif
        end
        $FAST_SYNTH_ROT CMTSOLUTION_$cmp ../STAT_LIST
	cd ../..
end
