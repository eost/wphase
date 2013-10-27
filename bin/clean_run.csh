#!/bin/csh -f

set RM   = "/bin/rm"

set LIST = "DATA_org DATA GF LOG SYNTH SYNTH_traces gs_SYNTH"
foreach dir ($LIST)
   if ( -e $dir ) ${RM} -r $dir
end

set LIST = "       i_wpinversion coeffs_rec_lut dec_bp_dat_fil_list i_seedfiles _ref"
set LIST = "$LIST  o_wpinversion upd_dec_bp_dat_fil_list rot_dec_bp_dat_fil_list"
set LIST = "$LIST  scr_dat_fil_list syn_fil_list coeffs_rec_lut  dec_bp_dat_fil_list"
set LIST = "$LIST  i_wpinversion  rdseed.err_log wtraces.pdf Gd wp_pages_6t.pdf  wp_pages.pdf"
foreach file ($LIST)
   if ( -e $file ) ${RM} $file
end

${RM} -f fort.15.*  
${RM} -f o_* 
${RM} -f p_*    
${RM} -f WCMTSOLUTION*
${RM} -f rdseed.err_log.* 
${RM} -f CWP_* 
${RM} -f *_tmp* 
${RM} -f page_*
${RM} -f *fort.15*
${RM} -f grid_search_*
${RM} -f ts_*
${RM} -f xy_*
${RM} -f gs_*
${RM} -f dp_*
