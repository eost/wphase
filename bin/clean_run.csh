#!/bin/csh -f

set RM = "/bin/rm -rf"
${RM} DATA GF LOG SYNTH
${RM} i_wpinversion coeffs_rec_lut dec_bp_dat_fil_list i_seedfiles
${RM} rot_dec_bp_dat_fil_list scr_dat_fil_list syn_fil_list
${RM} CWP_* 
${RM} *_tmp* 
${RM} o_* 
${RM} p_* page_* 
${RM} fort.15*
${RM} WCMTSOLUTION*
${RM} grid_search_*
${RM} ts_*
${RM} xy_*

