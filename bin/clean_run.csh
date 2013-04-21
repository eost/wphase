#!/bin/csh -f

set RM = "/bin/rm -rf"
${RM} DATA_org DATA GF LOG SYNTH SYNTH_traces gs_SYNTH
${RM} i_wpinversion coeffs_rec_lut dec_bp_dat_fil_list i_seedfiles _ref_o_wpinversion rdseed.err_log.*
${RM} upd_dec_bp_dat_fil_list rot_dec_bp_dat_fil_list scr_dat_fil_list syn_fil_list
${RM} CWP_*
${RM} wtraces.pdf
${RM} *_tmp* 
${RM} o_* 
${RM} p_* page_* 
${RM} *fort.15*
${RM} WCMTSOLUTION*
${RM} grid_search_*
${RM} ts_*
${RM} xy_*
${RM} gs_*
${RM} dp_*
${RM} Gd
