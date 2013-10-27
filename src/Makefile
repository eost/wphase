############################################################################
#
#	              W phase source inversion package 	            
#                               -------------
#
#        Main authors: Zacharie Duputel, Luis Rivera and Hiroo Kanamori
#                      
# (c) California Institute of Technology and Universite de Strasbourg / CNRS 
#                                  April 2013
#
#    Neither the name of the California Institute of Technology (Caltech) 
#    nor the names of its contributors may be used to endorse or promote 
#    products derived from this software without specific prior written 
#    permission
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################


# Compilers
CC          =  gcc          # C compiler
FC	    =  gfortran     # fortran compiler
			    #
# Some parameters
datalen     = 5000  # max. nb of samples
fnamelen    = 128   # max. nb of chars in file names
lineslen    = 256   # max. nb of chars per line in text files
idlen       = 16    # len of station id

OPTFLAG  = -fopenmp -D__GFS_01D__=1

CFLAGS   = -O2 -Wall -D__LEN_SIG__=$(datalen) \
	   -D__FSIZE__=$(fnamelen) -D__IDSIZE__=$(idlen)\
	   -D__LSIZE__=$(lineslen) $(OPTFLAG)

LDFLAGS  = -O2 -fsecond-underscore $(OPTFLAG)

FFLAGS   = -O2 -x f77-cpp-input -ffixed-line-length-132 -fsecond-underscore \
		  -D__LEN_SIG__=$(length) -std=legacy $(OPTFLAG)

##############################

BIN         = ../bin

export CC 
export FC

export CFLAGS
export LDFLAGS

EXE = 	$(BIN)/WP_HEADER.CSH		\
	$(BIN)/rec_dec_filt		\
	$(BIN)/syn_conv_filt		\
	$(BIN)/make_resp_lookup_table 	\
	$(BIN)/trim_sac_files		\
	$(BIN)/wpinversion		\
	$(BIN)/synth_v6			\
	$(BIN)/saclst                  	\
	$(BIN)/fast_synth 		\
	$(BIN)/fast_synth_only_Z  	\
	$(BIN)/fast_synth_rot 		\
	$(BIN)/prep_kernels_only_Z 	\
	$(BIN)/prep_kernels 		\
	$(BIN)/wpinversion_gs 		\
	$(BIN)/cmtascii			\
	$(BIN)/decim_one_sac_file_to_1sps

OBJS =  read_stats.o     \
	get_prefix.o     \
	rotate_cmt.o     \
	rotate_traces.o  \
	rotate_2_ortho_traces.o  \
	yyyymmdd2jjj.o   \
	save_sac.o       \
	get_length.o     \
	get_depths.o     \
	fast_synth_sub.o \
	sum_up.o         \
	tapering.o       \
	travel_times.o 	 \
	read_i_files.o 	 \
	distaz.o 	 \
        proto_alloc.o  	 \
	rwsacs.o	 \
	rwtextfiles.o  

all : $(EXE)

# Scripts header
$(BIN)/WP_HEADER.CSH:configure.csh
	@ ./configure.csh

###########
# LINKING

# optional routines:
$(BIN)/cmtascii:cmtascii.o proto_alloc.o rwtextfiles.o jacobi.o charplot.o read_i_files.o
	$(CC) $(LDFLAGS) cmtascii.o proto_alloc.o rwtextfiles.o jacobi.o charplot.o \
        read_i_files.o -o $(BIN)/cmtascii -lm

# Inversion
$(BIN)/wpinversion:wpinversion.o wpinversion_sub.o jacobi.o charplot.o \
	lsqenp2.o syn_conv_sub.o butterworth.o $(OBJS)
	$(FC) $(LDFLAGS) wpinversion.o wpinversion_sub.o jacobi.o charplot.o \
	lsqenp2.o syn_conv_sub.o butterworth.o $(OBJS) -o $(BIN)/wpinversion -lm

# Green's functions
$(BIN)/fast_synth: fast_synth.o $(OBJS)
	$(CC) ${LDFLAGS} -o $(BIN)/fast_synth  fast_synth.o $(OBJS) -lm 

$(BIN)/fast_synth_rot: fast_synth_rot.o $(OBJS)
	$(CC) ${LDFLAGS} -o $(BIN)/fast_synth_rot  fast_synth_rot.o $(OBJS) -lm 

$(BIN)/fast_synth_only_Z: fast_synth_only_Z.o $(OBJS)
	$(CC) ${LDFLAGS} -o $(BIN)/fast_synth_only_Z fast_synth_only_Z.o $(OBJS) -lm 

fast_synth.o: fast_synth.c
	${CC} ${CFLAGS} -c fast_synth.c  -o fast_synth.o 

$(BIN)/prep_kernels: prep_kernels.o butterworth.o sys_cmd.o syn_conv_sub.o $(OBJS)
	$(CC) ${LDFLAGS} -o $(BIN)/prep_kernels prep_kernels.o butterworth.o sys_cmd.o syn_conv_sub.o $(OBJS) -lm 

$(BIN)/prep_kernels_only_Z: prep_kernels_only_Z.o butterworth.o sys_cmd.o syn_conv_sub.o $(OBJS)
	$(CC) ${LDFLAGS} -o $(BIN)/prep_kernels_only_Z prep_kernels_only_Z.o butterworth.o sys_cmd.o syn_conv_sub.o $(OBJS) -lm 


$(BIN)/wpinversion_gs: wpinversion_gs.o wpinversion_sub.o \
	jacobi.o charplot.o butterworth.o syn_conv_sub.o $(OBJS) lsqenp2.o
	$(FC) $(LDFLAGS) wpinversion_gs.o wpinversion_sub.o jacobi.o \
	charplot.o butterworth.o syn_conv_sub.o $(OBJS) lsqenp2.o -o $(BIN)/wpinversion_gs -lm

# MT rotation
$(BIN)/rot_horiz_cmp: rot_horiz_cmp.o proto_alloc.o rwsacs.o \
	rwtextfiles.o read_i_files.o distaz.o
	$(CC) $(LDFLAGS) rot_horiz_cmp.o proto_alloc.o rwsacs.o \
	rwtextfiles.o read_i_files.o distaz.o -o $(BIN)/rot_horiz_cmp -lm 

# Filtering and deconvolution
$(BIN)/rec_dec_filt: rec_dec_filt.o butterworth.o      \
	rwsacs.o proto_alloc.o rwtextfiles.o read_i_files.o sort_tree.o 
	$(CC) ${LDFLAGS} rec_dec_filt.o butterworth.o  sort_tree.o \
	rwsacs.o proto_alloc.o rwtextfiles.o read_i_files.o -o  $(BIN)/rec_dec_filt -lm

$(BIN)/syn_conv_filt: syn_conv_filt.o butterworth.o rwsacs.o \
	proto_alloc.o rwtextfiles.o read_i_files.o syn_conv_sub.o
	$(CC) ${LDFLAGS} syn_conv_filt.o butterworth.o syn_conv_sub.o \
	rwsacs.o proto_alloc.o rwtextfiles.o read_i_files.o -o  $(BIN)/syn_conv_filt -lm


$(BIN)/synth_v6: synth_v6.o butterworth.o syn_conv_sub.o $(OBJS)
	$(CC) ${LDFLAGS} synth_v6.o butterworth.o syn_conv_sub.o $(OBJS) \
	-o  $(BIN)/synth_v6 -lm

$(BIN)/make_resp_lookup_table: make_resp_lookup_table.o lsqenp.o \
	proto_alloc.o rwtextfiles.o read_i_files.o
	$(FC) $(LDFLAGS) make_resp_lookup_table.o lsqenp.o proto_alloc.o rwtextfiles.o \
	read_i_files.o -o $(BIN)/make_resp_lookup_table


# Extract and data screening
$(BIN)/trim_sac_files: trim_sac_files.o distaz.o sort_tree.o rwsacs.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o
	$(CC) $(LDFLAGS) trim_sac_files.o distaz.o sort_tree.o rwsacs.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o -o $(BIN)/trim_sac_files -lm

${BIN}/saclst: saclst.o
	$(CC) ${CFLAGS} -o ${BIN}/saclst saclst.o

# Decimate
$(BIN)/decim_one_sac_file_to_1sps: decim_one_sac_file_to_1sps.o decimate.o proto_alloc.o rwsacs.o
	$(CC) $(LDFLAGS) decim_one_sac_file_to_1sps.o decimate.o proto_alloc.o rwsacs.o -o $(BIN)/decim_one_sac_file_to_1sps

########################
# COMPILE AND ASSEMBLE

lsqenp.o: lsqenp.f
	$(FC) $(FFLAGS) -c lsqenp.f 

lsqenp2.o: lsqenp2.f
	$(FC) $(FFLAGS) -c lsqenp2.f 


%.o: %.c
	${CC} ${CFLAGS} -c $<  -o $@ 


# Cleaning
clean:
	\rm -f *.o

cleanexe:
	\rm -f $(EXE)
