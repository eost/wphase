# Compilers
CC          =  cc           # C compiler
FC	    =  gfortran-4.3 # fortran compiler

# Some parameters
datalen     = 5000  # max. nb of samples
fnamelen    = 128   # max. nb of chars in file names
lineslen    = 256   # max. nb of chars per line in text files
idlen       = 16    # len of station id
safetydist  = 2.   
safetydelay = 10.

CFLAGS   = -g -Wall \
		  -D__LEN_SIG__=$(datalen) \
		  -D__FSIZE__=$(fnamelen) \
		  -D__IDSIZE__=$(idlen)\
	   	  -D__LSIZE__=$(lineslen) \
	 	  -D__SAFETY_DELAY__=$(safetydelay) \
		  -D__SAFETY_DIST__=$(safetydist)

LDFLAGS   = -g -fsecond-underscore

FFLAGS    = -x f77-cpp-input -g -ffixed-line-length-132 -fsecond-underscore \
		  -D__LEN_SIG__=$(length)

# Libmseed dependencies : type 'make with_miniseed' if miniseed reading is required
CINCLUDES   = -I/opt/libmseed
LDFLAGSSUP  = -L /opt/libmseed
LDLIBS      = -lmseed
##############################

BIN         = ../bin
SYNT_DIR    = FAST_SYNTH 

export CC 
export FC

export CFLAGS
export LDFLAGS

EXE = $(BIN)/rec_dec_filt_v1\
      $(BIN)/syn_conv_filt_v1\
      $(BIN)/make_resp_lookup_table\
      $(BIN)/trim_sac_files\
      $(BIN)/FAST\
      $(BIN)/wpinversion\
      $(BIN)/synth_v5\
      $(BIN)/rot_horiz_cmp\
      $(BIN)/saclst

all : $(EXE)

# enable libmseed dependencies
with_miniseed:  $(EXE) $(BIN)/trim_sac_files_qrt

# Inversion
$(BIN)/wpinversion:wpinversion.o rwsacs.o proto_alloc.o rwtextfiles.o \
	travel_times.o jacobi.o charplot.o read_i_files.o
	$(CC) $(LDFLAGS) wpinversion.o rwsacs.o proto_alloc.o rwtextfiles.o \
	travel_times.o jacobi.o charplot.o read_i_files.o -o $(BIN)/wpinversion -lm

wpinversion.o:wpinversion.c
	$(CC) $(CFLAGS) -c wpinversion.c

jacobi.o:jacobi.c
	$(CC) $(CFLAGS) -c jacobi.c

charplot.o:charplot.c
	$(CC) $(CFLAGS) -c charplot.c


# Green's functions

$(BIN)/FAST:
	@(cd $(SYNT_DIR) && $(MAKE))


# MT rotation

$(BIN)/rot_horiz_cmp: rot_horiz_cmp.o proto_alloc.o rwsacs.o \
	rwtextfiles.o read_i_files.o distaz.o 
	$(CC) $(LDFLAGS) rot_horiz_cmp.o proto_alloc.o rwsacs.o \
	rwtextfiles.o read_i_files.o distaz.o -o $(BIN)/rot_horiz_cmp -lm

rot_horiz_cmp.o:rot_horiz_cmp.c
	$(CC) $(CFLAGS) -c rot_horiz_cmp.c


# Recursive Filtering

$(BIN)/rec_dec_filt_v1: rec_dec_filt_v1.o butterworth_sin.o      \
	rwsacs.o proto_alloc.o rwtextfiles.o complex.o read_i_files.o sort_tree.o 
	${CC} ${LDFLAGS} rec_dec_filt_v1.o butterworth_sin.o  sort_tree.o \
	rwsacs.o proto_alloc.o rwtextfiles.o complex.o read_i_files.o -o  $(BIN)/rec_dec_filt_v1 -lm

$(BIN)/syn_conv_filt_v1: syn_conv_filt_v1.o butterworth_sin.o rwsacs.o \
	proto_alloc.o rwtextfiles.o complex.o read_i_files.o syn_conv_sub.o
	${CC} ${LDFLAGS} syn_conv_filt_v1.o butterworth_sin.o syn_conv_sub.o \
	rwsacs.o proto_alloc.o rwtextfiles.o complex.o read_i_files.o -o  $(BIN)/syn_conv_filt_v1 -lm


$(BIN)/synth_v5: synth_v5.o butterworth_sin.o rwsacs.o \
	proto_alloc.o rwtextfiles.o complex.o read_i_files.o syn_conv_sub.o
	${CC} ${LDFLAGS} synth_v5.o butterworth_sin.o syn_conv_sub.o \
	rwsacs.o proto_alloc.o rwtextfiles.o complex.o read_i_files.o -o  $(BIN)/synth_v5 -lm


rec_dec_filt_v1.o: rec_dec_filt_v1.c
	${CC} ${CFLAGS} -c rec_dec_filt_v1.c 

syn_conv_filt_v1.o: syn_conv_filt_v1.c
	${CC} ${CFLAGS} -c syn_conv_filt_v1.c 

syn_conv_sub.o: syn_conv_sub.c
	${CC} ${CFLAGS} -c syn_conv_sub.c 

synth_v5.o: synth_v5.c
	${CC} ${CFLAGS} -c synth_v5.c 

butterworth_sin.o: butterworth_sin.c
	${CC} ${CFLAGS} -c butterworth_sin.c


# Deconvolution parameters

$(BIN)/make_resp_lookup_table: make_resp_lookup_table.o lsqenp.o \
	proto_alloc.o rwtextfiles.o complex.o read_i_files.o
	$(FC) $(LDFLAGS) make_resp_lookup_table.o lsqenp.o proto_alloc.o rwtextfiles.o \
	complex.o read_i_files.o -o $(BIN)/make_resp_lookup_table

deconvolv.o: make_resp_lookup_table.c
	$(CC) $(CFLAGS) -c make_resp_lookup_table.c

lsqenp.o: lsqenp.f
	$(FC) $(FFLAGS) -c lsqenp.f 


# Extract and data screening

#(SEED)
$(BIN)/trim_sac_files: trim_sac_files.o distaz.o sort_tree.o rwsacs.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o
	$(CC) $(LDFLAGS) trim_sac_files.o distaz.o sort_tree.o rwsacs.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o -o $(BIN)/trim_sac_files -lm

trim_sac_files.o:trim_sac_files.c
	$(CC) $(CFLAGS) -c trim_sac_files.c

#(MINISEED)
$(BIN)/trim_sac_files_qrt: trim_sac_files_qrt.o distaz.o sort_tree.o rwsacs.o rmseed.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o 
	$(CC) $(LDFLAGS) trim_sac_files_qrt.o distaz.o sort_tree.o rwsacs.o rmseed.o \
	proto_alloc.o rwtextfiles.o travel_times.o read_i_files.o -o $(BIN)/trim_sac_files_qrt $(LDFLAGSSUP) $(LDLIBS) -lm 

trim_sac_files_qrt.o:trim_sac_files_qrt.c
	$(CC) $(CFLAGS) $(CINCLUDES) -c trim_sac_files_qrt.c

rmseed.o:rmseed.c
	$(CC) $(CFLAGS) $(CINCLUDES) -c rmseed.c

#(bin-tree sorting)
sort_tree.o:sort_tree.c
	$(CC) $(CFLAGS) -c sort_tree.c

# Shared routines

read_i_files.o:read_i_files.c
	$(CC) $(CFLAGS) -c read_i_files.c

proto_alloc.o: proto_alloc.c
	${CC} ${CFLAGS} -c proto_alloc.c

rwsacs.o: rwsacs.c
	${CC} ${CFLAGS} -c rwsacs.c

rwtextfiles.o:rwtextfiles.c
	$(CC) $(CFLAGS) -c rwtextfiles.c

complex.o: complex.c
	$(CC) $(CFLAGS) -c complex.c

${BIN}/saclst: saclst.c
	${CC} ${CFLAGS} -o ${BIN}/saclst saclst.c

distaz.o:distaz.c
	$(CC) $(CFLAGS) -c distaz.c 

travel_times.o:travel_times.c
	$(CC) $(CFLAGS) -c travel_times.c

# Cleaning
clean:
	\rm -f *.o
	cd FAST_SYNTH; make clean

cleanexe:
	\rm -f $(EXE) $(BIN)/trim_sac_files_qrt fast_synth fast_synth_only_Z
