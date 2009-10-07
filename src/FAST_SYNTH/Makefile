BIN    = ../../bin

##############################################

OBJS =  read_cmt.o   	  \
        read_stats.o 	  \
	../distaz.o   \
        ../proto_alloc.o  \
	../rwsacs.o	  \
	../rwtextfiles.o  \
	get_prefix.o 	  \
        rotate_cmt.o 	  \
	rotate_traces.o	  \
	yyyymmdd2jjj.o 	  \
        save_sac.o        \
        get_length.o	  \
	get_depths.o      \
	../travel_times.o \
	../read_i_files.o 

all:    $(BIN)/fast_synth $(BIN)/fast_synth_only_Z $(BIN)/fast_synth_rot


$(BIN)/fast_synth: fast_synth.o sum_up.o $(OBJS)
	@ $(CC) ${LDFLAGS} -o $(BIN)/fast_synth  fast_synth.o sum_up.o \
	$(OBJS) -lm 

$(BIN)/fast_synth_rot: fast_synth_rot.o sum_up.o $(OBJS)
	@ $(CC) ${LDFLAGS} -o $(BIN)/fast_synth_rot  fast_synth_rot.o sum_up.o \
	$(OBJS) -lm 

$(BIN)/fast_synth_only_Z: fast_synth_only_Z.o sum_up_only_Z.o $(OBJS)
	@ $(CC) ${LDFLAGS} -o $(BIN)/fast_synth_only_Z fast_synth_only_Z.o \
	sum_up_only_Z.o $(OBJS) -lm 

%.o: %.c
	@ ${CC} ${CFLAGS} -c $<  -o $@

clean: 
	\rm -f *.o 
