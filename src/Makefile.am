bin_PROGRAMS = pattern_sim integr_sim

AM_CFLAGS = -Wall -g @CFLAGS@

pattern_sim_SOURCES = pattern_sim.c diffraction.c utils.c image.c cell.c \
                      hdf5-file.c ewald.c detector.c sfac.c
pattern_sim_LDADD = @LIBS@

integr_sim_SOURCES = integr_sim.c diffraction.c utils.c ewald.c sfac.c \
                     statistics.c cell.c
integr_sim_LDADD = @LIBS@
