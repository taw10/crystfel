bin_PROGRAMS = pattern_sim process_hkl

AM_CFLAGS = -Wall

pattern_sim_SOURCES = pattern_sim.c diffraction.c utils.c image.c cell.c \
                      hdf5-file.c ewald.c detector.c sfac.c intensities.c
pattern_sim_LDADD = @LIBS@

process_hkl_SOURCES = process_hkl.c sfac.c statistics.c cell.c utils.c
process_hkl_LDADD = @LIBS@
