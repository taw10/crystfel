bin_PROGRAMS = pattern_sim process_hkl get_hkl

AM_CFLAGS = -Wall

pattern_sim_SOURCES = pattern_sim.c diffraction.c utils.c image.c cell.c \
                      hdf5-file.c ewald.c detector.c sfac.c intensities.c \
                      reflections.c
pattern_sim_LDADD = @LIBS@

process_hkl_SOURCES = process_hkl.c sfac.c statistics.c cell.c utils.c \
                      reflections.c
process_hkl_LDADD = @LIBS@

get_hkl_SOURCES = get_hkl.c sfac.c cell.c utils.c reflections.c
get_hkl_LDADD = @LIBS@
