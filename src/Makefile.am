bin_PROGRAMS = template_index simulate_patterns

template_index_SOURCES = main.c relrod.c utils.c image.c cell.c hdf5-file.c
template_index_LDADD = @LIBS@ -lm -lgsl -lgslcblas
AM_CFLAGS = -Wall -g

simulate_patterns_SOURCES = sim-main.c relrod.c utils.c image.c cell.c \
                            hdf5-file.c
simulate_patterns_LDADD = @LIBS@ -lm -lgsl -lgslcblas
