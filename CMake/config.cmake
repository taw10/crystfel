include(CheckIncludeFile)
include(CheckLibraryExists)

set(HAVE_CAIRO ${CAIRO_FOUND})
set(HAVE_TIFF ${TIFF_FOUND})
set(HAVE_GTK ${GTK2_FOUND})
set(HAVE_FFTW ${FFTW_FOUND})
set(HAVE_XGANDALF ${XGANDALF_FOUND})
set(HAVE_FDIP ${FDIP_FOUND})
set(HAVE_OPENCL ${OpenCL_FOUND})
set(HAVE_CBFLIB ${CBF_FOUND})


check_include_file(fcntl.h HAVE_FCNTL_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(unistd.h HAVE_UNISTD_H)
if(OpenCL_FOUND)
	check_include_file(CL/cl.h HAVE_CL_CL_H "-I${OpenCL_INCLUDE_DIRS}")
endif(OpenCL_FOUND)

check_library_exists(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)
