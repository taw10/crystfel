cmake_minimum_required(VERSION 3.6)
project(crystfel VERSION 0.6.3 LANGUAGES C)

include(GNUInstallDirs)

# Set CMAKE_MODULE_PATH to assist cmake in finding our own package definitions
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake")

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)
find_package(GTK2 2.0.0 COMPONENTS gtk)
find_package(GdkPixbuf2)
find_package(GDK)
find_package(FFTW)
find_package(CBF)
find_package(TIFF)
find_package(Cairo)
find_package(XGANDALF)
find_package(PINKINDEXER)
find_package(NBP)
find_package(FDIP)
find_package(OpenCL)
find_package(Curses)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H -DDATADIR=\"${CMAKE_INSTALL_FULL_DATADIR}\")
set(CMAKE_C_STANDARD 99)

# Make programs work without LD_LIBRARY_PATH when installed to a non system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Common include directories and libraries

set(COMMON_INCLUDES ${HDF5_INCLUDE_DIRS})
set(COMMON_LIBRARIES libcrystfel ${HDF5_C_LIBRARIES} Threads::Threads GSL::gsl)

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1)

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
set(HAVE_GDKPIXBUF ${GDKPIXBUF_FOUND})
set(HAVE_GDK ${GDK_FOUND})
set(HAVE_CURSES ${CURSES_FOUND})

set(PACKAGE_VERSION ${PROJECT_VERSION})

check_include_file(fcntl.h HAVE_FCNTL_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(unistd.h HAVE_UNISTD_H)

if (OpenCL_FOUND)
  check_include_file(CL/cl.h HAVE_CL_CL_H "-I${OpenCL_INCLUDE_DIRS}")
endif (OpenCL_FOUND)

# Check for nice clock function
check_symbol_exists(clock_gettime "time.h" HAVE_CLOCK_GETTIME)

# Check for CPU affinity functions (Linux-specific)
set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE")
check_symbol_exists(sched_setaffinity "sched.h" HAVE_CPU_AFFINITY)
unset(CMAKE_REQUIRED_DEFINITIONS)

configure_file(config.h.cmake.in config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the libcrystfel target
add_subdirectory(libcrystfel)

# Add tests
enable_testing()
add_subdirectory(tests)

# crystfel.pc
configure_file(libcrystfel/crystfel.pc.in crystfel.pc)
install(FILES ${CMAKE_BINARY_DIR}/crystfel.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# API docs
find_package(GtkDoc 1.25)
if (GTKDOC_FOUND)
  configure_file(doc/reference/libcrystfel/xml/gtkdocentities.ent.cmake.in
                 docs/xml/gtkdocentities.ent)
  configure_file(doc/reference/libcrystfel/CrystFEL-docs.sgml
                 docs/CrystFEL-docs.sgml COPYONLY)
  add_custom_target(documentation VERBATIM COMMAND
                    ${CMAKE_SOURCE_DIR}/doc/reference/libcrystfel/build-docs)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/html
          DESTINATION ${CMAKE_INSTALL_DOCDIR})
endif (GTKDOC_FOUND)

# ----------------------------------------------------------------------
# Build Target
# partial_sim

set(PARTIAL_SIM_SOURCES src/partial_sim.c)

add_executable(partial_sim ${PARTIAL_SIM_SOURCES})
target_include_directories(partial_sim PRIVATE ${COMMON_INCLUDES})
target_link_libraries(partial_sim ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# pattern_sim

set(PATTERN_SIM_SOURCES src/pattern_sim.c src/diffraction.c)

if(HAVE_OPENCL)
  set(PATTERN_SIM_OPENCL_SOURCES src/cl-utils.c src/diffraction-gpu.c)
  set(PATTERN_SIM_RESOURCE_FILES data/diffraction.cl)
endif(HAVE_OPENCL)

add_executable(pattern_sim ${PATTERN_SIM_SOURCES} ${PATTERN_SIM_OPENCL_SOURCES})
target_include_directories(pattern_sim PRIVATE ${COMMON_INCLUDES})
target_link_libraries(pattern_sim ${COMMON_LIBRARIES})

if(HAVE_OPENCL)
  target_include_directories(pattern_sim PRIVATE ${OpenCL_INCLUDE_DIRS})
  target_link_libraries(pattern_sim ${OpenCL_LIBRARIES})
endif(HAVE_OPENCL)

set_target_properties(pattern_sim PROPERTIES
                      RESOURCE "${PATTERN_SIM_RESOURCE_FILES}")

# ----------------------------------------------------------------------
# Build Target
# hdfsee

if (GTK2_FOUND)

  set(HDFSEE_SOURCES src/hdfsee.c src/dw-hdfsee.c src/hdfsee-render.c)
  set(HDFSEE_RESOURCE_FILES data/hdfsee.ui)
  add_executable(hdfsee ${HDFSEE_SOURCES})

  target_include_directories(hdfsee PRIVATE ${COMMON_INCLUDES} ${GTK2_INCLUDE_DIRS})
  target_link_libraries(hdfsee ${COMMON_LIBRARIES} ${GTK2_LIBRARIES})

  if (TIFF_FOUND)
    target_include_directories(hdfsee PRIVATE ${TIFF_INCLUDE_DIRS})
    target_link_libraries(hdfsee ${TIFF_LIBRARIES})
  endif (TIFF_FOUND)

  set_target_properties(hdfsee PROPERTIES RESOURCE "${HDFSEE_RESOURCE_FILES}")

  install (TARGETS hdfsee
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib
           RESOURCE DESTINATION share/crystfel)

endif (GTK2_FOUND)

# ----------------------------------------------------------------------
# Build Target
# cell_explorer

if (GTK2_FOUND)

  set(CELL_EXPLORER_SOURCES src/cell_explorer.c src/multihistogram.c)

  add_executable(cell_explorer ${CELL_EXPLORER_SOURCES})
  target_include_directories(cell_explorer PRIVATE ${COMMON_INCLUDES} ${GTK2_INCLUDE_DIRS})
  target_link_libraries(cell_explorer ${COMMON_LIBRARIES} ${GTK2_LIBRARIES})

  install (TARGETS cell_explorer
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)

endif (GTK2_FOUND)

# ----------------------------------------------------------------------
# Build Target
# render_hkl

if (CAIRO_FOUND)

  set(RENDER_HKL_SOURCES src/render_hkl.c)
  add_executable(render_hkl ${RENDER_HKL_SOURCES})

  target_include_directories(render_hkl PRIVATE ${COMMON_INCLUDES} ${CAIRO_INCLUDE_DIRS})
  target_link_libraries(render_hkl ${COMMON_LIBRARIES} ${CAIRO_LIBRARIES})

  install (TARGETS render_hkl
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)

endif (CAIRO_FOUND)

# ----------------------------------------------------------------------
# Build Target
# process_hkl

set(PROCESS_HKL_SOURCES src/process_hkl.c)
add_executable(process_hkl ${PROCESS_HKL_SOURCES})
target_include_directories(process_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(process_hkl ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# list_events

set(LIST_EVENTS_SOURCES src/list_events.c)
add_executable(list_events ${LIST_EVENTS_SOURCES})
target_include_directories(list_events PRIVATE ${COMMON_INCLUDES})
target_link_libraries(list_events ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# indexamajig

set(INDEXAMAJIG_SOURCES src/indexamajig.c src/im-sandbox.c src/process_image.c
                        src/time-accounts.c)

add_executable(indexamajig ${INDEXAMAJIG_SOURCES})
target_include_directories(indexamajig PRIVATE ${COMMON_INCLUDES})
target_link_libraries(indexamajig ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# get_hkl

set(GET_HKL_SOURCES src/get_hkl.c)
add_executable(get_hkl ${GET_HKL_SOURCES})
target_include_directories(get_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(get_hkl ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# compare_hkl

set(COMPARE_HKL_SOURCES src/compare_hkl.c)
add_executable(compare_hkl ${COMPARE_HKL_SOURCES})
target_include_directories(compare_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(compare_hkl ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# check_hkl

set(CHECK_HKL_SOURCES src/check_hkl.c)
add_executable(check_hkl ${CHECK_HKL_SOURCES})
target_include_directories(check_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(check_hkl ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# partialator

set(PARTIALATOR_SOURCES src/partialator.c src/post-refinement.c src/merge.c
                        src/rejection.c src/scaling.c)
add_executable(partialator ${PARTIALATOR_SOURCES})
target_include_directories(partialator PRIVATE ${COMMON_INCLUDES})
target_link_libraries(partialator ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# ambigator

set(AMBIGATOR_SOURCES src/ambigator.c)
add_executable(ambigator ${AMBIGATOR_SOURCES})
target_include_directories(ambigator PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ambigator ${COMMON_LIBRARIES})

# ----------------------------------------------------------------------
# Build Target
# geoptimiser

set(GEOPTIMISER_SOURCES src/geoptimiser.c src/hdfsee-render.c)
add_executable(geoptimiser ${GEOPTIMISER_SOURCES})
target_include_directories(geoptimiser PRIVATE ${COMMON_INCLUDES})
target_link_libraries(geoptimiser ${COMMON_LIBRARIES})

# Add features one by one so that #ifdef HAVE_XX//#include XX.h always works
# If Cairo, gdk-pixbuf and GDK are all found, geoptimiser will add PNG support
if (CAIRO_FOUND)
  target_include_directories(geoptimiser PRIVATE ${CAIRO_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${CAIRO_LIBRARIES})
endif (CAIRO_FOUND)

if (GDKPIXBUF_FOUND)
  target_include_directories(geoptimiser PRIVATE ${GDKPIXBUF_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${GDKPIXBUF_LIBRARIES})
endif (GDKPIXBUF_FOUND)

if (GDK_FOUND)
  target_include_directories(geoptimiser PRIVATE ${GDK_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${GDK_LIBRARIES})
endif (GDK_FOUND)

if (TIFF_FOUND)
  target_include_directories(geoptimiser PRIVATE ${TIFF_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${TIFF_LIBRARIES})
endif (TIFF_FOUND)

# ----------------------------------------------------------------------
# Build Target
# whirligig

set(WHIRLIGIG_SOURCES src/whirligig.c)
add_executable(whirligig ${WHIRLIGIG_SOURCES})
target_include_directories(whirligig PRIVATE ${COMMON_INCLUDES})
target_link_libraries(whirligig ${COMMON_LIBRARIES})


# ----------------------------------------------------------------------
# install targets

install(TARGETS
	partial_sim
	pattern_sim
	process_hkl
	list_events
	indexamajig
	get_hkl
	compare_hkl
	check_hkl
	partialator
	ambigator
	whirligig
	geoptimiser
	# hdfsee, cell_explorer, render_hkl only if dependencies are satisfied

	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	RESOURCE DESTINATION share/crystfel
)

install(FILES
        doc/man/ambigator.1
	doc/man/cell_explorer.1
	doc/man/check_hkl.1
	doc/man/compare_hkl.1
	doc/man/geoptimiser.1
	doc/man/get_hkl.1
	doc/man/hdfsee.1
	doc/man/indexamajig.1
	doc/man/list_events.1
	doc/man/partialator.1
	doc/man/partial_sim.1
	doc/man/pattern_sim.1
	doc/man/process_hkl.1
	doc/man/render_hkl.1
	doc/man/whirligig.1
	DESTINATION ${CMAKE_INSTALL_MANDIR}/man1
)

install(FILES
        doc/man/crystfel.7
	DESTINATION ${CMAKE_INSTALL_MANDIR}/man7
)

install(FILES
        doc/man/crystfel_geometry.5
	DESTINATION ${CMAKE_INSTALL_MANDIR}/man5
)
