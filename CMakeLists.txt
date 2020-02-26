cmake_minimum_required(VERSION 3.6)
set(CRYSTFEL_SHORT_VERSION 0.9.1)
set(CRYSTFEL_VERSION ${CRYSTFEL_SHORT_VERSION})
set(CRYSTFEL_API_VERSION 13)
project(crystfel VERSION ${CRYSTFEL_VERSION} LANGUAGES C)

cmake_policy(SET CMP0074 NEW)

include(GNUInstallDirs)

# Set CMAKE_MODULE_PATH to assist cmake in finding our own package definitions
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake")

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)
find_package(TIFF)
find_package(OpenCL)
find_package(PkgConfig)

# Request HDF5 1.10-style API (can't use 1.12-style yet)
add_definitions(-DH5_USE_110_API)

pkg_check_modules (GLIB glib-2.0)
include_directories (${GLIB_INCLUDE_DIRS})
link_directories (${GLIB_LIBRARY_DIRS})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

pkg_search_module(MSGPACK msgpack)
if (MSGPACK_FOUND)
  message(STATUS "Found Messagepack")
else ()
  message(STATUS "MessagePack not found.")
endif ()

pkg_search_module(ZMQ libzmq)
if (ZMQ_FOUND)
  message(STATUS "Found ZMQ")
else ()
  message(STATUS "ZMQ not found.")
endif ()

pkg_search_module(GTK gtk+-3.0)
if (NOT GTK_FOUND)

  # No GTK3.  Try GTK2?
  pkg_search_module(GTK gtk+-2.0)
  if (GTK_FOUND)
    message(STATUS "Using GTK 2")
    pkg_search_module(GDK gdk-2.0)
  endif ()

else ()

  message(STATUS "Using GTK 3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations")
  pkg_search_module(GDK gdk-3.0)

endif ()

# If no version of GTK was found, try for GDK
if (NOT GTK_FOUND)
  message(STATUS "GTK not found. cell_explorer will not be compiled")
  message(STATUS "Looking separately for GDK")
  pkg_search_module(GDK gdk-3.0)
  if (NOT GDK_FOUND)
    pkg_search_module(GDK gdk-2.0)
    if (GDK_FOUND)
      message(STATUS "Using GDK 2")
    endif()
  else ()
    message(STATUS "Using GDK 3")
  endif()
endif()
link_directories(${GTK_LIBRARY_DIRS})

pkg_search_module(CAIRO cairo)
pkg_search_module(GDKPIXBUF gdk-pixbuf-2.0)

include(CheckCCompilerFlag)
check_c_compiler_flag("-fdiagnostics-color=always" HAVE_DIAG_COLOR)
if (HAVE_DIAG_COLOR)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
endif (HAVE_DIAG_COLOR)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H)
set(CMAKE_C_STANDARD 99)

# Common include directories and libraries

set(COMMON_INCLUDES ${HDF5_INCLUDE_DIRS})
set(COMMON_LIBRARIES libcrystfel m GSL::gsl Threads::Threads)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  # If this is a Git checkout, add the revision ID
  message(STATUS "This is a Git checkout: adding revision ID to version number")
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)
  string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1)
  set(CRYSTFEL_VERSION "${CRYSTFEL_VERSION}+${GIT_SHA1}")
else ()
  if (DEFINED CRYSTFEL_SNAPSHOT_ID)
    message(STATUS "This is a packaged snapshot: adding snapshot ID to version number")
    set(CRYSTFEL_VERSION ${CRYSTFEL_VERSION}+${CRYSTFEL_SNAPSHOT_ID})
  else ()
    message(STATUS "This is a packaged release version")
  endif ()

endif ()

include(CheckIncludeFile)
include(CheckLibraryExists)

set(HAVE_CAIRO ${CAIRO_FOUND})
set(HAVE_TIFF ${TIFF_FOUND})
set(HAVE_GTK ${GTK_FOUND})
set(HAVE_OPENCL ${OpenCL_FOUND})
set(HAVE_GDKPIXBUF ${GDKPIXBUF_FOUND})
set(HAVE_GDK ${GDK_FOUND})
set(HAVE_MSGPACK ${MSGPACK_FOUND})
set(HAVE_ZMQ ${ZMQ_FOUND})

set(PACKAGE_VERSION ${PROJECT_VERSION})

check_include_file(fcntl.h HAVE_FCNTL_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(unistd.h HAVE_UNISTD_H)

# Workarounds for weird OpenCL quirks..
if (OpenCL_FOUND)

  # If include path ends with "OpenCL.framework", add the missing "/Headers"
  string(REGEX REPLACE "OpenCL.framework$" "OpenCL.framework/Headers"
         OpenCL_INCLUDE_DIRS ${OpenCL_INCLUDE_DIRS})

  # Check if the header is actually "CL/cl.h"
  check_include_file(CL/cl.h HAVE_CL_CL_H "-I${OpenCL_INCLUDE_DIRS}")

endif (OpenCL_FOUND)

# Check for nice clock function
check_symbol_exists(clock_gettime "time.h" HAVE_CLOCK_GETTIME)

# Check for argp
check_symbol_exists(argp_parse "argp.h" HAVE_ARGP)
if (!HAVE_ARGP)
  message(ERROR "argp not found.  If you're building on Mac OS using libraries"
                "from Homebrew, install argp-standalone.")
endif (!HAVE_ARGP)

# Add the libcrystfel target
add_subdirectory(libcrystfel)

configure_file(config.h.cmake.in config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add tests
enable_testing()
add_subdirectory(tests)

# Install scripts
install(DIRECTORY ${CMAKE_SOURCE_DIR}/scripts
        DESTINATION ${CMAKE_INSTALL_DOCDIR})

# ----------------------------------------------------------------------
# partial_sim

set(PARTIAL_SIM_SOURCES src/partial_sim.c)
add_executable(partial_sim ${PARTIAL_SIM_SOURCES})
target_include_directories(partial_sim PRIVATE ${COMMON_INCLUDES})
target_link_libraries(partial_sim ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES partial_sim)

# ----------------------------------------------------------------------
# pattern_sim

set(PATTERN_SIM_SOURCES src/pattern_sim.c src/diffraction.c)

if(HAVE_OPENCL)
  set(PATTERN_SIM_OPENCL_SOURCES src/cl-utils.c src/diffraction-gpu.c)
endif(HAVE_OPENCL)

add_executable(pattern_sim ${PATTERN_SIM_SOURCES} ${PATTERN_SIM_OPENCL_SOURCES})
target_include_directories(pattern_sim PRIVATE ${COMMON_INCLUDES})
target_link_libraries(pattern_sim ${COMMON_LIBRARIES} ${HDF5_C_LIBRARIES})

if(HAVE_OPENCL)
  target_include_directories(pattern_sim PRIVATE ${OpenCL_INCLUDE_DIRS})
  target_link_libraries(pattern_sim ${OpenCL_LIBRARIES})
endif(HAVE_OPENCL)

list(APPEND CRYSTFEL_EXECUTABLES pattern_sim)

# ----------------------------------------------------------------------
# cell_explorer

if (GTK_FOUND)

  set(CELL_EXPLORER_SOURCES src/cell_explorer.c src/multihistogram.c)

  add_executable(cell_explorer ${CELL_EXPLORER_SOURCES})
  target_include_directories(cell_explorer PRIVATE ${COMMON_INCLUDES} ${GTK_INCLUDE_DIRS})
  target_link_libraries(cell_explorer ${COMMON_LIBRARIES} ${GTK_LIBRARIES})

  list(APPEND CRYSTFEL_EXECUTABLES cell_explorer)

endif (GTK_FOUND)

# ----------------------------------------------------------------------
# render_hkl

if (CAIRO_FOUND)

  set(RENDER_HKL_SOURCES src/render_hkl.c)
  add_executable(render_hkl ${RENDER_HKL_SOURCES})

  target_include_directories(render_hkl PRIVATE ${COMMON_INCLUDES} ${CAIRO_INCLUDE_DIRS})
  target_link_libraries(render_hkl ${COMMON_LIBRARIES} ${CAIRO_LIBRARIES})

  list(APPEND CRYSTFEL_EXECUTABLES render_hkl)

endif (CAIRO_FOUND)

# ----------------------------------------------------------------------
# process_hkl

set(PROCESS_HKL_SOURCES src/process_hkl.c)
add_executable(process_hkl ${PROCESS_HKL_SOURCES})
target_include_directories(process_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(process_hkl ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES process_hkl)

# ----------------------------------------------------------------------
# list_events

set(LIST_EVENTS_SOURCES src/list_events.c)
add_executable(list_events ${LIST_EVENTS_SOURCES})
target_include_directories(list_events PRIVATE ${COMMON_INCLUDES})
target_link_libraries(list_events ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES list_events)

# ----------------------------------------------------------------------
# indexamajig

set(INDEXAMAJIG_SOURCES src/indexamajig.c src/im-sandbox.c src/process_image.c
                        src/time-accounts.c)

if ( ZMQ_FOUND AND MSGPACK_FOUND )
  list(APPEND INDEXAMAJIG_SOURCES src/im-zmq.c)
endif ()

add_executable(indexamajig ${INDEXAMAJIG_SOURCES})
target_include_directories(indexamajig PRIVATE ${COMMON_INCLUDES})
target_link_libraries(indexamajig ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES indexamajig)

if ( ZMQ_FOUND AND MSGPACK_FOUND )
  target_include_directories(indexamajig PRIVATE ${ZMQ_INCLUDE_DIR} ${MSGPACK_INCLUDE_DIR})
  target_link_libraries(indexamajig ${ZMQ_LIBRARIES} ${MSGPACK_LIBRARIES})
endif ()


# ----------------------------------------------------------------------
# get_hkl

set(GET_HKL_SOURCES src/get_hkl.c)
add_executable(get_hkl ${GET_HKL_SOURCES})
target_include_directories(get_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(get_hkl ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES get_hkl)

# ----------------------------------------------------------------------
# compare_hkl

set(COMPARE_HKL_SOURCES src/compare_hkl.c)
add_executable(compare_hkl ${COMPARE_HKL_SOURCES})
target_include_directories(compare_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(compare_hkl ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES compare_hkl)

# ----------------------------------------------------------------------
# check_hkl

set(CHECK_HKL_SOURCES src/check_hkl.c)
add_executable(check_hkl ${CHECK_HKL_SOURCES})
target_include_directories(check_hkl PRIVATE ${COMMON_INCLUDES})
target_link_libraries(check_hkl ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES check_hkl)

# ----------------------------------------------------------------------
# partialator

set(PARTIALATOR_SOURCES src/partialator.c src/post-refinement.c src/merge.c
                        src/rejection.c src/scaling.c)
add_executable(partialator ${PARTIALATOR_SOURCES})
target_include_directories(partialator PRIVATE ${COMMON_INCLUDES})
target_link_libraries(partialator ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES partialator)

# ----------------------------------------------------------------------
# ambigator

set(AMBIGATOR_SOURCES src/ambigator.c)
add_executable(ambigator ${AMBIGATOR_SOURCES})
target_include_directories(ambigator PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ambigator ${COMMON_LIBRARIES} ${HDF5_C_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES ambigator)

# ----------------------------------------------------------------------
# geoptimiser

set(GEOPTIMISER_SOURCES src/geoptimiser.c)
add_executable(geoptimiser ${GEOPTIMISER_SOURCES})
target_include_directories(geoptimiser PRIVATE ${COMMON_INCLUDES})
target_link_libraries(geoptimiser ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES geoptimiser)

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
# whirligig

set(WHIRLIGIG_SOURCES src/whirligig.c)
add_executable(whirligig ${WHIRLIGIG_SOURCES})
target_include_directories(whirligig PRIVATE ${COMMON_INCLUDES})
target_link_libraries(whirligig ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES whirligig)

# ----------------------------------------------------------------------
# make_pixelmap

set(MAKE_PIXELMAP_SOURCES src/make_pixelmap.c)
add_executable(make_pixelmap ${MAKE_PIXELMAP_SOURCES})
target_include_directories(make_pixelmap PRIVATE ${COMMON_INCLUDES})
target_link_libraries(make_pixelmap ${COMMON_LIBRARIES} ${HDF5_C_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES make_pixelmap)

# ----------------------------------------------------------------------
# cell_tool

set(CELL_TOOL_SOURCES src/cell_tool.c)
add_executable(cell_tool ${CELL_TOOL_SOURCES})
target_include_directories(cell_tool PRIVATE ${COMMON_INCLUDES})
target_link_libraries(cell_tool ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES cell_tool)

# ----------------------------------------------------------------------
# crystfel (main graphical user interface)

if (GTK_FOUND)

  set(CRYSTFEL_GUI_SOURCES src/crystfel_gui.c src/crystfelimageview.c)

  add_executable(crystfel ${CRYSTFEL_GUI_SOURCES})
  target_include_directories(crystfel PRIVATE ${COMMON_INCLUDES} ${GTK_INCLUDE_DIRS})
  target_link_libraries(crystfel ${COMMON_LIBRARIES} ${GTK_LIBRARIES})

  list(APPEND CRYSTFEL_EXECUTABLES crystfel)

endif (GTK_FOUND)

# ----------------------------------------------------------------------
# Install targets

set_target_properties(${CRYSTFEL_EXECUTABLES}
                      PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set_target_properties(${CRYSTFEL_EXECUTABLES}
                        PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR};${CMAKE_INSTALL_RPATH}")
endif("${isSystemDir}" STREQUAL "-1")

install(TARGETS ${CRYSTFEL_EXECUTABLES}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES
        doc/man/ambigator.1
	doc/man/cell_explorer.1
	doc/man/check_hkl.1
	doc/man/compare_hkl.1
	doc/man/geoptimiser.1
	doc/man/get_hkl.1
	doc/man/indexamajig.1
	doc/man/list_events.1
	doc/man/partialator.1
	doc/man/partial_sim.1
	doc/man/pattern_sim.1
	doc/man/process_hkl.1
	doc/man/render_hkl.1
	doc/man/whirligig.1
	doc/man/make_pixelmap.1
	doc/man/cell_tool.1
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
