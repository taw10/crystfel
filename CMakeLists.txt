cmake_minimum_required(VERSION 3.6)
set(CRYSTFEL_SHORT_VERSION 0.7.0)
set(CRYSTFEL_VERSION ${CRYSTFEL_SHORT_VERSION})
set(CRYSTFEL_API_VERSION 10)
project(crystfel VERSION ${CRYSTFEL_VERSION} LANGUAGES C)

include(GNUInstallDirs)

# Set CMAKE_MODULE_PATH to assist cmake in finding our own package definitions
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake")

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)
find_package(GTK2 2.0.0 COMPONENTS gtk)
find_package(TIFF)
find_package(OpenCL)
find_package(PkgConfig)
pkg_search_module(CAIRO cairo)
pkg_search_module(GDK gdk-2.0)
pkg_search_module(GDKPIXBUF gdk-pixbuf-2.0)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DHAVE_CONFIG_H -DDATADIR=\"${CMAKE_INSTALL_FULL_DATADIR}\")
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
set(HAVE_GTK ${GTK2_FOUND})
set(HAVE_OPENCL ${OpenCL_FOUND})
set(HAVE_GDKPIXBUF ${GDKPIXBUF_FOUND})
set(HAVE_GDK ${GDK_FOUND})

set(PACKAGE_VERSION ${PROJECT_VERSION})

check_include_file(fcntl.h HAVE_FCNTL_H)
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(unistd.h HAVE_UNISTD_H)

# If we supposedly have OpenCL, check if the header is actually "CL/cl.h"
if (OpenCL_FOUND)
  check_include_file(CL/cl.h HAVE_CL_CL_H "-I${OpenCL_INCLUDE_DIRS}")
endif (OpenCL_FOUND)

# Check for nice clock function
check_symbol_exists(clock_gettime "time.h" HAVE_CLOCK_GETTIME)

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

# API docs
option(BUILD_API_DOCS "Build API documentation" OFF)
if (BUILD_API_DOCS)
  find_package(GtkDoc 1.25)
  if (GTKDOC_FOUND)
    configure_file(doc/reference/libcrystfel/xml/gtkdocentities.ent.cmake.in
                   docs/xml/gtkdocentities.ent)
    configure_file(doc/reference/libcrystfel/CrystFEL-docs.sgml
                   docs/CrystFEL-docs.sgml COPYONLY)
    add_custom_target(documentation ALL VERBATIM COMMAND
                      ${CMAKE_SOURCE_DIR}/doc/reference/libcrystfel/build-docs)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/html
            DESTINATION ${CMAKE_INSTALL_DOCDIR})
  else (GTKDOC_FOUND)
    message(SEND_ERROR "gtk-doc not found")
  endif (GTKDOC_FOUND)
endif (BUILD_API_DOCS)

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
  set(PATTERN_SIM_RESOURCE_FILES data/diffraction.cl)
endif(HAVE_OPENCL)

add_executable(pattern_sim ${PATTERN_SIM_SOURCES} ${PATTERN_SIM_OPENCL_SOURCES})
target_include_directories(pattern_sim PRIVATE ${COMMON_INCLUDES})
target_link_libraries(pattern_sim ${COMMON_LIBRARIES} ${HDF5_C_LIBRARIES})

if(HAVE_OPENCL)
  target_include_directories(pattern_sim PRIVATE ${OpenCL_INCLUDE_DIRS})
  target_link_libraries(pattern_sim ${OpenCL_LIBRARIES})
endif(HAVE_OPENCL)

set_target_properties(pattern_sim PROPERTIES
                      RESOURCE "${PATTERN_SIM_RESOURCE_FILES}")

list(APPEND CRYSTFEL_EXECUTABLES pattern_sim)

# ----------------------------------------------------------------------
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

  list(APPEND CRYSTFEL_EXECUTABLES hdfsee)

endif (GTK2_FOUND)

# ----------------------------------------------------------------------
# cell_explorer

if (GTK2_FOUND)

  set(CELL_EXPLORER_SOURCES src/cell_explorer.c src/multihistogram.c)

  add_executable(cell_explorer ${CELL_EXPLORER_SOURCES})
  target_include_directories(cell_explorer PRIVATE ${COMMON_INCLUDES} ${GTK2_INCLUDE_DIRS})
  target_link_libraries(cell_explorer ${COMMON_LIBRARIES} ${GTK2_LIBRARIES})

  list(APPEND CRYSTFEL_EXECUTABLES cell_explorer)

endif (GTK2_FOUND)

# ----------------------------------------------------------------------
# render_hkl

if (CAIRO_FOUND)

  set(RENDER_HKL_SOURCES src/render_hkl.c)
  add_executable(render_hkl ${RENDER_HKL_SOURCES})

  target_include_directories(render_hkl PRIVATE ${COMMON_INCLUDES} ${CAIRO_INCLUDE_DIRS})
  target_link_libraries(render_hkl ${COMMON_LIBRARIES} ${CAIRO_LDFLAGS})

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

add_executable(indexamajig ${INDEXAMAJIG_SOURCES})
target_include_directories(indexamajig PRIVATE ${COMMON_INCLUDES})
target_link_libraries(indexamajig ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES indexamajig)

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

set(GEOPTIMISER_SOURCES src/geoptimiser.c src/hdfsee-render.c)
add_executable(geoptimiser ${GEOPTIMISER_SOURCES})
target_include_directories(geoptimiser PRIVATE ${COMMON_INCLUDES})
target_link_libraries(geoptimiser ${COMMON_LIBRARIES})
list(APPEND CRYSTFEL_EXECUTABLES geoptimiser)

# Add features one by one so that #ifdef HAVE_XX//#include XX.h always works
# If Cairo, gdk-pixbuf and GDK are all found, geoptimiser will add PNG support
if (CAIRO_FOUND)
  target_include_directories(geoptimiser PRIVATE ${CAIRO_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${CAIRO_LDFLAGS})
endif (CAIRO_FOUND)

if (GDKPIXBUF_FOUND)
  target_include_directories(geoptimiser PRIVATE ${GDKPIXBUF_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${GDKPIXBUF_LDFLAGS})
endif (GDKPIXBUF_FOUND)

if (GDK_FOUND)
  target_include_directories(geoptimiser PRIVATE ${GDK_INCLUDE_DIRS})
  target_link_libraries(geoptimiser ${GDK_LDFLAGS})
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
# Install targets

set_target_properties(${CRYSTFEL_EXECUTABLES}
                      PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set_target_properties(${CRYSTFEL_EXECUTABLES}
                        PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR};${CMAKE_INSTALL_RPATH}")
endif("${isSystemDir}" STREQUAL "-1")

install(TARGETS ${CRYSTFEL_EXECUTABLES}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	RESOURCE DESTINATION share/crystfel)

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
