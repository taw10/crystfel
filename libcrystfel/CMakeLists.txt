project(libcrystfel VERSION ${CRYSTFEL_SHORT_VERSION} LANGUAGES C)

find_package(Curses)
find_package(CBF)
find_package(XGANDALF)
find_package(PINKINDEXER)
find_package(NBP)
find_package(FDIP)
pkg_search_module(FFTW fftw3)

set(HAVE_CURSES ${CURSES_FOUND})
set(HAVE_FFTW ${FFTW_FOUND})
set(HAVE_XGANDALF ${XGANDALF_FOUND})
set(HAVE_FDIP ${FDIP_FOUND})
set(HAVE_CBFLIB ${CBF_FOUND})

# Check for CPU affinity functions (Linux-specific)
set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE")
check_symbol_exists(sched_setaffinity "sched.h" HAVE_CPU_AFFINITY)
unset(CMAKE_REQUIRED_DEFINITIONS)

# Check for nice clock function
check_symbol_exists(clock_gettime "time.h" HAVE_CLOCK_GETTIME)

# Find out where forkpty() is declared
set(CMAKE_REQUIRED_LIBRARIES "-lutil")
check_symbol_exists(forkpty "pty.h" HAVE_FORKPTY_PTY_H)
check_symbol_exists(forkpty "util.h" HAVE_FORKPTY_UTIL_H)
unset(CMAKE_REQUIRED_LIBRARIES)
if(HAVE_FORKPTY_PTY_H)
  message(STATUS "Found forkpty() in pty.h")
elseif(HAVE_FORKPTY_UTIL_H)
  message(STATUS "Found forkpty() in util.h")
else()
  message(SEND_ERROR "Couldn't find forkpty()")
endif()

configure_file(config.h.cmake.in config.h)

set(LIBCRYSTFEL_SOURCES
    src/reflist.c
    src/utils.c
    src/cell.c
    src/detector.c
    src/thread-pool.c
    src/image.c
    src/hdf5-file.c
    src/geometry.c
    src/peakfinder8.c
    src/statistics.c
    src/symmetry.c
    src/stream.c
    src/peaks.c
    src/reflist-utils.c
    src/filters.c
    src/render.c
    src/index.c
    src/dirax.c
    src/mosflm.c
    src/cell-utils.c
    src/integer_matrix.c
    src/crystal.c
    src/xds.c
    src/integration.c
    src/predict-refine.c
    src/histogram.c
    src/events.c
    src/felix.c
    src/peakfinder8.c
    src/taketwo.c
    src/xgandalf.c
)

if (HAVE_FFTW)
	set(LIBCRYSTFEL_FFTW_SOURCES src/asdf.c)
endif (HAVE_FFTW)

set(LIBCRYSTFEL_HEADERS
    src/hdf5-file.h
    src/reflist.h
    src/symmetry.h
    src/cell.h
    src/reflist-utils.h
    src/thread-pool.h
    src/statistics.h
    src/utils.h
    src/detector.h
    src/geometry.h
    src/peakfinder8.h
    src/peaks.h
    src/stream.h
    src/render.h
    src/index.h
    src/image.h
    src/filters.h
    src/dirax.h
    src/mosflm.h
    src/cell-utils.h
    src/integer_matrix.h
    src/crystal.h
    src/xds.h
    src/predict-refine.h
    src/integration.h
    src/histogram.h
    src/events.h
    src/asdf.h
    src/felix.h
    src/peakfinder8.h
    src/taketwo.h
    src/xgandalf.h
)

add_library(${PROJECT_NAME} SHARED
	${LIBCRYSTFEL_SOURCES}
	${LIBCRYSTFEL_FFTW_SOURCES}
        ${LIBCRYSTFEL_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${CRYSTFEL_API_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${CRYSTFEL_SHORT_VERSION})

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_CONFIG_H)

target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE util ${HDF5_C_LIBRARIES} Threads::Threads GSL::gsl m)

if (XGANDALF_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${XGANDALF_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${XGANDALF_LIBRARIES})
endif (XGANDALF_FOUND)

if (FDIP_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FDIP_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${FDIP_LIBRARIES})
endif (FDIP_FOUND)

if (PINKINDEXER_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${PINKINDEXER_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${PINKINDEXER_LIBRARIES})
endif (PINKINDEXER_FOUND)

if (NBP_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${NBP_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${NBP_LIBRARIES})
endif (NBP_FOUND)

if (FFTW_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LDFLAGS})
endif (FFTW_FOUND)

if (CBF_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${CBF_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${CBF_LIBRARIES})
endif (CBF_FOUND)

if (CURSES_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${CURSES_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${CURSES_LIBRARIES})
endif (CURSES_FOUND)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall)
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${LIBCRYSTFEL_HEADERS}")

install(TARGETS libcrystfel
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/crystfel
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# crystfel.pc
configure_file(crystfel.pc.in crystfel.pc)
install(FILES ${CMAKE_BINARY_DIR}/libcrystfel/crystfel.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
