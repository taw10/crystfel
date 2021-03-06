project(libcrystfel VERSION ${CRYSTFEL_SHORT_VERSION} LANGUAGES C)

find_package(Curses)
pkg_check_modules (XGANDALF xgandalf)
pkg_check_modules (PINKINDEXER pinkIndexer)
pkg_check_modules (FDIP fdip)
find_package(ZLIB REQUIRED)
find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)
find_package(Doxygen)
pkg_search_module(FFTW fftw3)

configure_file(src/libcrystfel-version.c.cmake.in libcrystfel-version.c)

set(HAVE_CURSES ${CURSES_FOUND})
set(HAVE_FFTW ${FFTW_FOUND})
set(HAVE_XGANDALF ${XGANDALF_FOUND})
set(HAVE_PINKINDEXER ${PINKINDEXER_FOUND})
set(HAVE_FDIP ${FDIP_FOUND})

# Recent enough version of zlib?
set(CMAKE_REQUIRED_LIBRARIES "-lz")
check_symbol_exists(gzbuffer "zlib.h" HAVE_GZBUFFER)
unset(CMAKE_REQUIRED_LIBRARIES)

configure_file(config.h.cmake.in config.h)

bison_target(symopp src/symop.y ${CMAKE_CURRENT_BINARY_DIR}/symop-parse.c COMPILE_FLAGS --report=all)
flex_target(symopl src/symop.l ${CMAKE_CURRENT_BINARY_DIR}/symop-lex.c
	DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/symop-lex.h)
add_flex_bison_dependency(symopl symopp)
include_directories(${PROJECT_SOURCE_DIR}/src)

set(LIBCRYSTFEL_SOURCES
    src/reflist.c
    src/utils.c
    src/cell.c
    src/thread-pool.c
    src/image.c
    src/geometry.c
    src/peakfinder8.c
    src/symmetry.c
    src/stream.c
    src/peaks.c
    src/reflist-utils.c
    src/filters.c
    src/index.c
    src/cell-utils.c
    src/integer_matrix.c
    src/crystal.c
    src/integration.c
    src/predict-refine.c
    src/peakfinder8.c
    src/rational.c
    src/spectrum.c
    src/datatemplate.c
    src/colscale.c
    src/detgeom.c
    src/image-cbf.c
    src/image-hdf5.c
    src/fom.c
    ${BISON_symopp_OUTPUTS}
    ${FLEX_symopl_OUTPUTS}
    src/indexers/dirax.c
    src/indexers/mosflm.c
    src/indexers/xds.c
    src/indexers/felix.c
    src/indexers/taketwo.c
    src/indexers/asdf.c
    src/indexers/xgandalf.c
    src/indexers/pinkindexer.c
)

set(LIBCRYSTFEL_HEADERS
    src/reflist.h
    src/symmetry.h
    src/cell.h
    src/reflist-utils.h
    src/thread-pool.h
    src/utils.h
    src/geometry.h
    src/peakfinder8.h
    src/peaks.h
    src/stream.h
    src/index.h
    src/image.h
    src/filters.h
    src/cell-utils.h
    src/integer_matrix.h
    src/crystal.h
    src/predict-refine.h
    src/integration.h
    src/peakfinder8.h
    src/rational.h
    src/spectrum.h
    src/datatemplate.h
    src/colscale.h
    src/detgeom.h
    src/image-msgpack.h
    src/fom.h
)

if (MSGPACK_FOUND)
  set(LIBCRYSTFEL_SOURCES ${LIBCRYSTFEL_SOURCES} src/image-msgpack.c)
endif (MSGPACK_FOUND)

if (DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/doc/index.md index.md)
  set(DOXYGEN_SHOW_INCLUDE_FILES NO)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
  set(DOXYGEN_PREDEFINED HAVE_FFTW)
  doxygen_add_docs(api-docs ${PROJECT_SOURCE_DIR}/src
                   ${CMAKE_CURRENT_BINARY_DIR}/index.md
                   ${PROJECT_SOURCE_DIR}/doc/coding.md)
endif (DOXYGEN_FOUND)

add_library(${PROJECT_NAME} SHARED
	${LIBCRYSTFEL_SOURCES}
	${CMAKE_CURRENT_BINARY_DIR}/libcrystfel-version.c
	${LIBCRYSTFEL_FFTW_SOURCES}
        ${LIBCRYSTFEL_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${CRYSTFEL_API_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${CRYSTFEL_SHORT_VERSION})

# Suppress extra "lib" in front of "libcrystfel"
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_CONFIG_H)

target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE util ${HDF5_C_LIBRARIES} ${ZLIB_LIBRARIES}
                      Threads::Threads GSL::gsl m)

if (XGANDALF_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${XGANDALF_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${XGANDALF_LDFLAGS})
endif (XGANDALF_FOUND)

if (FDIP_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FDIP_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${FDIP_LDFLAGS})
endif (FDIP_FOUND)

if (PINKINDEXER_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${PINKINDEXER_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${PINKINDEXER_LDFLAGS})
endif (PINKINDEXER_FOUND)

if (FFTW_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${FFTW_LDFLAGS})
endif (FFTW_FOUND)

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
