project(libcrystfel)

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
)

add_library(${PROJECT_NAME} SHARED
	${LIBCRYSTFEL_SOURCES}
	${LIBCRYSTFEL_FFTW_SOURCES}
        ${LIBCRYSTFEL_HEADERS})

set_target_properties(${PROJECT_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH 1)

set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION ${CRYSTFEL_API_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${CRYSTFEL_VERSION})

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC util ${HDF5_C_LIBRARIES})

if (XGANDALF_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${XGANDALF_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${XGANDALF_LIBRARIES})
endif (XGANDALF_FOUND)

if (FDIP_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FDIP_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${FDIP_LIBRARIES})
endif (FDIP_FOUND)

if (PINKINDEXER_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${PINKINDEXER_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${PINKINDEXER_LIBRARIES})
endif (PINKINDEXER_FOUND)

if (NBP_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${NBP_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${NBP_LIBRARIES})
endif (NBP_FOUND)

if (FFTW_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${FFTW_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${FFTW_LIBRARIES})
endif (FFTW_FOUND)

if (CBF_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${CBF_INCLUDES})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${CBF_LIBRARIES})
endif (CBF_FOUND)

if (CURSES_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${CURSES_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} PUBLIC ${CURSES_LIBRARIES})
endif (CURSES_FOUND)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall)

install (TARGETS libcrystfel
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# crystfel.pc
configure_file(crystfel.pc.in crystfel.pc)
install(FILES ${CMAKE_BINARY_DIR}/libcrystfel/crystfel.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
