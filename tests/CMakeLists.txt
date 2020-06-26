add_test(NAME process_hkl_check_1
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/process_hkl_check_1 $<TARGET_FILE:process_hkl>)
add_test(NAME process_hkl_check_2
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/process_hkl_check_2 $<TARGET_FILE:process_hkl>)
add_test(NAME process_hkl_check_3
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/process_hkl_check_3 $<TARGET_FILE:process_hkl>)
add_test(NAME process_hkl_check_4
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/process_hkl_check_4 $<TARGET_FILE:process_hkl>)
add_test(NAME partialator_merge_check_1
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/partialator_merge_check_1 $<TARGET_FILE:partialator>)
add_test(NAME partialator_merge_check_2
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/partialator_merge_check_2 $<TARGET_FILE:partialator>)
add_test(NAME partialator_merge_check_3
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/partialator_merge_check_3 $<TARGET_FILE:partialator>)

add_executable(ambi_check ambi_check.c)
target_include_directories(ambi_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ambi_check ${COMMON_LIBRARIES})
add_test(ambi_check ambi_check)

add_executable(cell_check cell_check.c)
target_include_directories(cell_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(cell_check ${COMMON_LIBRARIES})
add_test(cell_check cell_check)

add_executable(centering_check centering_check.c)
target_include_directories(centering_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(centering_check ${COMMON_LIBRARIES})
add_test(centering_check centering_check)

add_executable(integration_check integration_check.c histogram.c)
target_include_directories(integration_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(integration_check PRIVATE ${COMMON_LIBRARIES})
if (CURSES_FOUND)
  target_include_directories(integration_check PRIVATE ${CURSES_INCLUDE_DIRS})
  target_link_libraries(integration_check PRIVATE ${CURSES_LIBRARIES})
endif (CURSES_FOUND)
add_test(integration_check integration_check)

add_executable(list_check list_check.c)
target_include_directories(list_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(list_check ${COMMON_LIBRARIES})
add_test(list_check list_check)

add_executable(prediction_gradient_check prediction_gradient_check.c)
target_include_directories(prediction_gradient_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(prediction_gradient_check ${COMMON_LIBRARIES})
add_test(prediction_gradient_check prediction_gradient_check)

add_executable(prof2d_check prof2d_check.c histogram.c)
target_include_directories(prof2d_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(prof2d_check PRIVATE ${COMMON_LIBRARIES})
if (CURSES_FOUND)
  target_include_directories(prof2d_check PRIVATE ${CURSES_INCLUDE_DIRS})
  target_link_libraries(prof2d_check PRIVATE ${CURSES_LIBRARIES})
endif (CURSES_FOUND)
add_test(prof2d_check prof2d_check)

add_executable(ring_check ring_check.c)
target_include_directories(ring_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ring_check PRIVATE ${COMMON_LIBRARIES})
if (FDIP_FOUND)
  target_include_directories(ring_check PRIVATE ${FDIP_INCLUDES})
  target_link_libraries(ring_check PRIVATE ${FDIP_LIBRARIES})
endif (FDIP_FOUND)
add_test(ring_check ring_check)

add_executable(symmetry_check symmetry_check.c)
target_include_directories(symmetry_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(symmetry_check ${COMMON_LIBRARIES})
add_test(symmetry_check symmetry_check)

add_executable(transformation_check transformation_check.c)
target_include_directories(transformation_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(transformation_check ${COMMON_LIBRARIES})
add_test(transformation_check transformation_check)

if (HAVE_OPENCL)
  add_executable(gpu_sim_check gpu_sim_check.c ../src/diffraction.c
                 ../src/diffraction-gpu.c ../src/cl-utils.c)
  target_include_directories(gpu_sim_check PRIVATE ${COMMON_INCLUDES} ${OpenCL_INCLUDE_DIRS})
  target_link_libraries(gpu_sim_check ${COMMON_LIBRARIES} ${OpenCL_LIBRARIES})
  add_test(gpu_sim_check gpu_sim_check)
endif (HAVE_OPENCL)

add_executable(rational_check rational_check.c)
target_include_directories(rational_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(rational_check ${COMMON_LIBRARIES})
add_test(rational_check rational_check)

add_executable(spectrum_check spectrum_check.c)
target_include_directories(spectrum_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(spectrum_check ${COMMON_LIBRARIES})
add_test(spectrum_check spectrum_check)

add_executable(cellcompare_check cellcompare_check.c)
target_include_directories(cellcompare_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(cellcompare_check ${COMMON_LIBRARIES})
add_test(cellcompare_check cellcompare_check)

add_executable(polarisation_check polarisation_check.c)
target_include_directories(polarisation_check PRIVATE ${COMMON_INCLUDES})
target_link_libraries(polarisation_check ${COMMON_LIBRARIES})
add_test(polarisation_check polarisation_check)

add_executable(evparse1 evparse1.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse1 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse1 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse1 evparse1)

add_executable(evparse2 evparse2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse2 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse2 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse2 evparse2)

add_executable(evparse3 evparse3.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse3 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse3 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse3 evparse3)

add_executable(evparse4 evparse4.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse4 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse4 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse4 evparse4)

add_executable(evparse5 evparse5.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse5 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse5 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse5 evparse5)

add_executable(evparse6 evparse6.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse6 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse6 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse6 evparse6)

add_executable(evparse7 evparse7.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(evparse7 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(evparse7 ${COMMON_LIBRARIES} -lhdf5)
add_test(evparse7 evparse7)

add_executable(ev_enum1 ev_enum1.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(ev_enum1 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ev_enum1 ${COMMON_LIBRARIES} -lhdf5)
add_test(NAME ev_enum1
         COMMAND ev_enum1 ${CMAKE_CURRENT_SOURCE_DIR}/ev_enum1.h5)

add_executable(ev_enum2 ev_enum2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/../libcrystfel/src/image-hdf5.c)
target_include_directories(ev_enum2 PRIVATE ${COMMON_INCLUDES})
target_link_libraries(ev_enum2 ${COMMON_LIBRARIES} -lhdf5)
add_test(NAME ev_enum2
  COMMAND ev_enum2 ${CMAKE_CURRENT_SOURCE_DIR}/ev_enum2.h5
  ${CMAKE_CURRENT_SOURCE_DIR}/ev_enum2.geom)
