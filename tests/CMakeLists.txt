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

add_executable(integration_check integration_check.c)
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

add_executable(prof2d_check prof2d_check.c)
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

