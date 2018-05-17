# - Find NBP - numericalBraggPrediction
# Find the native NBP includes and library
#
#  NBP_INCLUDES    - where to find IndexerBase.h
#  NBP_LIBRARIES   - List of libraries when using NBP.
#  NBP_FOUND       - True if NBP found.

if (NBP_INCLUDES)
  # Already in cache, be silent
  set (NBP_FIND_QUIETLY TRUE)
endif (NBP_INCLUDES)

find_path (NBP_INCLUDES numericalBraggPrediction/ProjectionCalculation.h
           PATHS
		   ${CMAKE_INSTALL_PREFIX}/include)

find_library (NBP_LIBRARIES numericalBraggPrediction
              PATHS
		      ${CMAKE_INSTALL_PREFIX}/lib)

# handle the QUIETLY and REQUIRED arguments and set NBP_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NBP DEFAULT_MSG NBP_LIBRARIES NBP_INCLUDES)
