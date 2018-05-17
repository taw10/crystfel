# - Find XGANDALF
# Find the native XGANDALF includes and library
#
#  XGANDALF_INCLUDES    - where to find IndexerBase.h
#  XGANDALF_LIBRARIES   - List of libraries when using XGANDALF.
#  XGANDALF_FOUND       - True if XGANDALF found.

if (XGANDALF_INCLUDES)
  # Already in cache, be silent
  set (XGANDALF_FIND_QUIETLY TRUE)
endif (XGANDALF_INCLUDES)

find_path (XGANDALF_INCLUDES xgandalf/IndexerBase.h
           PATHS
		   ${CMAKE_INSTALL_PREFIX}/include)

find_library (XGANDALF_LIBRARIES xgandalf
              PATHS
		      ${CMAKE_INSTALL_PREFIX}/lib)

# handle the QUIETLY and REQUIRED arguments and set XGANDALF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (XGANDALF DEFAULT_MSG XGANDALF_LIBRARIES XGANDALF_INCLUDES)
