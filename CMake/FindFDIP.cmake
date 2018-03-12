# - Find FDIP
# Find the native FDIP includes and library
#
#  FDIP_INCLUDES    - where to find streakFinder.h
#  FDIP_LIBRARIES   - List of libraries when using FDIP.
#  FDIP_FOUND       - True if FDIP found.

if (FDIP_INCLUDES)
  # Already in cache, be silent
  set (FDIP_FIND_QUIETLY TRUE)
endif (FDIP_INCLUDES)

find_path (FDIP_INCLUDES streakFinder.h
           PATHS
		   ${CMAKE_INSTALL_PREFIX}/include
           PATH_SUFFIXES fastDiffractionImageProcessing)

find_library (FDIP_LIBRARIES fastDiffractionImageProcessing
              PATHS
		      ${CMAKE_INSTALL_PREFIX}/lib)

# handle the QUIETLY and REQUIRED arguments and set FDIP_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FDIP DEFAULT_MSG FDIP_LIBRARIES FDIP_INCLUDES)

mark_as_advanced (FDIP_LIBRARIES FDIP_INCLUDES)
