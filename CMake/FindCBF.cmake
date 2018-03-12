# - Find CBF
# Find the native CBF includes and library
#
#  CBF_INCLUDES    - where to find cbf.h
#  CBF_LIBRARIES   - List of libraries when using CBF.
#  CBF_FOUND       - True if CBF found.

if (CBF_INCLUDES)
  # Already in cache, be silent
  set (CBF_FIND_QUIETLY TRUE)
endif (CBF_INCLUDES)

find_path (CBF_INCLUDES cbf.h PATH_SUFFIXES libcbf cbf)

find_library (CBF_LIBRARIES NAMES cbf)

# handle the QUIETLY and REQUIRED arguments and set CBF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CBF DEFAULT_MSG CBF_LIBRARIES CBF_INCLUDES)

mark_as_advanced (CBF_LIBRARIES CBF_INCLUDES)
