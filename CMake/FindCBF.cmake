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

find_path (CBF_CBF_INCLUDES cbf/cbf.h)
find_path (CBFLIB_CBF_INCLUDES cbflib/cbf.h)

if (NOT CBF_CBF_INCLUDES MATCHES "NOTFOUND")
  message(STATUS "Found cbf/cbf.h in ${CBF_CBF_INCLUDES}")
  set(CBF_INCLUDES ${CBF_CBF_INCLUDES})
  set(HAVE_CBF_CBF_H ON)
elseif (NOT CBFLIB_CBF_INCLUDES MATCHES "NOTFOUND")
  message(STATUS "Found cbflib/cbf.h in ${CBFLIB_CBF_INCLUDES}")
  set(CBF_INCLUDES ${CBFLIB_CBF_INCLUDES})
  set(HAVE_CBFLIB_CBF_H ON)
endif (NOT CBF_CBF_INCLUDES MATCHES "NOTFOUND")

find_library (CBF_LIBRARIES NAMES cbf)

# handle the QUIETLY and REQUIRED arguments and set CBF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CBF DEFAULT_MSG CBF_LIBRARIES CBF_INCLUDES)
