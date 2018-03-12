# - Find PINKINDEXER
# Find the native PINKINDEXER includes and library
#
#  PINKINDEXER_INCLUDES    - where to find IndexerBase.h
#  PINKINDEXER_LIBRARIES   - List of libraries when using PINKINDEXER.
#  PINKINDEXER_FOUND       - True if PINKINDEXER found.

if (PINKINDEXER_INCLUDES)
  # Already in cache, be silent
  set (PINKINDEXER_FIND_QUIETLY TRUE)
endif (PINKINDEXER_INCLUDES)

find_path (PINKINDEXER_INCLUDES pinkIndexer/PinkIndexer.h
           PATHS
		   ${CMAKE_INSTALL_PREFIX}/include)

find_library (PINKINDEXER_LIBRARIES pinkIndexer
              PATHS
		      ${CMAKE_INSTALL_PREFIX}/lib)

# handle the QUIETLY and REQUIRED arguments and set PINKINDEXER_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PINKINDEXER DEFAULT_MSG PINKINDEXER_LIBRARIES PINKINDEXER_INCLUDES)

mark_as_advanced (PINKINDEXER_LIBRARIES PINKINDEXER_INCLUDES)
