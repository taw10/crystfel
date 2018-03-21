# - Try to find the GDK library
# Once done this will define
#
#  GDK_FOUND - system has gdk
#  GDK_INCLUDE_DIRS - the gdk include directory
#  GDK_LIBRARIES - Link these to use gdk
#
# Define GDK_MIN_VERSION for which version desired.
#

include(FindPkgConfig)

if (GDK_FIND_REQUIRED)
  set(_pkgconfig_REQUIRED "REQUIRED")
else (GDK_FIND_REQUIRED)
  set(_pkgconfig_REQUIRED "")
endif (GDK_FIND_REQUIRED)

if (GDK_MIN_VERSION)
  pkg_search_module(GDK ${_pkgconfig_REQUIRED} gdk-2.0>=${GDK_MIN_VERSION})
else (GDK_MIN_VERSION)
  pkg_search_module(GDK ${_pkgconfig_REQUIRED} gdk-2.0)
endif (GDK_MIN_VERSION)

if (GDK_FOUND)
  message(STATUS "Found GDK (using pkg-config): ${GDK_LIBRARIES}")
endif (GDK_FOUND)

# Backup option if we don't have pkg-config
if (NOT GDK_FOUND AND NOT PKG_CONFIG_FOUND)

  find_path(GDK_INCLUDE_DIRS gdk/gdk.h)
  find_library(GDK_LIBRARIES gdk)

  if (GDK_LIBRARIES AND GDK_INCLUDE_DIRS)
    set(GDK_FOUND 1)
    if (NOT GDK_FIND_QUIETLY)
      message(STATUS "Found GDK_LIBRARIES}")
    endif(NOT GDK_FIND_QUIETLY)
  else (GDK_LIBRARIES AND GDK_INCLUDE_DIRS)
    if (GDK_FIND_REQUIRED)
      message(SEND_ERROR "Could not find GDK")
    else (GDK_FIND_REQUIRED)
      if (NOT GDK_FIND_QUIETLY)
        message(STATUS "Could not find GDK")
      endif (NOT GDK_FIND_QUIETLY)
    endif (GDK_FIND_REQUIRED)
  endif (GDK_LIBRARIES AND GDK_INCLUDE_DIRS)

endif (NOT GDK_FOUND AND NOT PKG_CONFIG_FOUND)

