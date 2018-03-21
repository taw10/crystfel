# - Try to find the gdk-pixbuf-2.0 library
# Once done this will define
#
#  GDKPIXBUF_FOUND - system has gdk-pixbuf
#  GDKPIXBUF_INCLUDE_DIRS - the gdk-pixbuf include directory
#  GDKPIXBUF_LIBRARIES - Link these to use gdk-pixbuf
#
# Define GDKPIXBUF_MIN_VERSION for which version desired.
#

include(FindPkgConfig)

if (GdkPixbuf2_FIND_REQUIRED)
  set(_pkgconfig_REQUIRED "REQUIRED")
else (GdkPixbuf2_FIND_REQUIRED)
  set(_pkgconfig_REQUIRED "")
endif (GdkPixbuf2_FIND_REQUIRED)

if (GDKPIXBUF_MIN_VERSION)
  pkg_search_module(GDKPIXBUF ${_pkgconfig_REQUIRED} gdk-pixbuf-2.0>=${GDKPIXBUF_MIN_VERSION})
else (GDKPIXBUF_MIN_VERSION)
  pkg_search_module(GDKPIXBUF ${_pkgconfig_REQUIRED} gdk-pixbuf-2.0)
endif (GDKPIXBUF_MIN_VERSION)

if (GDKPIXBUF_FOUND)
  message(STATUS "Found GdkPixbuf2 (using pkg-config): ${GDKPIXBUF_LIBRARIES}")
endif (GDKPIXBUF_FOUND)

# Backup option if we don't have pkg-config
if (NOT GDKPIXBUF_FOUND AND NOT PKG_CONFIG_FOUND)

  find_path(GDKPIXBUF_INCLUDE_DIRS gdk-pixbuf/gdk-pixbuf.h)
  find_library(GDKPIXBUF_LIBRARIES gdk-pixbuf)

  if (GDKPIXBUF_LIBRARIES AND GDKPIXBUF_INCLUDE_DIRS)
    set(GDKPIXBUF_FOUND 1)
    if (NOT GdkPixbuf2_FIND_QUIETLY)
      message(STATUS "Found GdkPixbuf2: ${GDKPIXBUF_LIBRARIES}")
    endif(NOT GdkPixbuf2_FIND_QUIETLY)
  else (GDKPIXBUF_LIBRARIES AND GDKPIXBUF_INCLUDE_DIRS)
    if (GdkPixbuf2_FIND_REQUIRED)
      message(SEND_ERROR "Could not find GdkPixbuf2")
    else (GdkPixbuf2_FIND_REQUIRED)
      if (NOT GdkPixbuf2_FIND_QUIETLY)
        message(STATUS "Could not find GdkPixbuf2")
      endif (NOT GdkPixbuf2_FIND_QUIETLY)
    endif (GdkPixbuf2_FIND_REQUIRED)
  endif (GDKPIXBUF_LIBRARIES AND GDKPIXBUF_INCLUDE_DIRS)

endif (NOT GDKPIXBUF_FOUND AND NOT PKG_CONFIG_FOUND)

