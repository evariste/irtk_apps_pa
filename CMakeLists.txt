
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# IF(COMMAND cmake_policy)
#   CMAKE_POLICY(SET CMP0003 NEW)
# ENDIF(COMMAND cmake_policy)
if(COMMAND cmake_policy)
  # policy settings ...
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


PROJECT(AppsPA)

# IRTK

if (NOT EXISTS ${IRTK_BINARY_DIR})
  SET(IRTK_BINARY_DIR  NOT_FOUND CACHE PATH "Where the IRTK binaries are")
  message("Please set IRTK_BINARY_DIR")
endif()

SUBDIRS(src applications)

INCLUDE_DIRECTORIES(${AppsPA_SOURCE_DIR}/include)

IF(APPLE)
  # Need an edited version of the original file generated in the IRTK
  # binary directory.  It needs to have quotes put around framework
  # listings in the libraries so that they do not get mangled by cmake
  # and we end up with things like ' -lCarbon ' on the link command
  # line.
  SET(IRTK_INCLUDE_FILE ${AppsPA_SOURCE_DIR}/IRTKInclude.cmake)

  IF (EXISTS ${IRTK_INCLUDE_FILE})
    INCLUDE (${IRTK_INCLUDE_FILE})
  ELSE (EXISTS ${IRTK_INCLUDE_FILE})
    SET(MSG1 "
Unable to find an IRTK include file in local source directory.

   ${IRTK_INCLUDE_FILE}

Please copy over from IRTK binary directory

   ${IRTK_BINARY_DIR}/lib/IRTKInclude.cmake

And insert quotes in link library lists so that -framework XXXX does not
become -framework -lXXXX after inclusion.

XXXX is typically 'Carbon', 'Cocoa', 'Glut' or 'ApplicationServices'
")
    MESSAGE (${MSG1})

  ENDIF (EXISTS ${IRTK_INCLUDE_FILE})

ELSE(APPLE)
  INCLUDE(${IRTK_BINARY_DIR}/lib/IRTKInclude.cmake)
ENDIF(APPLE)



LINK_LIBRARIES(pa++)

if (EXISTS ${IRTK_BINARY_DIR})
  SET(EXECUTABLE_OUTPUT_PATH ${IRTK_BINARY_DIR}/bin CACHE PATH
    "Where to put the executables")
  SET(LIBRARY_OUTPUT_PATH ${IRTK_BINARY_DIR}/lib CACHE PATH
    "Where to put the libraries")
endif()


# if (APPLE)
#    FIND_LIBRARY(CARBON_LIBRARY Carbon)
# endif()
