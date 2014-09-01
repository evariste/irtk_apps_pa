# Compatible with CMake 2.6.x
GET_FILENAME_COMPONENT(CMAKE_CURRENT_LIST_DIR C:/Users/pa12/packages/irtk/CMakeLists.txt PATH)

# Prefer project-own Find<Package>.cmake modules for FIND_PACKAGE().
SET(CMAKE_MODULE_PATH C:/Users/pa12/packages/irtk)

# Option to produce condor executables
OPTION(BUILD_CONDOR_EXE "Build condor executables" OFF)
MARK_AS_ADVANCED(BUILD_CONDOR_EXE)
IF (BUILD_CONDOR_EXE)
   SET (CMAKE_X_LIBS "-lGL -L/usr/X11R6/lib -lm -lXext -lXt -lSM -lX11  -lICE -ldl -lnsl")
   SET (CMAKE_MODULE_LINK_FLAGS "-static")
   SET (CMAKE_SHLIB_LINK_FLAGS "-static")
ENDIF (BUILD_CONDOR_EXE)

# Finding GNU scientific library GSL
FIND_PACKAGE(GSL REQUIRED)
IF (GSL_FOUND)
  INCLUDE_DIRECTORIES(C:/Users/pa12/packages/gsl/build)
  LINK_DIRECTORIES()
  LINK_LIBRARIES(C:/Users/pa12/packages/gsl/build/Release/gsl.lib;C:/Users/pa12/packages/gsl/build/Release/gslcblas.lib)
ELSE(GSL_FOUND)
  INCLUDE_DIRECTORIES(C:/Users/pa12/packages/gsl/build)
  LINK_DIRECTORIES()
  LINK_LIBRARIES(C:/Users/pa12/packages/gsl/build/Release/gsl.lib;C:/Users/pa12/packages/gsl/build/Release/gslcblas.lib)
ENDIF (GSL_FOUND)

# add boost dependencies
find_package( Boost 1.46 )

if ( NOT Boost_FOUND )
message( STATUS "Boost could not be found." )
   set( BOOST_ROOT C:/Users/pa12/packages/boost_1_55_0 CACHE PATH "Please enter path to Boost include folder." FORCE )
endif ()

message( STATUS "Boost_INCLUDE_DIRS : '" C:/Users/pa12/packages/boost_1_55_0 "'" )
if ( Boost_FOUND )
	INCLUDE_DIRECTORIES(C:/Users/pa12/packages/boost_1_55_0)
endif()

# Option to build with OpenCV or not.
OPTION(BUILD_WITH_OPENCV "Build using OPENCV" OFF)

IF (BUILD_WITH_OPENCV)
  FIND_PACKAGE(OpenCV REQUIRED)
  ADD_DEFINITIONS(-DHAS_OPENCV)
  INCLUDE_DIRECTORIES()
  LINK_LIBRARIES()
ENDIF (BUILD_WITH_OPENCV)

# Option to produce multi-threaded executables using TBB
OPTION(BUILD_WITH_TBB "Build multi-threaded executables using TBB" OFF)

IF (BUILD_WITH_TBB)
  FIND_PACKAGE(TBB REQUIRED)
  IF (TBB_FOUND)
    # Attention: DO NOT define TBB_DEPRECATED by default or before including the
    #            other TBB header files, in particular parallel_for. The deprecated
    #            behavior of parallel_for is to not choose the chunk size (grainsize)
    #            automatically!
    #
    # http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
    ADD_DEFINITIONS(-DHAS_TBB)
    INCLUDE_DIRECTORIES()
    LINK_DIRECTORIES()
    LINK_LIBRARIES()
  ENDIF (TBB_FOUND)
ENDIF(BUILD_WITH_TBB)

INCLUDE(C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindZLIB.cmake)

IF (ZLIB_FOUND)
  ADD_DEFINITIONS(-DHAS_ZLIB -DHAVE_ZLIB)
  INCLUDE_DIRECTORIES(C:/Users/pa12/packages/zlib123-dll/include)
  LINK_LIBRARIES(C:/Users/pa12/packages/zlib123-dll/lib/zdll.lib)
ENDIF (ZLIB_FOUND)

INCLUDE(C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindPNG.cmake)

FIND_PACKAGE(FLTK REQUIRED)

IF (FLTK_FOUND)
  INCLUDE_DIRECTORIES(C:/Users/pa12/packages/fltk-1.1.10/build/)
  LINK_LIBRARIES (C:/Users/pa12/packages/fltk-1.1.10/build/bin/Release/fltk_images.lib;C:/Users/pa12/packages/fltk-1.1.10/build/bin/Release/fltk.lib;C:/Users/pa12/packages/fltk-1.1.10/build/bin/Release/fltk_gl.lib;C:/Users/pa12/packages/fltk-1.1.10/build/bin/Release/fltk_forms.lib;wsock32;comctl32)
ENDIF (FLTK_FOUND)

IF (APPLE)
  FIND_PROGRAM(APPLE_REZ Rez /Developer/Tools)
  FIND_FILE(FLTK_REZ_FILE mac.r C:/Users/pa12/packages/fltk-1.1.10/build/)
  MARK_AS_ADVANCED(APPLE_REZ)
  MARK_AS_ADVANCED(FLTK_REZ_FILE)
ENDIF (APPLE)

IF (NOT BUILD_CONDOR_EXE)

   INCLUDE (C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindGLUT.cmake)

   IF (GLUT_INCLUDE_PATH)
      INCLUDE_DIRECTORIES(C:/Users/pa12/packages/glut/glut-3.7.6-bin)
   ENDIF(GLUT_INCLUDE_PATH)

   IF (GLUT_LIBRARY)
      LINK_LIBRARIES (C:/Users/pa12/packages/glut/glut-3.7.6-bin/glut32.lib)
   ENDIF (GLUT_LIBRARY)

   INCLUDE (C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindOpenGL.cmake)

   IF (GLU_LIBRARY)
      LINK_LIBRARIES ()
   ENDIF (GLU_LIBRARY)
 
   IF (OPENGL_INCLUDE_PATH)
      INCLUDE_DIRECTORIES()
   ENDIF(OPENGL_INCLUDE_PATH)

   IF (OPENGL_LIBRARY)
      LINK_LIBRARIES (glu32;opengl32)
   ENDIF (OPENGL_LIBRARY)

ENDIF (NOT BUILD_CONDOR_EXE)

INCLUDE (C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/Dart.cmake)

IF (WIN32)
  SET(CMAKE_CXX_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /TP /Ze /W0 /TP /Ze /W0")
  ADD_DEFINITIONS(-DvtkCommon_EXPORTS)
ELSE (WIN32)
  SET(CMAKE_CXX_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /TP /Ze /W0 -Wall -Wno-deprecated -Wno-write-strings")
ENDIF (WIN32)

# Option to wrap python.
OPTION(WRAP_PYTHON "Generate python wrapper libraries" OFF)

IF (WRAP_PYTHON)
  SET(ENABLE_WRAPPING TRUE)
  INCLUDE(C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindPythonLibs.cmake)
ENDIF (WRAP_PYTHON)

# Find SWIG.
IF (ENABLE_WRAPPING)
  INCLUDE(C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindSWIG.cmake)
  IF (SWIG_FOUND)
    INCLUDE()  
    IF (UNIX)
      SET(CMAKE_CXX_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /TP /Ze /W0 -fPIC")
      SET(CMAKE_C_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /TP /Ze /W0 -fPIC")
    ENDIF (UNIX)
  ENDIF (SWIG_FOUND)
ENDIF (ENABLE_WRAPPING)

OPTION(WRAP_CYTHON "Generate cython wrapper libraries" OFF)
IF (WRAP_CYTHON)
  IF (UNIX)
    SET(CMAKE_CXX_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /TP /Ze /W0 -fPIC")
    SET(CMAKE_C_FLAGS " /DWIN32 /D_WINDOWS /W3 /Zm1000 -fPIC")
  ENDIF (UNIX)
ENDIF (WRAP_CYTHON)

ADD_DEFINITIONS(-DIMPERIAL -DANSI -DHAS_CONTRIB -DNO_BOUNDS -DENABLE_UNIX_COMPRESS)

INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/common++/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/geometry++/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/image++/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/contrib++/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/packages/transformation/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/packages/registration/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/packages/registration2/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/packages/segmentation/include)
INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/external/gco-v3.0)

LINK_DIRECTORIES(C:/Users/pa12/packages/irtk/build-vtk/lib) 

# Option to build with PNG or not.
OPTION(BUILD_WITH_PNG "Build using PNG" OFF)

IF (BUILD_WITH_PNG)
  ADD_DEFINITIONS(-DHAS_PNG)
  INCLUDE_DIRECTORIES(C:/Users/pa12/packages/libpng-1.2.37/libpng;C:/Users/pa12/packages/zlib123-dll/include)
  LINK_LIBRARIES(C:/Users/pa12/packages/libpng-1.2.37/Release/libpng.lib;C:/Users/pa12/packages/zlib123-dll/lib/zdll.lib)
ENDIF (BUILD_WITH_PNG)

# Option to build with VTK or not.
OPTION(BUILD_WITH_VTK "Build using VTK" OFF)

IF (BUILD_WITH_VTK)
   # Add VTK
   INCLUDE(C:/Program Files (x86)/CMake 2.8/share/cmake-2.8/Modules/FindVTK.cmake)

   IF (VTK_FOUND)
      ADD_DEFINITIONS(-DHAS_VTK)
      INCLUDE_DIRECTORIES(C:/Users/pa12/packages/VTK/build;C:/Users/pa12/packages/VTK/build/Common;C:/Users/pa12/packages/VTK/build/Utilities;C:/Users/pa12/packages/VTK/build/VolumeRendering;C:/Users/pa12/packages/VTK/build/Rendering;C:/Users/pa12/packages/VTK/build/Charts;C:/Users/pa12/packages/VTK/build/Chemistry;C:/Users/pa12/packages/VTK/build/Utilities/vtkalglib;C:/Users/pa12/packages/VTK/Infovis;C:/Users/pa12/packages/VTK/Geovis;C:/Users/pa12/packages/VTK/Views;C:/Users/pa12/packages/VTK/VolumeRendering;C:/Users/pa12/packages/VTK/Hybrid;C:/Users/pa12/packages/VTK/Widgets;C:/Users/pa12/packages/VTK/Rendering;C:/Users/pa12/packages/VTK/Charts;C:/Users/pa12/packages/VTK/Chemistry;C:/Users/pa12/packages/VTK/Rendering/Testing/Cxx;C:/Users/pa12/packages/VTK/IO;C:/Users/pa12/packages/VTK/Imaging;C:/Users/pa12/packages/VTK/Graphics;C:/Users/pa12/packages/VTK/GenericFiltering;C:/Users/pa12/packages/VTK/Filtering;C:/Users/pa12/packages/VTK/Common;C:/Users/pa12/packages/VTK/Utilities;C:/Users/pa12/packages/VTK/Common/Testing/Cxx;C:/Users/pa12/packages/VTK/build/Utilities/vtknetcdf/include;C:/Users/pa12/packages/VTK/Utilities/vtknetcdf/include;C:/Users/pa12/packages/VTK/build/Utilities/vtklibproj4;C:/Users/pa12/packages/VTK/Utilities/vtklibproj4;C:/Users/pa12/packages/VTK/build/Utilities/DICOMParser;C:/Users/pa12/packages/VTK/Utilities/DICOMParser;C:/Users/pa12/packages/VTK/build/Utilities/vtkfreetype/include;C:/Users/pa12/packages/VTK/Utilities/vtkfreetype/include;C:/Users/pa12/packages/VTK/build/Utilities/LSDyna;C:/Users/pa12/packages/VTK/Utilities/LSDyna;C:/Users/pa12/packages/VTK/build/Utilities/MaterialLibrary;C:/Users/pa12/packages/VTK/Utilities/MaterialLibrary;C:/Users/pa12/packages/VTK/build/Utilities/vtkmetaio;C:/Users/pa12/packages/VTK/Utilities/vtkmetaio;C:/Users/pa12/packages/VTK/build/Utilities/verdict;C:/Users/pa12/packages/VTK/Utilities/verdict;C:/Users/pa12/packages/VTK/build/Utilities/vtkhdf5;C:/Users/pa12/packages/VTK/Utilities/vtkhdf5;C:/Users/pa12/packages/VTK/build/Utilities/vtkhdf5/src;C:/Users/pa12/packages/VTK/Utilities/vtkhdf5/src;C:/Users/pa12/packages/VTK/build/Utilities/vtkhdf5/hl/src;C:/Users/pa12/packages/VTK/Utilities/vtkhdf5/hl/src;C:/Users/pa12/packages/VTK/Utilities/utf8/source;C:/Users/pa12/packages/VTK/Infovis;C:/Users/pa12/packages/VTK/Utilities/vtkalglib;C:/Users/pa12/packages/VTK/Geovis;C:/Users/pa12/packages/VTK/Views)
      LINK_DIRECTORIES(C:/Users/pa12/packages/VTK/build/bin)

      # Add patented library if available
      IF (VTK_KITS MATCHES "PATENTED")
         ADD_DEFINITIONS(-DHAS_VTK_PATENTED)
	  LINK_LIBRARIES (vtkPatented)
      ENDIF (VTK_KITS MATCHES "PATENTED")

       # Add patented library if available
      IF (VTK_KITS MATCHES "HYBRID")
         ADD_DEFINITIONS(-DHAS_VTK_HYBRID)
	 LINK_LIBRARIES (vtkHybrid)
      ENDIF (VTK_KITS MATCHES "HYBRID")

     LINK_LIBRARIES (vtkRendering vtkImaging
      vtkGraphics vtkFiltering vtkIO vtkCommon)
   ENDIF (VTK_FOUND)
ENDIF (BUILD_WITH_VTK)


LINK_LIBRARIES(segmentation++ registration2++ registration++ transformation++ contrib++
image++ geometry++ common++)

IF (UNIX) 
   LINK_LIBRARIES( pthread)
ENDIF ()

# Options to build with nifti, znz and possibly fslio
OPTION(BUILD_WITH_NIFTI "Build using NIFTI support" ON)
IF (BUILD_WITH_NIFTI)
   ADD_DEFINITIONS(-DHAS_NIFTI)
   INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/nifti/niftilib)
   INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/nifti/znzlib)
   LINK_LIBRARIES(znz)
   LINK_LIBRARIES(niftiio)
ENDIF (BUILD_WITH_NIFTI)

# Option to build with cardiac spatial temporal correction, segmentation and motion tracking toolbox.
OPTION(BUILD_CARDIAC "Build with cardiac tool box" OFF)
IF (BUILD_CARDIAC)
   ADD_DEFINITIONS(-DHAS_CARDIAC)
   INCLUDE_DIRECTORIES(C:/Users/pa12/packages/irtk/packages/cardiac/include)
   LINK_LIBRARIES(cardiac++)
ENDIF (BUILD_CARDIAC)

