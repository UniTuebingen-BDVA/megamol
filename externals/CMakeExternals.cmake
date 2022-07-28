# MegaMol
# Copyright (c) 2020, MegaMol Dev Team
# All rights reserved.
#

# Require git
find_package(Git REQUIRED)

# Clone external script
if (NOT EXISTS "${CMAKE_BINARY_DIR}/script-externals")
  message(STATUS "Downloading external scripts")
  execute_process(COMMAND
    ${GIT_EXECUTABLE} clone -b v2.6 https://github.com/UniStuttgart-VISUS/megamol-cmake-externals.git script-externals --depth 1
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

# Include external script
include("${CMAKE_BINARY_DIR}/script-externals/cmake/External.cmake")

# Commonly needed for path setup
include(GNUInstallDirs)

#
# Centralized function to require externals to add them once by invoking
# require_external(<EXTERNAL_TARGET>).
#
# Think of this function as a big switch, testing for the name and presence
# of the external target to guard against duplicated targets.
#
function(require_external NAME)
  set(FETCHCONTENT_QUIET ON CACHE BOOL "")

  # ###########################################################################
  # ### Header-only libraries #################################################
  # ###########################################################################

  # asmjit
  if (NAME STREQUAL "asmjit")
    if (TARGET asmjit)
      return()
    endif ()

    add_external_headeronly_project(asmjit INTERFACE
      GIT_REPOSITORY https://github.com/asmjit/asmjit.git
      GIT_TAG "8474400e82c3ea65bd828761539e5d9b25f6bd83")

  # Delaunator
  elseif (NAME STREQUAL "Delaunator")
    if (TARGET Delaunator)
      return()
    endif ()

    add_external_headeronly_project(Delaunator
      GIT_REPOSITORY https://github.com/delfrrr/delaunator-cpp.git
      GIT_TAG "v0.4.0"
      INCLUDE_DIR "include")

  # Eigen
  elseif (NAME STREQUAL "Eigen")
    if (TARGET Eigen)
      return()
    endif ()

    add_external_headeronly_project(Eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG "3.3.7")

  # glowl
  elseif (NAME STREQUAL "glowl")
    if (TARGET glowl)
      return()
    endif ()

    add_external_headeronly_project(glowl
      GIT_REPOSITORY https://github.com/invor/glowl.git
      GIT_TAG "dafee75f11c5d759df30ff651d6763e4e674dd0e"
      INCLUDE_DIR "include")
    target_compile_definitions(glowl INTERFACE GLOWL_OPENGL_INCLUDE_GLAD2)

  # libcxxopts
  elseif (NAME STREQUAL "libcxxopts")
    if (TARGET libcxxopts)
      return()
    endif ()

    add_external_headeronly_project(libcxxopts
      GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
      # we are waiting for v3 which brings allowing unrecognized options
      #GIT_TAG "v2.1.1"
      GIT_TAG "dd45a0801c99d62109aaa29f8c410ba8def2fbf2"
      INCLUDE_DIR "include")

  # mmpld_io
  elseif (NAME STREQUAL "mmpld_io")
    if (TARGET mmpld_io)
      return()
    endif ()

    add_external_headeronly_project(mmpld_io
      GIT_REPOSITORY https://github.com/UniStuttgart-VISUS/mmpld_io.git
      GIT_TAG 0002c64e0be4dddc137e4fe37db4b96361bc79bd
      INCLUDE_DIR "include")

  # nanoflann
  elseif (NAME STREQUAL "nanoflann")
    if (TARGET nanoflann)
      return()
    endif ()

    add_external_headeronly_project(nanoflann
      GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
      GIT_TAG "v1.3.0"
      INCLUDE_DIR "include")

  # tinygltf
  elseif (NAME STREQUAL "tinygltf")
    if (TARGET tinygltf)
      return()
    endif ()

    add_external_headeronly_project(tinygltf
      GIT_REPOSITORY https://github.com/syoyo/tinygltf.git
      GIT_TAG "v2.5.0")
    target_compile_definitions(tinygltf INTERFACE TINYGLTF_NO_INCLUDE_JSON)

  # sim_sort
  elseif (NAME STREQUAL "sim_sort")
    if (TARGET sim_sort)
      return()
    endif ()

    add_external_headeronly_project(sim_sort
      GIT_REPOSITORY https://github.com/alexstraub1990/simultaneous-sort.git
      GIT_TAG 220fdf37fec2d9d3e3f7674194544ee70eb93ee7 # master on 2021-07-26, because nothing was specified here.
      INCLUDE_DIR "include")

  # ###########################################################################
  # ### Built libraries #######################################################
  # ###########################################################################

  # adios2
  elseif (NAME STREQUAL "adios2")
    if (TARGET adios2)
      return()
    endif ()

    if (WIN32)
      set(ADIOS2_LIB "lib/adios2.lib")
    else ()
      set(ADIOS2_LIB "${CMAKE_INSTALL_LIBDIR}/libadios2.a")
    endif ()

    add_external_project(adios2 STATIC
      GIT_REPOSITORY https://github.com/ornladios/ADIOS2.git
      GIT_TAG "v2.5.0"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${ADIOS2_LIB}"
      CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DADIOS2_BUILD_EXAMPLES=OFF
        -DADIOS2_BUILD_TESTING=OFF
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DADIOS2_USE_BZip2=OFF
        -DADIOS2_USE_Fortran=OFF
        -DADIOS2_USE_HDF5=OFF
        -DADIOS2_USE_PNG=OFF
        -DADIOS2_USE_Profiling=OFF
        -DADIOS2_USE_Python=OFF
        -DADIOS2_USE_SST=OFF
        -DADIOS2_USE_SZ=OFF
        -DADIOS2_USE_SysVShMem=OFF
        -DADIOS2_USE_ZFP=OFF
        -DADIOS2_USE_ZeroMQ=OFF
        -DMPI_GUESS_LIBRARY_NAME=${MPI_GUESS_LIBRARY_NAME})

    add_external_library(adios2
      LIBRARY ${ADIOS2_LIB})

  # bhtsne
  elseif (NAME STREQUAL "bhtsne")
    if (TARGET bhtsne)
      return()
    endif ()

    if (WIN32)
      set(BHTSNE_LIB "lib/bhtsne.lib")
    else ()
      set(BHTSNE_LIB "lib/libbhtsne.a")
    endif ()

    add_external_project(bhtsne STATIC
      GIT_REPOSITORY https://github.com/lvdmaaten/bhtsne.git
      GIT_TAG "36b169c88250d0afe51828448dfdeeaa508f13bc"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${BHTSNE_LIB}"
      PATCH_COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/externals/bhtsne/CMakeLists.txt"
        "<SOURCE_DIR>/CMakeLists.txt")

    add_external_library(bhtsne
      LIBRARY ${BHTSNE_LIB})

  # blend2d
  elseif (NAME STREQUAL "blend2d")
    if (TARGET blend2d)
      return()
    endif ()

    if (WIN32)
      set(BLEND2D_LIB "lib/blend2d.lib")
    else ()
      set(BLEND2D_LIB "lib/libblend2d.a")
    endif ()

    require_external(asmjit)
    external_get_property(asmjit SOURCE_DIR)

    add_external_project(blend2d STATIC
      GIT_REPOSITORY https://github.com/blend2d/blend2d.git
      GIT_TAG "8aeac6cb34b00898ae725bd76eb3bb2c7cffcf86"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${BLEND2D_IMPORT_LIB}" "<INSTALL_DIR>/${BLEND2D_LIB}"
      CMAKE_ARGS
        -DASMJIT_DIR=${SOURCE_DIR}
        -DBLEND2D_STATIC=ON)

    add_external_library(blend2d
      DEPENDS asmjit
      INCLUDE_DIR "include"
      LIBRARY ${BLEND2D_LIB})

  # chemfiles
  elseif(NAME STREQUAL "chemfiles")
    if (TARGET chemfiles)
      return()
    endif()

    if (WIN32)
      set(CHEMFILES_LIB "lib/chemfiles.lib")
    else ()
      set(CHEMFILES_LIB "lib/libchemfiles.a")
    endif ()

    add_external_project(chemfiles STATIC
      GIT_REPOSITORY https://github.com/chemfiles/chemfiles.git
      GIT_TAG "0.10.2"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${CHEMFILES_LIB}"
    )

    add_external_library(chemfiles
      INCLUDE_DIR "include"
      LIBRARY ${CHEMFILES_LIB})

  # Corsair CUE SDK
  elseif (NAME STREQUAL "CUESDK")
    if (TARGET CUESDK)
      return()
    endif ()

    FetchContent_Declare(
      cuesdk_archive
      URL https://github.com/CorsairOfficial/cue-sdk/releases/download/v3.0.378/CUESDK_3.0.378.zip)
    FetchContent_GetProperties(cuesdk_archive)
    if (NOT cuesdk_archive_POPULATED)
      FetchContent_Populate(cuesdk_archive)
      add_library(CUESDK SHARED IMPORTED GLOBAL)
      set_target_properties(CUESDK PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${cuesdk_archive_SOURCE_DIR}/include"
        IMPORTED_CONFIGURATIONS "Release"
        IMPORTED_LOCATION "${cuesdk_archive_SOURCE_DIR}/redist/x64/CUESDK.x64_2017.dll"
        IMPORTED_IMPLIB "${cuesdk_archive_SOURCE_DIR}/lib/x64/CUESDK.x64_2017.lib")
      install(DIRECTORY "${cuesdk_archive_SOURCE_DIR}/redist/x64/" DESTINATION "bin" FILES_MATCHING PATTERN "*2017.dll")
    endif ()

  # expat
  elseif (NAME STREQUAL "expat")
    if (TARGET expat)
      return()
    endif ()

    if (WIN32)
      set(EXPAT_LIB "lib/expat<SUFFIX>.lib")
    else ()
      set(EXPAT_LIB "lib/libexpat.a")
    endif ()

    # Files in core were originally at 64f3cf982a156a62c1fdb44d864144ee5871159e
    # This seems to be master at 07.06.2017, somewhere between 2.2.0 and 2.2.1
    add_external_project(expat STATIC
      GIT_REPOSITORY https://github.com/libexpat/libexpat
      GIT_TAG "R_2_2_1"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${EXPAT_LIB}"
      SOURCE_SUBDIR "expat"
      DEBUG_SUFFIX "d"
      CMAKE_ARGS
        -DBUILD_doc=OFF
        -DBUILD_examples=OFF
        -DBUILD_shared=OFF
        -DBUILD_tests=OFF
        -DBUILD_tools=OFF)

    add_external_library(expat
      LIBRARY ${EXPAT_LIB}
      DEBUG_SUFFIX "d")

  # glad
  elseif (NAME STREQUAL "glad")
    if (TARGET glad)
      return()
    endif ()

    if (WIN32)
      set(GLAD_LIB "lib/glad.lib")
    else ()
      set(GLAD_LIB "${CMAKE_INSTALL_LIBDIR}/libglad.a")
    endif ()

    add_external_project(glad STATIC
      SOURCE_DIR glad
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${GLAD_LIB}")

    add_external_library(glad
      PROJECT glad
      LIBRARY ${GLAD_LIB})

  # IceT
  elseif (NAME STREQUAL "IceT")
    if (TARGET IceTCore)
      return()
    endif ()

    if (WIN32)
      set(ICET_CORE_LIB "lib/IceTCore.lib")
      set(ICET_GL_LIB "lib/IceTGL.lib")
      set(ICET_MPI_LIB "lib/IceTMPI.lib")
    else ()
      set(ICET_CORE_LIB "lib/libIceTCore.a")
      set(ICET_GL_LIB "lib/libIceTGL.a")
      set(ICET_MPI_LIB "lib/libIceTMPI.a")
    endif ()

    add_external_project(IceT STATIC
      GIT_REPOSITORY https://gitlab.kitware.com/icet/icet.git
      GIT_TAG abf5bf2b92c0531170c8db2621b375065c7da7c4 # master on 2021-07-26, because nothing was specified here.
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${ICET_CORE_LIB}" "<INSTALL_DIR>/${ICET_GL_LIB}" "<INSTALL_DIR>/${ICET_MPI_LIB}"
      CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DICET_BUILD_TESTING=OFF
        -DMPI_GUESS_LIBRARY_NAME=${MPI_GUESS_LIBRARY_NAME})

    add_external_library(IceTCore
      PROJECT IceT
      LIBRARY ${ICET_CORE_LIB})

    add_external_library(IceTGL
      PROJECT IceT
      LIBRARY ${ICET_GL_LIB})

    add_external_library(IceTMPI
      PROJECT IceT
      LIBRARY ${ICET_MPI_LIB})

  # libigl
  elseif (NAME STREQUAL "libigl")
    if (TARGET libigl)
      return()
    endif ()

    if (WIN32)
      set(LIBIGL_LIB "")
    else ()
      set(LIBIGL_LIB "")
    endif ()

    add_external_headeronly_project(libigl
      GIT_REPOSITORY https://github.com/libigl/libigl.git
      GIT_TAG "v2.1.0"
      INCLUDE_DIR "include")

  # lua
  elseif (NAME STREQUAL "lua")
    if (TARGET lua)
      return()
    endif ()

    if (WIN32)
      set(LUA_LIB "lib/lua.lib")
    else ()
      set(LUA_LIB "lib/liblua.a")
    endif ()

    add_external_project(lua STATIC
      GIT_REPOSITORY https://github.com/lua/lua.git
      GIT_TAG v5.3.5
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${LUA_LIB}"
      PATCH_COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/externals/lua/CMakeLists.txt"
        "<SOURCE_DIR>/CMakeLists.txt"
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/externals/lua/lua.hpp"
        "<SOURCE_DIR>/lua.hpp")

    add_external_library(lua
      LIBRARY ${LUA_LIB})

  # obj-io
  elseif (NAME STREQUAL "obj-io")
    if (TARGET obj-io)
      return()
    endif ()

    add_external_headeronly_project(obj-io INTERFACE
      GIT_REPOSITORY https://github.com/thinks/obj-io.git
      GIT_TAG bfe835200fdff49b45a6de4561741203f85ad028 # master on 2021-07-26, because nothing was specified here.
      INCLUDE_DIR "include/thinks")

  # qhull
  elseif (NAME STREQUAL "qhull")
    if (TARGET qhull)
      return()
    endif ()

    if (WIN32)
      set(QHULL_LIB "lib/qhull<SUFFIX>.lib")
    else ()
      set(QUHULL_LIB "lib/libqhull<SUFFIX>.a")
    endif ()

    add_external_project(qhull STATIC
      GIT_REPOSITORY https://github.com/qhull/qhull.git
      GIT_TAG "v7.3.2"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${QHULL_LIB}"
      DEBUG_SUFFIX _d
      PATCH_COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/externals/qhull/CMakeLists.txt"
        "<SOURCE_DIR>/CMakeLists.txt")

    add_external_library(qhull
      INCLUDE_DIR "include"
      LIBRARY ${QHULL_LIB}
      DEBUG_SUFFIX _d)

  # quickhull
  elseif (NAME STREQUAL "quickhull")
    if (TARGET quickhull)
      return()
    endif ()

    if (WIN32)
      set(QUICKHULL_LIB "lib/quickhull.lib")
    else ()
      set(QUICKHULL_LIB "lib/libquickhull.a")
    endif ()

    add_external_project(quickhull STATIC
      GIT_REPOSITORY https://github.com/akuukka/quickhull.git
      GIT_TAG 4f65e0801b8f60c9a97da2dadbe63c2b46397694 # master on 2021-07-26, because nothing was specified here.
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${QUICKHULL_LIB}"
      PATCH_COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/externals/quickhull/CMakeLists.txt"
        "<SOURCE_DIR>/CMakeLists.txt")

    add_external_library(quickhull
      LIBRARY ${QUICKHULL_LIB})

  # tinyply
  elseif (NAME STREQUAL "tinyply")
    if (TARGET tinyply)
      return()
    endif ()

    if (WIN32)
      set(TNY_LIB "${CMAKE_INSTALL_LIBDIR}/tinyply<SUFFIX>.lib")
    else ()
      set(TNY_LIB "${CMAKE_INSTALL_LIBDIR}/libtinyply<SUFFIX>.a")
    endif ()

    add_external_project(tinyply STATIC
      GIT_REPOSITORY https://github.com/ddiakopoulos/tinyply.git
      GIT_TAG "2.1"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${TNY_LIB}"
      DEBUG_SUFFIX d
      CMAKE_ARGS
        -DSHARED_LIB=OFF)

    add_external_library(tinyply
      LIBRARY ${TNY_LIB}
      DEBUG_SUFFIX d)

  # tracking
  elseif (NAME STREQUAL "tracking")
    if (TARGET tracking)
      return()
    endif ()

    if (NOT WIN32)
      message(WARNING "External 'tracking' requested, but not available on non-Windows systems")
    endif ()

    set(TRACKING_LIB "lib/tracking.lib")
    set(TRACKING_NATNET_LIB "lib/NatNetLib.lib")

    add_external_project(tracking STATIC
      GIT_REPOSITORY https://github.com/UniStuttgart-VISUS/mm-tracking.git
      GIT_TAG "v2.0"
      BUILD_BYPRODUCTS
        "<INSTALL_DIR>/${TRACKING_LIB}"
        "<INSTALL_DIR>/${TRACKING_NATNET_LIB}"
      CMAKE_ARGS
        -DCREATE_TRACKING_TEST_PROGRAM=OFF)

    add_external_library(tracking
      LIBRARY ${TRACKING_LIB})

    add_external_library(natnet
      PROJECT tracking
      LIBRARY ${TRACKING_NATNET_LIB})

    external_get_property(tracking SOURCE_DIR)
    set(tracking_files "${SOURCE_DIR}/tracking/conf/tracking.conf" PARENT_SCOPE)

  # vtkm
  elseif (NAME STREQUAL "vtkm")
    if (TARGET vtkm)
      return()
    endif ()

    set(VTKM_VER 1.4)
    set(LIB_VER 1)

    if (WIN32)
      set(VTKM_CONT_LIB "lib/vtkm_cont-${VTKM_VER}.lib")
      set(VTKM_RENDERER_LIB "lib/vtkm_rendering-${VTKM_VER}.lib")
      set(VTKM_WORKLET_LIB "lib/vtkm_worklet-${VTKM_VER}.lib")
    else ()
      set(VTKM_CONT_LIB "${CMAKE_INSTALL_LIBDIR}/libvtkm_cont-${VTKM_VER}.a")
      set(VTKM_RENDERER_LIB "${CMAKE_INSTALL_LIBDIR}/libvtkm_rendering-${VTKM_VER}.a")
      set(VTKM_WORKLET_LIB "${CMAKE_INSTALL_LIBDIR}/libvtkm_worklet-${VTKM_VER}.a")
    endif ()

    add_external_project(vtkm STATIC
      GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk-m.git
      GIT_TAG "v${VTKM_VER}.0"
      BUILD_BYPRODUCTS
        "<INSTALL_DIR>/${VTKM_CONT_LIB}"
        "<INSTALL_DIR>/${VTKM_RENDERER_LIB}"
        "<INSTALL_DIR>/${VTKM_WORKLET_LIB}"
      CMAKE_ARGS
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DVTKm_ENABLE_CUDA:BOOL=${vtkm_ENABLE_CUDA}
        -DVTKm_ENABLE_TESTING:BOOL=OFF
        -DVTKm_ENABLE_DEVELOPER_FLAGS:BOOL=OFF
        -DVTKm_ENABLE_EXAMPLES:BOOL=OFF
        -DVTKm_INSTALL_ONLY_LIBRARIES:BOOL=ON
        -DVTKm_USE_64BIT_IDS:BOOL=OFF
        #-DCMAKE_BUILD_TYPE=Release
      )

    add_external_library(vtkm
      PROJECT vtkm
      LIBRARY ${VTKM_CONT_LIB})

    add_external_library(vtkm_renderer
      PROJECT vtkm
      LIBRARY ${VTKM_RENDERER_LIB})

    add_external_library(vtkm_worklet
      PROJECT vtkm
      LIBRARY ${VTKM_WORKLET_LIB})

  # zfp
  elseif (NAME STREQUAL "zfp")
    if (TARGET zfp)
      return()
    endif ()

    if (WIN32)
      set(ZFP_LIB "lib/zfp.lib")
    else ()
      set(ZFP_LIB "${CMAKE_INSTALL_LIBDIR}/libzfp.a")
    endif ()

    add_external_project(zfp STATIC
      GIT_REPOSITORY https://github.com/LLNL/zfp.git
      GIT_TAG "0.5.2"
      BUILD_BYPRODUCTS "<INSTALL_DIR>/${ZFP_LIB}"
      CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_UTILITIES=OFF
        -DBUILD_TESTING=OFF
        -DZFP_WITH_ALIGNED_ALLOC=ON
        -DZFP_WITH_CACHE_FAST_HASH=ON
        -DCMAKE_BUILD_TYPE=Release)

    add_external_library(zfp
      LIBRARY ${ZFP_LIB})

  # vr interop mwk-mint
  elseif(NAME STREQUAL "mwk-mint")
    if(TARGET mwk-mint)
      return()
    endif()

    if (MSVC_IDE)
      set(MSVC_TOOLSET "-${CMAKE_VS_PLATFORM_TOOLSET}")
    else ()
      set(MSVC_TOOLSET "")
    endif ()

    if(WIN32)
      set(MWKMint_LIB "${CMAKE_INSTALL_LIBDIR}/interop.lib")
      set(MWKMint_Spout_LIB "${CMAKE_INSTALL_LIBDIR}/Spout2.lib")
      set(MWKMint_ZMQ_LIB "${CMAKE_INSTALL_LIBDIR}/libzmq${MSVC_TOOLSET}-mt-sgd-4_3_5.lib")
    else()
      set(MWKMint_LIB "")
    endif()

    add_external_project(mwk-mint STATIC
      GIT_REPOSITORY https://github.com/UniStuttgart-VISUS/MWK-mint/
      GIT_TAG "master"
      BUILD_BYPRODUCTS
        "<INSTALL_DIR>/${MWKMint_LIB}"
        "<INSTALL_DIR>/${MWKMint_Spout_LIB}"
        "<INSTALL_DIR>/${MWKMint_ZMQ_LIB}"
    )

    add_external_library(interop
      PROJECT mwk-mint
      LIBRARY ${MWKMint_LIB}
    )

    add_external_library(Spout2
      PROJECT mwk-mint
      LIBRARY ${MWKMint_Spout_LIB}
    )

  else ()
    message(FATAL_ERROR "Unknown external required \"${NAME}\"")
  endif ()

  mark_as_advanced(FORCE FETCHCONTENT_BASE_DIR)
  mark_as_advanced(FORCE FETCHCONTENT_FULLY_DISCONNECTED)
  mark_as_advanced(FORCE FETCHCONTENT_QUIET)
  mark_as_advanced(FORCE FETCHCONTENT_UPDATES_DISCONNECTED)
endfunction(require_external)
