cmake_minimum_required (VERSION 3.5.1)
project( umap-apps LANGUAGES CXX C )
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_C_STANDARD 99)

# Set a default build type if none was specified
set( default_build_type "Release" )
if ( EXISTS "${CMAKE_SOURCE_DIR}/.git" )
  set( default_build_type "Debug" )
endif()

if ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES )
  message(STATUS 
    "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
      "MinSizeRel" "RelWithDebInfo")
endif()

# Make sure we can find the umap installation
if ( NOT UMAP_INSTALL_DIR )
  message( FATAL_ERROR 
    "UMAP_INCLUDE_DIR not specified, please make sure to include -DUMAP_INSTALL_DIR=\"<umap install dir>\"")
endif()

if ( NOT EXISTS ${UMAP_INSTALL_DIR}/include )
  message( FATAL_ERROR
    "${UMAP_INSTALL_DIR}/include directory does not exist")
endif()

if ( NOT EXISTS ${UMAP_INSTALL_DIR}/lib )
  message( FATAL_ERROR
    "${UMAP_INSTALL_DIR}/lib directory does not exist")
endif()

set (UMAPINCLUDEDIRS "${UMAP_INSTALL_DIR}/include")
set (UMAPLIBDIR "${UMAP_INSTALL_DIR}/lib")

add_subdirectory(src)

