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
if ( NOT UMAP_INSTALL_PATH )
  message( FATAL_ERROR 
    "UMAP_INSTALL_PATH not specified, please make sure to include -DUMAP_INSTALL_PATH=\"<umap install dir>\"")
endif()

if ( NOT EXISTS ${UMAP_INSTALL_PATH}/include )
  message( FATAL_ERROR
    "${UMAP_INSTALL_PATH}/include directory does not exist")
endif()

if ( NOT EXISTS ${UMAP_INSTALL_PATH}/lib )
  message( FATAL_ERROR
    "${UMAP_INSTALL_PATH}/lib directory does not exist")
endif()

set (UMAPINCLUDEDIRS "${UMAP_INSTALL_PATH}/include")
set (UMAPLIBDIR "${UMAP_INSTALL_PATH}/lib")

add_subdirectory(src)

