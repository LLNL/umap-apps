project(median_calculation)

FIND_PACKAGE(CURL QUIET)

find_library( BZ2_LIBRARY libbz2.so )
find_library( CFITS_LIBRARY libcfitsio.a PATHS ${CFITS_LIBRARY_PATH} )
find_path( CFITS_INCLUDE_DIR fitsio.h PATHS ${CFITS_INCLUDE_PATH} )

if ( CFITS_LIBRARY AND CFITS_INCLUDE_DIR )
  set(CFITS_LIBRARIES ${CFITS_LIBRARY} ${CURL_LIBRARIES} ${BZ2_LIBRARY} )
  include_directories( ${CFITS_INCLUDE_DIR} )

  FIND_PACKAGE( OpenMP REQUIRED )
  if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

      include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${UMAPINCLUDEDIRS})

      add_executable(test_median_calculation test_median_calculation.cpp)
      target_link_libraries(test_median_calculation ${UMAPLIBDIR}/libumap.a ${CFITS_LIBRARIES})
      install(TARGETS test_median_calculation
              LIBRARY DESTINATION lib
              ARCHIVE DESTINATION lib/static
              RUNTIME DESTINATION bin )

      add_executable(run_random_vector run_random_vector.cpp)
      target_link_libraries(run_random_vector ${UMAPLIBDIR}/libumap.a ${CFITS_LIBRARIES})
      install(TARGETS run_random_vector
              LIBRARY DESTINATION lib
              ARCHIVE DESTINATION lib/static
              RUNTIME DESTINATION bin )

  else()
    message("Skipping median_calculation, OpenMP required")
  endif()
else()
  message( STATUS 
    "Skipping median_calculation, (missing CFITS_LIBRARY_PATH/CFITS_INCLUDE_PATH)"
  )
endif()

add_subdirectory(debug_program)