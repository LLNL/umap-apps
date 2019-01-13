FIND_PACKAGE( OpenMP REQUIRED )
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

  include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

  add_executable(debug_median_calculation debug_median_calculation.cpp)
  install(TARGETS debug_median_calculation
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib/static
          RUNTIME DESTINATION bin )

else()
  message("Skipping median_calculation/debug_median_calculation, OpenMP required")
endif()