add_executable(ClosedCube 
  closed_cube.cpp ${SAMPLE_HPP} )
target_link_libraries(ClosedCube  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
