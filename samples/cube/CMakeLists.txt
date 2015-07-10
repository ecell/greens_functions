add_executable(Cube 
  cube.cpp ${SAMPLE_HPP} )
target_link_libraries(Cube  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
