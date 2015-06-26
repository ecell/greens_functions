add_executable(Infplane 
  Infplane.cpp ${SAMPLE_HPP} )
target_link_libraries(Infplane  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
