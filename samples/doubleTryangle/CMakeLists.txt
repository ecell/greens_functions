add_executable(DoubleTryangle 
  DoubleTryangle.cpp ${SAMPLE_HPP} )
target_link_libraries(DoubleTryangle  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
