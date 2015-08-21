add_executable(DoubleTriangle 
  DoubleTriangle.cpp ${SAMPLE_HPP} )
target_link_libraries(DoubleTriangle  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
