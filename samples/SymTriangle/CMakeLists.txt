add_executable(SymTriangle 
  SymTriangle.cpp ${SAMPLE_HPP} )
target_link_libraries(SymTriangle  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
