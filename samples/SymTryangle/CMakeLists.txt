add_executable(SymTryangle 
  SymTryangle.cpp ${SAMPLE_HPP} )
target_link_libraries(SymTryangle  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
