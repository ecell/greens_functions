add_executable(Cube 
  cube.cpp ${SAMPLE_HPP} )
target_link_libraries(Cube  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )

add_executable(Octa 
  octahedron.cpp ${SAMPLE_HPP} )
target_link_libraries(Octa  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )

add_executable(Conc 
  concave.cpp ${SAMPLE_HPP} )
target_link_libraries(Conc  
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )
