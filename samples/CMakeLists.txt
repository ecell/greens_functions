
add_executable(Readertest
  Readertest.cpp )
target_link_libraries(Readertest
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )

add_executable(sim
    Simulator.cpp )
target_link_libraries(sim
  greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} )

add_subdirectory(infinite_plane)
add_subdirectory(SymTriangle)
add_subdirectory(doubleTriangle)
add_subdirectory(cube)
add_subdirectory(closed_cube)
