add_executable(make_sjy_table make_sjy_table.cpp)
add_executable(make_cjy_table make_cjy_table.cpp)

target_link_libraries(make_sjy_table gslcblas)
target_link_libraries(make_cjy_table gsl gslcblas)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/greens_functions/SphericalBesselTable.hpp
    COMMAND make_sjy_table ${CMAKE_BINARY_DIR}/greens_functions/SphericalBesselTable.hpp
    DEPENDS make_sjy_table)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/greens_functions/CylindricalBesselTable.hpp
    COMMAND make_cjy_table ${CMAKE_BINARY_DIR}/greens_functions/CylindricalBesselTable.hpp
    DEPENDS make_cjy_table)

add_custom_target(
    BesselTables ALL DEPENDS
    ${CMAKE_BINARY_DIR}/greens_functions/SphericalBesselTable.hpp
    ${CMAKE_BINARY_DIR}/greens_functions/CylindricalBesselTable.hpp)
