cmake_minimum_required(VERSION 2.6)

# XXX remove `cmake_policy set old` later
cmake_policy(SET CMP0048 OLD)
project(greens_functions)
enable_testing()

set (CPP_FILES
    GreensFunction1DAbsAbs.cpp GreensFunction1DAbsSinkAbs.cpp GreensFunction1DRadAbs.cpp
    GreensFunction2DAbsSym.cpp GreensFunction2DRadAbs.cpp GreensFunction2DAbs.cpp
    GreensFunction2DRefWedgeAbs.cpp
    GreensFunction3D.cpp GreensFunction3DAbs.cpp GreensFunction3DAbsSym.cpp
    GreensFunction3DRadAbs.cpp GreensFunction3DRadAbsBase.cpp GreensFunction3DRadInf.cpp
    GreensFunction3DSym.cpp SphericalBesselGenerator.cpp
    findRoot.cpp funcSum.cpp freeFunctions.cpp
)

set (HPP_FILES
    compat.h
    factorial.hpp Defs.hpp GreensFunction.hpp GreensFunction1DAbsAbs.hpp
    GreensFunction1DAbsSinkAbs.hpp GreensFunction1DRadAbs.hpp
    GreensFunction2DAbsSym.hpp GreensFunction2DRadAbs.hpp GreensFunction2DAbs.hpp
    GreensFunction2DRefWedgeAbs.cpp
    GreensFunction3D.hpp GreensFunction3DAbs.hpp GreensFunction3DAbsSym.hpp
    GreensFunction3DRadAbs.hpp GreensFunction3DRadAbsBase.hpp
    GreensFunction3DRadInf.hpp GreensFunction3DSym.hpp SphericalBesselGenerator.hpp
    findRoot.hpp funcSum.hpp freeFunctions.hpp PairGreensFunction.hpp
)

add_definitions(-DHAVE_CONFIG_H)

# check some functions are supported or not...
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <cmath>
    int main() {
        double a, b;
        sincos(0., &a, &b);
        return 0;
    }" ECELL_GREENS_FUNCTIONS_HAVE_SINCOS)

configure_file(${PROJECT_SOURCE_DIR}/config.h.in
               ${CMAKE_BINARY_DIR}/greens_functions/config.h)
install(FILES ${CMAKE_BINARY_DIR}/greens_functions/config.h DESTINATION "include/greens_functions")
message(STATUS "config.h was generated: ${CMAKE_BINARY_DIR}/greens_functions/config.h")

include_directories(${PROJECT_BINARY_DIR})
# include_directories(${PROJECT_BINARY_DIR}/greens_functions)

option(GREENS_FUNCTIONS_LINK_BOOST_LIBRARY "link boost.unit_test_framework library" OFF)
if(GREENS_FUNCTIONS_LINK_BOOST_LIBRARY)
    find_package(Boost REQUIRED COMPONENTS unit_test_framework)
else()
    find_package(Boost REQUIRED)
endif()
include_directories(${Boost_INCLUDE_DIRS})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

add_library(greens_functions STATIC
    ${CPP_FILES} ${HPP_FILES} "${CMAKE_BINARY_DIR}/greens_functions/config.h")
target_link_libraries(greens_functions INTERFACE GSL::gsl INTERFACE GSL::gslcblas)

if (NOT NO_BESSEL_TABLE)
    message(STATUS "prepairing to generate BesselTable.")
    add_subdirectory("${PROJECT_SOURCE_DIR}/tablegen")

    message(STATUS "adding dependency to BesselTable.")
    add_dependencies(greens_functions BesselTables)

    install(
        FILES
        "${CMAKE_BINARY_DIR}/greens_functions/SphericalBesselTable.hpp"
        "${CMAKE_BINARY_DIR}/greens_functions/CylindricalBesselTable.hpp"
        DESTINATION "include/greens_functions")
else()
    message(STATUS "NO_BESSEL_TABLE is defined. stop generating BesselTable.")
endif()

install(TARGETS greens_functions DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/greens_functions")

# # python wrapper
#
# find_library(BOOST_PYTHON boost_python)
# include_directories("/usr/include/python2.7")
#
# add_library( _greens_functions SHARED
#     python/greens_functions.cpp ${CPP_FILES} ${HPP_FILES} "${CMAKE_CURRENT_BINARY_DIR}/config.h")
# target_link_libraries( _greens_functions
#     ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} ${BOOST_PYTHON})

add_subdirectory("${PROJECT_SOURCE_DIR}/tests")
