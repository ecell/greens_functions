cmake_minimum_required(VERSION 2.8)
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
    factorial.hpp Defs.hpp GreensFunction.hpp GreensFunction1DAbsAbs.hpp
    GreensFunction1DAbsSinkAbs.hpp GreensFunction1DRadAbs.hpp
    GreensFunction2DAbsSym.hpp GreensFunction2DRadAbs.hpp GreensFunction2DAbs.hpp
    GreensFunction2DRefWedgeAbs.cpp
    GreensFunction3D.hpp GreensFunction3DAbs.hpp GreensFunction3DAbsSym.hpp
    GreensFunction3DRadAbs.hpp GreensFunction3DRadAbsBase.hpp
    GreensFunction3DRadInf.hpp GreensFunction3DSym.hpp SphericalBesselGenerator.hpp
    findRoot.hpp funcSum.hpp freeFunctions.hpp PairGreensFunction.hpp
)

add_definitions("-DHAVE_DECL_INFINITY -Wall")

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

if (NOT NO_BESSEL_TABLE)
    add_subdirectory(tablegen)
endif()

#add_executable(main
#    main.cpp ${CPP_FILES} ${HPP_FILES} "compat.h")

#target_link_libraries(main
#    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

add_library(greens_functions SHARED
    ${CPP_FILES} ${HPP_FILES} "compat.h")

target_link_libraries(greens_functions
    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

if (NOT NO_BESSEL_TABLE)
    add_dependencies(greens_functions BesselTables)

    install(FILES ${HPP_FILES}
        "${CMAKE_CURRENT_BINARY_DIR}/SphericalBesselTable.hpp"
        "${CMAKE_CURRENT_BINARY_DIR}/CylindricalBesselTable.hpp"
        DESTINATION "include/greens_functions")
endif()

include_directories(${PROJECT_BINARY_DIR}/greens_functions)

install(TARGETS greens_functions DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/greens_functions")

# python wrapper

find_library(BOOST_PYTHON boost_python)
include_directories("/usr/include/python2.7")

add_library( _greens_functions SHARED
    python/greens_functions.cpp ${CPP_FILES} ${HPP_FILES} "compat.h")
target_link_libraries( _greens_functions
    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} ${BOOST_PYTHON})

add_subdirectory(tests)
