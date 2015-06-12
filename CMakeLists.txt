cmake_minimum_required(VERSION 2.8)
enable_testing()

set (CPP_FILES
    GreensFunction1DAbsAbs.cpp GreensFunction1DAbsSinkAbs.cpp GreensFunction1DRadAbs.cpp GreensFunction2DAbsSym.cpp GreensFunction2DRadAbs.cpp
    GreensFunction3D.cpp GreensFunction3DAbs.cpp GreensFunction3DAbsSym.cpp GreensFunction3DRadAbs.cpp GreensFunction3DRadAbsBase.cpp
    GreensFunction3DRadInf.cpp GreensFunction3DSym.cpp SphericalBesselGenerator.cpp findRoot.cpp funcSum.cpp freeFunctions.cpp
)

set (HPP_FILES
    factorial.hpp Defs.hpp GreensFunction.hpp GreensFunction1DAbsAbs.hpp GreensFunction1DAbsSinkAbs.hpp GreensFunction1DRadAbs.hpp GreensFunction2DAbsSym.hpp
    GreensFunction2DRadAbs.hpp GreensFunction3D.hpp GreensFunction3DAbs.hpp GreensFunction3DAbsSym.hpp GreensFunction3DRadAbs.hpp GreensFunction3DRadAbsBase.hpp
    GreensFunction3DRadInf.hpp GreensFunction3DSym.hpp SphericalBesselGenerator.hpp findRoot.hpp funcSum.hpp freeFunctions.hpp PairGreensFunction.hpp
)

add_definitions("-DHAVE_DECL_INFINITY")

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

#add_executable(main
#    main.cpp ${CPP_FILES} ${HPP_FILES} "compat.h")

#target_link_libraries(main
#    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

add_library(greens_functions SHARED
    ${CPP_FILES} ${HPP_FILES} "compat.h")

target_link_libraries(greens_functions
    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

#python wrapper

#set (PYTHON_CPP
#GreensFunction1DAbsAbs.cpp GreensFunction1DAbsSinkAbs.cpp GreensFunction1DRadAbs.cpp GreensFunction3D.cpp GreensFunction3DAbs.cpp GreensFunction3DSym.cpp GreensFunction3DAbsSym.cpp GreensFunction3DRadAbs.cpp GreensFunction3DRadAbsBase.cpp GreensFunction3DRadInf.cpp freeFunctions.cpp funcSum.cpp findRoot.cpp SphericalBesselGenerator.cpp
#)

#set (PYTHON_HPP
#    factorial.hpp Defs.hpp GreensFunction.hpp GreensFunction1DAbsAbs.hpp GreensFunction1DAbsSinkAbs.hpp GreensFunction1DRadAbs.hpp
#    GreensFunction3D.hpp GreensFunction3DAbs.hpp GreensFunction3DAbsSym.hpp GreensFunction3DRadAbs.hpp GreensFunction3DRadAbsBase.hpp
#    GreensFunction3DRadInf.hpp GreensFunction3DSym.hpp SphericalBesselGenerator.hpp findRoot.hpp funcSum.hpp freeFunctions.hpp PairGreensFunction.hpp
#)

#find_library(BOOST_PYTHON boost_python)
#include_directories("/usr/include/python2.7")

#add_library( _greens_functions SHARED
#    greens_functions.cpp ${PYTHON_CPP} ${PYTHON_HPP} "compat.h")
#target_link_libraries( _greens_functions
#    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} ${BOOST_PYTHON})

add_subdirectory(tests)
