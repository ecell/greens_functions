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

add_definitions("-DHAVE_DECL_INFINITY -pg")

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

#add_executable(Main
#    main.cpp ${CPP_FILES} ${HPP_FILES} "compat.h")

#target_link_libraries(Main
#    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

add_library(greens_functions SHARED
    ${CPP_FILES} ${HPP_FILES} "compat.h")

target_link_libraries(greens_functions
    ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

add_subdirectory(tests)
