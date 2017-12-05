set(TEST_NAMES
    GreensFunction3D_test GreensFunction3DAbs_test GreensFunction3DAbsSym_test
    GreensFunction3DRadAbs_test GreensFunction3DRadInf_test GreensFunction3DSym_test
    GreensFunction1DAbsAbs_test GreensFunction1DRadAbs_test
    GreensFunction2DAbs_test GreensFunction2DRefWedgeAbs_test
)

set(test_library_dependencies "")
find_package(Boost COMPONENTS unit_test_framework)
if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
	add_definitions(-DBOOST_TEST_DYN_LINK)
	add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
    link_directories(${Boost_LIBRARY_DIRS})
    set(test_library_dependencies ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} greens_functions ${test_library_dependencies} )
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)
