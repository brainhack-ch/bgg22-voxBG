cmake_minimum_required(VERSION 3.22)
project(graph_analysis)
include(cmake-scripts/code-coverage.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Ensures we use all optimizations.
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse4.1 -march=native")

set(LINK_FLAGS "LINKER:--no-undefined,-rpath,/lib/usr/local/include/")

# Required to  python lib
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

set(PYTHON_VERSION 3.10)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#set(Boost_DEBUG ON)
#find_package(Boost 1.79 COMPONENTS ALL)
#find_package(Boost 1.79 COMPONENTS program_options system filesystem REQUIRED PATHS)
find_package(PythonLibs ${PYTHON_VERSION} REQUIRED)
find_package(Python ${PYTHON_VERSION} COMPONENTS Interpreter Development REQUIRED)

add_subdirectory(pybind11)
#find_package(pybind11 2.10 REQUIRED PATHS lib)

#include_directories(${Boost_INCLUDE_DIR})
#message("boost lib: ${Boost_LIBRARIES}, inc:${Boost_INCLUDE_DIR}")

#link_libraries('boost_filesystem')
#link_libraries('boost_system')


set(SOURCES src/Python/Octree.h src/Python/Octree.cpp src/Python/ConstructionNode.cpp src/Python/ConstructionNode.h)
set(SOURCES_TEST src/Python/Octree.h src/Python/Octree.cpp src/Python/Octree.h
        src/Python/ConstructionNode.cpp src/Python/ConstructionNode.h
        test/Python/Octree_test.cpp test/Python/main_test.cpp test/Python/ConstructionNode_test.cpp test/Python/utils.h test/Python/utils.cpp)


# For main program
#add_executable(graphPruner ${SOURCES})
#target_include_directories(graphPruner PRIVATE eigen)
#target_link_options(graphPruner PRIVATE ${LINK_FLAGS})



# For unit tests
add_executable(unit_test ${SOURCES_TEST})
target_link_libraries(unit_test PRIVATE gtest_main)
target_link_libraries(unit_test PRIVATE gmock_main)
target_link_options(unit_test PRIVATE ${LINK_FLAGS})

target_include_directories(unit_test PRIVATE src)
target_include_directories(unit_test PRIVATE eigen)
target_code_coverage(unit_test target_code_coverage_PUBLIC)
#target_compile_features(unit_test, PUBLIC, cxx_std_17)

include(GoogleTest)
gtest_discover_tests(unit_test)


# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

add_library(graphPrunerLib MODULE ${SOURCES})
target_compile_definitions(graphPrunerLib PUBLIC PYTHON_BIND=1)
target_link_libraries(graphPrunerLib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} pybind11::module
        pybind11::lto pybind11::headers pybind11::python_link_helper pybind11::embed pybind11::opt_size
        pybind11::windows_extras)
set(PYBIND11_PYTHON_VERSION 3.10 CACHE STRING "")
target_include_directories(graphPrunerLib PUBLIC ${PYTHON_INCLUDE_DIRS} eigen)
target_link_options(graphPrunerLib PRIVATE ${LINK_FLAGS})

