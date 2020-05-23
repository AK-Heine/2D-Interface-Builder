cmake_minimum_required(VERSION 3.17.2)

# set the project name
project(2D-Interface-Builder)

# add the dynamic library
add_library(cxxbackend MODULE "interfacebuilder/backend/backend.cpp")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Python3 COMPONENTS Interpreter Development)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_subdirectory(interfacebuilder/backend/external/pybind11)

target_link_libraries(cxxbackend PUBLIC ${PYTHON_LIBRARIES} pybind11)

pybind11_add_module(backend "interfacebuilder/backend/backend.cpp")
set_target_properties( backend
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
)

