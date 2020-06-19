cmake_minimum_required(VERSION 3.17.2)

# set the project name
project(2D-Interface-Builder)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")

# specify sources
file(GLOB SOURCES    
    interfacebuilder/backend/*.cpp
)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES Release))
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
  # add the dynamic library
  set(TARGET cxxbackend)
  add_library(${TARGET} MODULE ${SOURCES})  
  message("TARGET: ${TARGET}")
endif()

if(CMAKE_BUILD_TYPE MATCHES TESTING)
  message("BUILD TYPE: ${CMAKE_BUILD_TYPE}")
  set(TARGET backend.exe)
  add_executable(${TARGET} ${SOURCES})  
  message("TARGET: ${TARGET}")
endif()


find_package(Python3 COMPONENTS Interpreter Development)
#find_package(PythonLibs REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

add_subdirectory(interfacebuilder/backend/external/pybind11)
target_link_libraries(${TARGET} PUBLIC ${Python3_LIBRARY} pybind11)

if(CMAKE_BUILD_TYPE MATCHES Release)
  pybind11_add_module(backend ${SOURCES})
  set_target_properties(backend
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/interfacebuilder"
  )
endif()

