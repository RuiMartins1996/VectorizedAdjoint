cmake_minimum_required(VERSION 3.1)
project(Vanderpol VERSION 1.0.0 DESCRIPTION "Implements sensitivity analysis for the VanDerPol equation")
set(CMAKE_CXX_STANDARD 17)

# Print project description
message(STATUS "Project Description: ${PROJECT_DESCRIPTION}")

# Set sources (if any)
set(SOURCES main.cpp)

# Set executable name
set(EXECUTABLE_NAME vanderpol)

# Enable AVX2 instructions and suppress warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wno-ignored-attributes")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx2 -Wno-ignored-attributes")

# Build the executable
add_executable(${EXECUTABLE_NAME} ${SOURCES})

# Include directories for lib (assuming it's in the parent directory)
include_directories(${CMAKE_SOURCE_DIR}/../../lib/include)

# Include directories for AADC (assuming it's in the parent directory)
include_directories(${CMAKE_SOURCE_DIR}/../../aadc/include)

# Set AADC source directory
set(AADC_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../../aadc)

# Determine AADC library directory and type
if (NOT DEFINED AADC_LIB_DIR)
    set(AADC_LIB_DIR "lib")
endif()

if (NOT DEFINED USE_SHARED_AADC)
    option(USE_SHARED_AADC "Use shared AADC lib" ON)
endif()

# Define AADC library import type
if (USE_SHARED_AADC)
    add_library(aadc SHARED IMPORTED)
    set_target_properties(aadc PROPERTIES IMPORTED_LOCATION ${AADC_SOURCE_DIR}/${AADC_LIB_DIR}/libaadc.so)
else()
    add_library(aadc STATIC IMPORTED)
    set_target_properties(aadc PROPERTIES IMPORTED_LOCATION ${AADC_SOURCE_DIR}/${AADC_LIB_DIR}/libaadc_static.a)
endif()

# Link against AADC library
target_link_libraries(${EXECUTABLE_NAME} PUBLIC aadc)

# Link Boost (assuming it's required for the project)
find_package(Boost CONFIG REQUIRED)
target_link_libraries(${EXECUTABLE_NAME} PUBLIC Boost::boost)
