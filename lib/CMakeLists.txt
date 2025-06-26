cmake_minimum_required(VERSION 3.9)

project(backpropagation VERSION 1.0.0 DESCRIPTION "Library to implement discrete adjoint sensitivity analysis for optimisation problems involving ODEs enabled with features such as automatic differentiation and SIMD vectorisation.")

set(CMAKE_CXX_STANDARD 17)

# Enable AVX2 instructions and suppress warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wno-ignored-attributes")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx2 -Wno-ignored-attributes")

# Include directories
set(PUBLIC_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include/)

# Specify path to aadc library
set(AADC_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../aadc)

# Choose between shared or static linking
option(USE_SHARED_AADC "Use shared AADC lib" ON)

if (USE_SHARED_AADC)
    # Linking with shared library
    add_library(aadc SHARED IMPORTED)
    SET_TARGET_PROPERTIES(aadc PROPERTIES
        IMPORTED_LOCATION ${AADC_SOURCE_DIR}/lib/libaadc.so  # Adjust path as necessary
    )
else()
    # Linking with static library
    add_library(aadc STATIC IMPORTED)
    SET_TARGET_PROPERTIES(aadc PROPERTIES
        IMPORTED_LOCATION ${AADC_SOURCE_DIR}/lib/libaadc_static.a  # Adjust path as necessary
    )
endif()

# Build the library as INTERFACE (header-only)
add_library(backpropagation INTERFACE)

# Set include directories
target_include_directories(backpropagation INTERFACE
    ${PUBLIC_INCLUDE_DIRECTORIES}
    ${AADC_SOURCE_DIR}/include  # Include aadc headers
)

# Link Boost to the library
find_package(Boost CONFIG REQUIRED)
target_link_libraries(backpropagation INTERFACE Boost::boost)

# Link aadc library
target_link_libraries(backpropagation INTERFACE aadc)

# Set library properties
set_target_properties(backpropagation
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PUBLIC_INCLUDE_DIRECTORIES}"
)

# Installation
install(TARGETS backpropagation
    EXPORT ${PROJECT_NAME}Targets
    PUBLIC_HEADER DESTINATION include
)

install(FILES ${INCLUDE_FILES} DESTINATION include)
