cmake_minimum_required(VERSION 3.5)

# Set project info
project(VanDerPolAdjoint
    VERSION 1.0
    DESCRIPTION "Van der Pol sensitivity analysis using PETSc's TSAdjoint"
    LANGUAGES C)

# Ensure PETSc determines the compilers
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

execute_process(COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
set(CMAKE_C_COMPILER ${C_COMPILER})

# Enable RPATH linking
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Output compiled binary to top-level build/debug folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# Include PETSc and find it with pkg-config
find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

# Add source files
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.c")

# Define executable
add_executable(vanderpol ${SOURCES})

# Link PETSc
target_link_libraries(vanderpol PRIVATE PkgConfig::PETSC)

# Optionally: include PETSc headers explicitly (often not needed with pkg-config)
include_directories($ENV{PETSC_DIR}/include)
include_directories($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)

# Enable testing support (optional)
enable_testing()
