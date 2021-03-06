cmake_minimum_required(VERSION 3.14.5)
project(mlcoalsim
        VERSION 2.0.0
        LANGUAGES C)

# Set compliance with C11 standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Do not request compiler specific extensions
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Wextra")

# MPI is required
find_package(MPI 3.1 REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

enable_testing()

# Include testing library
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/munit/munit.c")
    find_program(GIT git) 
    if(GIT AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
        execute_process(COMMAND ${GIT} submodule update --init --recursive)
    else()
        message (FATAL_ERROR "It looks like you don't have submodules checked out.
If you are compiling from git, make sure to run `git submodule update --init --recursive', which will invoke that command for you.")
    endif()
endif()

# Include third-party libs
add_subdirectory("${CMAKE_SOURCE_DIR}/lib")

# Include ZF_LOG
add_subdirectory(${CMAKE_SOURCE_DIR}/zf_log)

# Include the Benchmark library
add_subdirectory(${CMAKE_SOURCE_DIR}/benchmark)

# Include string module
add_subdirectory("${CMAKE_SOURCE_DIR}/source/libstr")
target_compile_options(libstr PRIVATE -Wall -Wextra)

# Include the test module
add_subdirectory("${CMAKE_SOURCE_DIR}/test")

file(GLOB SOURCES "source/*.c")

# This is handy so we can enable MPI/Benchmark in our IDE of choice by just defining an environment variable
if (DEFINED ENV{WITH_MPI})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DinMPI")
endif()
if (DEFINED ENV{WITH_BENCHMARK_ENABLED})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DBENCHMARK_ENABLED")
endif()

# Set default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'RelWithDebInfo' as none was specified.")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

####################
### Executables  ###
####################
# No MPI / No ZnS       : mlcoalsimX
# No MPI / With ZnS     : mlcoalsimX_ZnS
# With MPI / No ZnS     : mlcoalsimXmpi
# With MPI / With ZnS   : mlcoalsimXmpi_ZnS
# mlcoalsimB is same as mlcoalsimX but with benchmarking enabled
set(ALL_EXECUTABLES mlcoalsimX mlcoalsimX_ZnS mlcoalsimXmpi mlcoalsimXmpi_ZnS mlcoalsimB mlcoalsimBmpi)
foreach(EXECUTABLE ${ALL_EXECUTABLES})
    add_executable(${EXECUTABLE} ${SOURCES})

    # With MPI
    if (${EXECUTABLE} MATCHES ".*mpi")
        target_link_libraries(${EXECUTABLE} ${MPI_LIBRARIES})
        set_target_properties(${EXECUTABLE} PROPERTIES COMPILE_FLAGS "-DinMPI")
    endif ()

    # With ZnS
    if (${EXECUTABLE} MATCHES ".*ZnS")
        get_property(CURRENT_FLAGS TARGET ${EXECUTABLE} PROPERTY COMPILE_FLAGS)
        set_target_properties(${EXECUTABLE} PROPERTIES COMPILE_FLAGS "${CURRENT_FLAGS} -DZNS_ACTIVE")
    endif ()

    # With Benchmarks
    if (${EXECUTABLE} MATCHES "mlcoalsimB.*")
        get_property(CURRENT_FLAGS TARGET ${EXECUTABLE} PROPERTY COMPILE_FLAGS)
        set_target_properties(${EXECUTABLE} PROPERTIES COMPILE_FLAGS "${CURRENT_FLAGS} -DBENCHMARK_ENABLED")
    endif ()

    target_link_libraries(${EXECUTABLE} zf_log benchmark libstr)
    install(TARGETS ${EXECUTABLE} DESTINATION ${CMAKE_INSTALL_PREFIX})
endforeach()
