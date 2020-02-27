cmake_minimum_required(VERSION 3.14.5)
project(mlcoalsim
        VERSION 2.0.0
        LANGUAGES C)

# Set compliance with C11 standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Do not request compiler specific extensions
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic")

# MPI is required
find_package(MPI 3 REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

# Include ZF_LOG
add_subdirectory(zf_log)

# Include the Benchmark library
add_subdirectory(benchmark)

file(GLOB SOURCES "source/*.c")

# This is handy so we can enable MPI in our IDE of choice by just defining an environment variable
if (DEFINED ENV{WITH_MPI})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DinMPI")
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
set(ALL_EXECUTABLES mlcoalsimX mlcoalsimX_ZnS mlcoalsimXmpi mlcoalsimXmpi_ZnS)
foreach(EXECUTABLE ${ALL_EXECUTABLES})
    add_executable(${EXECUTABLE} ${SOURCES})

    # With MPI
    if (${EXECUTABLE} MATCHES ".*mpi")
        target_link_libraries(${EXECUTABLE} ${MPI_LIBRARIES} -DinMPI)
    endif ()

    # With ZnS
    if (${EXECUTABLE} MATCHES ".*ZnS")
        set_target_properties(${EXECUTABLE} PROPERTIES COMPILE_FLAGS "-DZNS_ACTIVE")
    endif ()
    target_link_libraries(${EXECUTABLE} zf_log)
    install(TARGETS ${EXECUTABLE} DESTINATION ${CMAKE_INSTALL_PREFIX})
endforeach()

