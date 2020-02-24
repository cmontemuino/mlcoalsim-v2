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
# No MPI / No ZnS
add_executable(mlcoalsimX ${SOURCES})
set_target_properties(mlcoalsimX PROPERTIES COMPILE_FLAGS "-Wall -pedantic")
target_link_libraries(mlcoalsimX zf_log)
install(TARGETS mlcoalsimX DESTINATION ${CMAKE_INSTALL_PREFIX})

# No MPI / With ZnS
add_executable(mlcoalsimX_ZnS ${SOURCES})
set_target_properties(mlcoalsimX_ZnS PROPERTIES COMPILE_FLAGS "-DZNS_ACTIVE")
target_link_libraries(mlcoalsimX_ZnS zf_log)
install(TARGETS mlcoalsimX_ZnS DESTINATION ${CMAKE_INSTALL_PREFIX})

# With MPI / No ZnS
add_executable(mlcoalsimXmpi ${SOURCES})
target_link_libraries(mlcoalsimXmpi zf_log)
target_link_libraries(mlcoalsimXmpi ${MPI_LIBRARIES} -DinMPI)
install(TARGETS mlcoalsimXmpi DESTINATION ${CMAKE_INSTALL_PREFIX})

# With MPI / With ZnS
add_executable(mlcoalsimXmpi_ZnS ${SOURCES})
set_target_properties(mlcoalsimX_ZnS PROPERTIES COMPILE_FLAGS "-DZNS_ACTIVE")
target_link_libraries(mlcoalsimXmpi_ZnS zf_log)
target_link_libraries(mlcoalsimXmpi_ZnS ${MPI_LIBRARIES} -DinMPI)
install(TARGETS mlcoalsimXmpi_ZnS DESTINATION ${CMAKE_INSTALL_PREFIX})
