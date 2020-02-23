cmake_minimum_required(VERSION 3.14.5)
project(mlcoalsimv2
        VERSION 2.0.0
        LANGUAGES C)

# Set compliance with C11 standard
set(CMAKE_C_STANDARD 11)
# Do not request compiler specific extensions
set(CMAKE_C_EXTENSIONS False)

set(COMPILER_FLAGS, "-O3 -I. -lm -Wall -pedantic")
set(CMAKE_C_FLAGS "${COMPILER_FLAGS}")

# MPI is required
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

# Include ZF_LOG
add_subdirectory(zf_log)

set(SOURCE_FILES
        source/distrib.c
        source/get_mod_sm.c
        source/input_hh3_sm.c
        source/main_sm2.c
        source/ms2_sm.c
        source/neut_tests.c
        source/ran1.c
        source/streec2_sm.c)

# This is handy so we can enable MPI in our IDE of choice by just defining an environment variable
if (DEFINED ENV{WITH_MPI})
    set(CMAKE_C_FLAGS "${COMPILER_FLAGS} -DinMPI")
endif()

####################
### Executables  ###
####################
# No MPI / No ZnS
add_executable(mlcoalsimX ${SOURCE_FILES})
set_target_properties(mlcoalsimX PROPERTIES COMPILE_FLAGS "-Wall -pedantic")
target_link_libraries(mlcoalsimX zf_log)
install(TARGETS mlcoalsimX DESTINATION ${CMAKE_INSTALL_PREFIX})

# No MPI / With ZnS
add_executable(mlcoalsimX_ZnS ${SOURCE_FILES})
set_target_properties(mlcoalsimX_ZnS PROPERTIES COMPILE_FLAGS "-DZNS_ACTIVE")
target_link_libraries(mlcoalsimX_ZnS zf_log)
install(TARGETS mlcoalsimX_ZnS DESTINATION ${CMAKE_INSTALL_PREFIX})

# With MPI / No ZnS
add_executable(mlcoalsimXmpi ${SOURCE_FILES})
target_link_libraries(mlcoalsimXmpi zf_log)
target_link_libraries(mlcoalsimXmpi ${MPI_LIBRARIES} -DinMPI)
install(TARGETS mlcoalsimXmpi DESTINATION ${CMAKE_INSTALL_PREFIX})

# With MPI / With ZnS
add_executable(mlcoalsimXmpi_ZnS ${SOURCE_FILES})
set_target_properties(mlcoalsimX_ZnS PROPERTIES COMPILE_FLAGS "-DZNS_ACTIVE")
target_link_libraries(mlcoalsimXmpi_ZnS zf_log)
target_link_libraries(mlcoalsimXmpi_ZnS ${MPI_LIBRARIES} -DinMPI)
install(TARGETS mlcoalsimXmpi_ZnS DESTINATION ${CMAKE_INSTALL_PREFIX})
