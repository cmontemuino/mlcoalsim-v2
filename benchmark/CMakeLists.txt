cmake_minimum_required(VERSION 3.14.5)

project(benchmark)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic")

set(HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(HEADERS benchmark.h benchmark_types.h)
set(SOURCES benchmark.c)

add_library(${PROJECT_NAME} ${HEADERS} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${HEADERS_DIR}>)
target_link_libraries(${PROJECT_NAME} zf_log)
