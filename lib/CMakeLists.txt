add_library(munit SHARED "${CMAKE_CURRENT_SOURCE_DIR}/munit/munit.c")
target_include_directories(munit INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/munit")
