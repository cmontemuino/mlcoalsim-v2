set(CMAKE_THREAD_PREFER_PTHREAD TRUE)

set(BINARY ${CMAKE_PROJECT_NAME}_test)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.c)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})

target_link_libraries(${BINARY} munit libstr)

#add_test(NAME str COMMAND ${BINARY} /str)
add_test(NAME ${BINARY} COMMAND ${BINARY})
