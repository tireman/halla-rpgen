cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project (NpolRates)
message("Configuring project NpolRates.")

set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

# include ROOT
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
message("    -> ROOT Include directory is: ${ROOT_INCLUDE_DIR}")
message("    -> ROOT Library Directory is: ${ROOT_LIBRARY_DIR}")
include_directories(${ROOT_INCLUDE_DIR})

file(COPY ${PROJECT_SOURCE_DIR}/scripts DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

add_executable(NpolRates ${sources})
target_link_libraries(NpolRates NpolClasses ${ROOT_LIBRARIES})
