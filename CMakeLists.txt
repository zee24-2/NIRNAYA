cmake_minimum_required(VERSION 3.20)
project(SovereignSolver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-O3 -fno-fast-math -Wall -Wextra)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(sov_c_api STATIC src/driver/c_api.cpp)

add_executable(sovsolve src/driver/cli_main.cpp)
target_link_libraries(sovsolve PRIVATE sov_c_api)

add_executable(checker tools/checker.cpp)
add_executable(gen_refinery_models tools/gen_refinery_models.cpp)
add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests PRIVATE sov_c_api)
