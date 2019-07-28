cmake_minimum_required(VERSION 3.10.2)
project(filius)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CXX_DBG "-g")
set(CXX_FLG "--warn-unused-variable -unknown-pragmas -Wall -Werror")
set(CXX_THD "-pthread")
set(CXX_OPT "-Ofast")
set(CXX_MCR "-D__DEBUG__ -D__REGEN__")
list(APPEND CMAKE_CXX_FLAGS "${CXX_DBG} ${CXX_FLAG} ${CXX_THD} ${CXX_OPT} ${CXX_MCR}")

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(G2PP_SRC "${SRC_DIR}/G2PP.cpp")
set(OPTM_SRC "${SRC_DIR}/Optimization.cpp")

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

add_executable(TestSimulation ${SOURCES})

# PHONY TARGETS
add_custom_target(run "./TestSimulation")
add_custom_target(build make -j)
add_custom_target(memcheck memcheck "./TestSimulation")

add_library(G2++ SHARED ${G2PP_SRC})
add_library(OPTM SHARED ${OPTM_SRC})
