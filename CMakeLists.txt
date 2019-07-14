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

add_library(G2++ SHARED ${G2PP_SRC})
add_library(OPT SHARED ${OPTM_SRC})

add_executable(TestSimulation ${SOURCES})

link_directories(${CMAKE_BINARY_DIR})
set (PROJECT_LINK_LIBS ${PROJECT_BINARY_DIR}/libOPT.so)
target_link_libraries(TestSimulation ${PROJECT_LINK_LIBS})
