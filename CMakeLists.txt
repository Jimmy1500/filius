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

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

add_executable(TestSimulation ${SOURCES})
