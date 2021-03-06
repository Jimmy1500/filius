PROJECT(filius CXX)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.2)

SET(CXX_DBG "-g")
SET(CXX_FLG "--warn-unused-variable -unknown-pragmas -Wall -Werror")
SET(CXX_THD "-pthread")
SET(CXX_OPT "-Ofast")
SET(CXX_MCR "-D__DEBUG__ -D__REGEN__")
LIST(APPEND CMAKE_CXX_FLAGS "${CXX_DBG} ${CXX_FLAG} ${CXX_THD} ${CXX_OPT}")

if (MSVC)
elseif(APPLE)
elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wextra")
else()
    message( FATAL_ERROR "OS not recognized!" )
endif()

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
FILE(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
ADD_EXECUTABLE(TestSimulation ${SOURCES})

# PHONY TARGETS
ADD_CUSTOM_TARGET(mk make -j4)
ADD_CUSTOM_TARGET(rn TestSimulation)
ADD_CUSTOM_TARGET(mc memcheck ./TestSimulation)

SET(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
# ADD_LIBRARY(G2++ SHARED ${SRC_DIR}/G2PP.cpp)
# ADD_LIBRARY(OPTM SHARED ${SRC_DIR}/Optimization.cpp)
