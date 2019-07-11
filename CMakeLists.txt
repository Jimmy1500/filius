cmake_minimum_required(VERSION 3.10.2)
project(filius)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

add_executable(TestSimulation ${SOURCES})
target_link_libraries(TestSimulation Threads::Threads)
