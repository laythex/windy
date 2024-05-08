cmake_minimum_required(VERSION 3.22)

set(PROJECT_NAME windy)
project(${PROJECT_NAME})

find_package(DOLFIN REQUIRED)
find_package(Boost COMPONENTS filesystem REQUIRED)

include(${DOLFIN_USE_FILE})

file(GLOB SOURCES *.cpp)
add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} dolfin)
target_link_libraries(${PROJECT_NAME} mshr)
target_link_libraries(${PROJECT_NAME} Boost::filesystem)
