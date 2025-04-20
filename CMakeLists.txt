cmake_minimum_required(VERSION 3.16)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# uncomment to do full static linking including glibc
# usually not recommended
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
set(BUILD_SHARED_LIBS OFF)

project(FlameFractalRenderer VERSION 0.1 LANGUAGES CXX)
include_directories(src)

# boost (link statically)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
set(Boost_USE_STATIC_LIBS ON)

# png (link statically)
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
set(PNG_STATIC TRUE)

# shared library for fractal renderer
add_library(shared_util STATIC src/rng/flame_rng.cpp src/utils/image.cpp src/utils/json.cpp)
target_link_libraries(shared_util PRIVATE PNG::PNG)

# ffr-buf
add_executable(ffr-buf src/ffr_buf.cpp)
target_link_libraries(ffr-buf PRIVATE Boost::program_options shared_util)

# ffr-img
add_executable(ffr-img src/ffr_img.cpp)
target_link_libraries(ffr-img PRIVATE Boost::program_options shared_util)

# cmake commands (run in /build)
# cmake .. -DCMAKE_BUILD_TYPE=Release
# make VERBOSE=1

# finding files automatically to build a target
# add_library(library_name)
# file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")
# target_sources(library_name PRIVATE ${SOURCES})
