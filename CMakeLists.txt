cmake_minimum_required(VERSION 3.10)
project(mach2h5)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_VERBOSE_MAKEFILE ON)


set(CMAKE_CXX_STANDARD 17)

add_library(blosc_filter src/blosc_filter.c include/blosc_filter.h)
add_library(zstd_filter src/zstd_h5plugin.c include/zstd/zstd.h include/zstd/zstd_h5plugin.h)

add_executable(mach2h5
        include/EigenH5.h
        src/main.cpp
        include/cxxopts.hpp
        include/highfive/highfive.hpp
        include/highfive/H5Filter.hpp
        )

target_include_directories(mach2h5 PRIVATE include)
target_include_directories(zstd_filter PRIVATE include)
target_include_directories(blosc_filter PRIVATE include)
find_package(ZSTD)
find_package(BLOSC)

target_include_directories(zstd_filter PUBLIC ${ZSTD_INCLUDE_DIR})
target_include_directories(blosc_filter PUBLIC ${BLOSC_INCLUDE_DIR})

target_link_libraries (blosc_filter ${BLOSC_LIBRARY})
target_link_libraries (zstd_filter ${ZSTD_LIBRARY})

find_package(HDF5)
find_package(Boost REQUIRED COMPONENTS system iostreams)

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#add_executable (example example.cpp)
target_link_libraries (mach2h5 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} blosc_filter zstd_filter ${BLOSC_LIBRARY} ${ZSTD_LIBRARY})


target_link_libraries (mach2h5 ${Boost_LIBRARIES})
target_link_libraries (mach2h5 stdc++fs)
target_link_libraries (mach2h5 ${BLOSC_LIBRARY})
target_link_libraries (mach2h5 ${ZSTD_LIBRARY})


target_include_directories(mach2h5 PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(mach2h5 PRIVATE ${ZSTD_INCLUDE_DIR})
target_include_directories(mach2h5 PRIVATE ${BLOSC_INCLUDE_DIR})
target_include_directories(mach2h5 PRIVATE ${BOOST_INCLUDEDIR})



