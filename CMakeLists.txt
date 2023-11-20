cmake_minimum_required(VERSION 3.10)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

option(NON_ADD_LIB "disable add any lib" ON)
option(ADD_STATIC_LIB "enable add static lib" OFF)
option(ADD_SHARED_LIB "enable add shared lib" OFF)



# set(OpenCV_DIR "/home/jiao/repo/3rdparty/opencv/opencv-4.7.0/build")


project(sgbm)


message(STATUS "opencv dir:" ${OpenCV_DIR})


message(STATUS "CMAKE_BUILD_TYPE :  ${CMAKE_BUILD_TYPE}")
message(STATUS "Runtime platform: ${CMAKE_SYSTEM_NAME}")
# 编译器路径打印
message(STATUS "C Complier: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Complier: ${CMAKE_CXX_COMPILER}")


if(MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /arch:AVX2 ")
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]" OR CMAKE_C_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif (ANDROID) # arrch64
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -lm -llog")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11 -lm -llog")
else ()
        set(CMAKE_C_FLAGS " -fPIC -O3 -mavx2")
        set(CMAKE_CXX_FLAGS "-fPIC -std=c++11 -O3 -mavx2 -fpermissive")
endif ()

find_package(Threads)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED opencv_core opencv_highgui opencv_imgcodecs opencv_calib3d)

# 打印库信息
message(STATUS "OpenCV_INCLUDE_DIRS:" ${OpenCV_INCLUDE_DIRS})
message(STATUS "OpenCV_LIBS:" ${OpenCV_LIBS})


# 头文件路径
include_directories(${OpenCV_INCLUDE_DIRS}  ${CMAKE_SOURCE_DIR}/inc)


add_subdirectory(src)
add_subdirectory(example)

file(GLOB demo example/main.cpp)
install(FILES  ${demo} DESTINATION demo)