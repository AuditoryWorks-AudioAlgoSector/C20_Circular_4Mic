cmake_minimum_required(VERSION 2.8)

# 项目名称
set(PROJECT_NAME "RNNoise16")
project(${PROJECT_NAME})
set(CMAKE_CXX_STANDARD 14)

set(LINUX_PLATFORM 1)

if(LINUX_PLATFORM)
MESSAGE(STATUS "This is intel build~~~~!")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-g -ggdb -O0 -Wall -fPIC") #-D_LARGEFILE64_SOURCE")
else()

MESSAGE(STATUS "This is arm build~~~~!")
set(CMAKE_FIND_ROOT_PATH "/opt/hisi-linux/x86-arm/arm-himix200-linux")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_C_COMPILER "/opt/hisi-linux/x86-arm/arm-himix200-linux/bin/arm-himix200-linux-gcc")
set(CMAKE_CXX_COMPILER "/opt/hisi-linux/x86-arm/arm-himix200-linux/bin/arm-himix200-linux-g++")
set(CMAKE_ASM_COMPILER "/opt/hisi-linux/x86-arm/arm-himix200-linux/bin/arm-himix200-linux-gcc")
set(CMAKE_CROSS "/opt/hisi-linux/x86-arm/arm-himix200-linux/bin/arm-himix200-linux-")
set(CMAKE_C_FLAGS "-Wall -O2 -D_REENTRANT -mcpu=cortex-a53 -mfloat-abi=softfp -mfpu=neon-vfpv4 -mno-unaligned-access -fno-aggressive-loop-optimizations -ffunction-sections -fdata-sections -DAPP_GLIBC")
set(CMAKE_CXX_FLAGS "-Wall -O2 -D_REENTRANT -mcpu=cortex-a53 -mfloat-abi=softfp -mfpu=neon-vfpv4 -mno-unaligned-access -fno-aggressive-loop-optimizations  -ffunction-sections -fdata-sections -DDAPP_GLIBC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-gc-sections")
endif()


#add_definitions("-DBAND_EXPAND=0")
add_definitions("-DHEXMS_FRAME=1")
#add_definitions("-DUSE_SNR=1")
#add_definitions("-DUSE_LOW_PASS=1")
#add_definitions("-DTRAINING=1")
add_definitions("-DDR_WAV_IMPLEMENTATION=1")


if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    MESSAGE(STATUS "current platform: Linux ")
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows")
    add_definitions("-D_MBCS=1")
    MESSAGE(STATUS "current platform: Windows")
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
    MESSAGE(STATUS "current platform: FreeBSD")
ELSE ()
    MESSAGE(STATUS "other platform: ${CMAKE_SYSTEM_NAME}")
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")

# 统一构建目标输出目录
# 将所有构建目标输出目录定位在构建目录下，按类型分类子目录
#     静态库输出目录 => 构建目录/lib
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
#     动态库输出目录 => 构建目录/bin
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
#     可执行程序输出目录 => 构建目录/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
#统一输出库格式
set(CMAKE_DEBUG_POSTFIX "_d")

# 添加头文件搜索路径，这样其他构建目标可以引用其头文件
include_directories(include)
link_directories(${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
# ---------- 库 ----------
# staticlib 库
add_subdirectory(src)
# webrtc 库
add_subdirectory(webrtc-native-code)
# rnnoise_agc_wrapper
add_subdirectory(wrapper)
# nn_agc_wrapper
add_subdirectory(agc_demo)

# ---------- 应用程序 ----------
# testapp 程序
add_subdirectory(examples)

