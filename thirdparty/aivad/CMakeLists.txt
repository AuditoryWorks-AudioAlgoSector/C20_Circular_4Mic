
#set(FFTW_PATH "/home/xiezc/lcd/fftw-3.3.8/fftw_ubuntu/")
#set(RNNOISE_AGC_DYNAMIC_LIB "ragc")
#set(CMAKE_CXX_FLAGS "-g -ggdb -O3 -Wall -fPIC") #-D_LARGEFILE64_SOURCE")

if(LINUX_PLATFORM)
MESSAGE(STATUS "agc demo This is intel build~~~~!")
set(FFTW_INCLUDE_PATH "/home/xiezc/lcd/fftw-3.3.8/fftw_ubuntu/include")
set(FFTW_LIB_PATH "/home/xiezc/lcd/fftw-3.3.8/fftw_ubuntu/lib")
set(NNVAD_AGC_STATIC_LIB "nnvagc")
else()
MESSAGE(STATUS "agc demo This is arm build~~~~!")
set(FFTW_INCLUDE_PATH "/home/xiezc/lcd/fftw-3.3.8/fftw_arm/include")
set(FFTW_LIB_PATH "/home/xiezc/lcd/fftw-3.3.8/fftw_arm/lib")
set(NNVAD_AGC_STATIC_LIB "nnvagc_arm")
endif()

include_directories(./)

include_directories(../../fftw-3.3.8/fftw_ubuntu/include)

aux_source_directory(../utility DIR_UTILITY_SRCS)

set(AGC_SRCS
        ./agc.c 
        ./awi_ai_vad.c 
        ./awi_ai_vad_cfg.c 
        ./forward.c 
        ./mfcc.c 
        ./nnvad_params.c
        ./ringbuf.c 
        ./nn_vad_agc.cpp
        )

add_library(${NNVAD_AGC_STATIC_LIB} STATIC ${AGC_SRCS} ${DIR_UTILITY_SRCS})

#add_executable(nnagc_demo main.c)
#target_link_libraries(nnagc_demo -lnnvagc_d -lm -lfftw3f)
#target_link_directories(nnagc_demo PUBLIC ${FFTW_LIB_PATH})


add_executable(nnvad_agc_demo nnvad_agc_demo.cpp ${AGC_SRCS})
target_link_libraries(nnvad_agc_demo libfftw3f.a -lm)
target_link_directories(nnvad_agc_demo PUBLIC ${FFTW_LIB_PATH})