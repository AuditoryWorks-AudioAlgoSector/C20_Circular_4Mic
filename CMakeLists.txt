cmake_minimum_required(VERSION 3.5)
project(asp)

set(CMAKE_CXX_STANDARD 11)
set(DEBUG 1)

if(DEBUG)
    set(CMAKE_BUILD_TYPE DEBUG)
else(DEBUG)
    set(CMAKE_BUILD_TYPE RELEASE)
endif(DEBUG)

message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYSTEM_LINUX)
endif()
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYSTEM_DARWIN)
endif()

add_definitions(-DHAVE_ATLAS)
add_definitions(-Wno-deprecated-declarations)
add_definitions(-fPIC)
add_definitions(-Wall)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")

SET(CHANNEL_CONVERT
        src/core/channel-convert/awi_channel_convert.c)

SET(BF
        src/algo/beam_fusion/awi_beam_fusion_cfg.c
        src/algo/beam_fusion/awi_beam_fusion.c)

SET(COMMON
        src/core/common/awi_filterbankframe.c
        src/core/common/awi_filterbank.c
        src/core/common/awi_utils.c)

SET(FBF
        src/algo/beamforming/awi_fbf.c
        src/algo/beamforming/awi_fbf_cfg.c
        src/algo/beamforming/fbf_params.c
        src/algo/beamforming/awi_abf.c
        src/algo/beamforming/awi_abf_cfg.c)


SET(MODULE
        src/awi_algo.c
        src/algo/aec/awi_subband_aec.c
        src/algo/aec/awi_subband_aec_cfg.c
        src/algo/aes/awi_aes.c
        src/algo/aes/awi_aes_cfg.c
        src/algo/cng/awi_cng.c
        src/algo/cng/awi_cng_cfg.c
        src/algo/nse/awi_nse_cfg.c
        src/algo/nse/awi_nse.c
        src/algo/ns/awi_ns_cfg.c
        src/algo/ns/awi_ns.c
        src/algo/agc/awi_agc_cfg.c
        src/algo/agc/awi_agc.c
        src/algo/vad/awi_simple_vad_cfg.c
        src/algo/vad/awi_simple_vad.c
        src/algo/limiter/awi_limiter_cfg.c
        src/algo/limiter/awi_limiter.c
        src/algo/noise_gate/awi_noise_gate_cfg.c
        src/algo/noise_gate/awi_noise_gate.c
        src/algo/pre_de_emphasis/awi_emphasis_cfg.c
        src/algo/pre_de_emphasis/awi_emphasis.c
        src/algo/wpe/awi_eigen_wpe.cpp
        src/algo/wpe/awi_wpe_cfg.c
        )

SET(FFT
        src/core/fft/awi_fft.c
        src/core/fft/kiss_fft.c)

SET(FILTER_BANK
        src/core/filter-bank/awi_dft_analysis_filter_bank.c
        src/core/filter-bank/awi_dft_synthesis_filter_bank.c)

SET(AIVAD
        ./thirdparty/aivad/agc.c 
        ./thirdparty/aivad/awi_ai_vad.c 
        ./thirdparty/aivad/awi_ai_vad_cfg.c 
        ./thirdparty/aivad/forward.c 
        ./thirdparty/aivad/mfcc.c 
        ./thirdparty/aivad/params.c 
        ./thirdparty/aivad/ringbuf.c 
        ./thirdparty/aivad/nn_vad_agc.c)

SET(RNNOISE_HEXFRAME
        ./thirdparty/rnnoise_hexframe/src/kiss_fft_opus.c
        ./thirdparty/rnnoise_hexframe/src/celt_lpc.c
        ./thirdparty/rnnoise_hexframe/src/denoise.c
        ./thirdparty/rnnoise_hexframe/src/pitch.c
        ./thirdparty/rnnoise_hexframe/src/rnn_data.c
        ./thirdparty/rnnoise_hexframe/src/rnn_reader.c
        ./thirdparty/rnnoise_hexframe/src/rnn.c)

ADD_LIBRARY(aspmodule SHARED
        ${COMMON}
        ${FILTER_BANK}
        ${FFT}
        ${MODULE}
        ${CHANNEL_CONVERT}
        ${BF}
        ${CORE}
        ${FBF}
        ${RNNOISE_HEXFRAME}
        ${AIVAD}
        )

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_executable(test_pipeline test/base/test_pipeline.c)
    add_executable(test_algo test/base/test_algo.c)
endif()



include_directories("include/")
include_directories("include/core/")
include_directories("include/core/channel-convert")
include_directories("include/core/fft")
include_directories("include/core/common")
include_directories("include/core/filter-bank")
include_directories("include/algo/beamforming")
include_directories("include/algo/ns")
include_directories("include/algo/doa")
include_directories("include/algo/bs")
include_directories("include/algo/beam_fusion")
include_directories("include/algo/nse")
include_directories("include/algo/aec")
include_directories("include/algo/aes")
include_directories("include/algo/agc")
include_directories("include/algo/cng")
include_directories("include/algo/vad")
include_directories("include/algo/limiter")
include_directories("include/algo/noise_gate")
include_directories("include/algo/pre_de_emphasis")
include_directories("/usr/include")
include_directories("thirdparty/fftw/include")
include_directories("thirdparty/aivad")
include_directories("thirdparty/rnnoise_hexframe/include")
include_directories("include/algo/wpe")
include_directories("thirdparty/eigen3")

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

    target_link_libraries(test_pipeline aspmodule 
                        ${CMAKE_SOURCE_DIR}/thirdparty/fftw/lib/libfftw3f.so
                        /usr/lib/x86_64-linux-gnu/libsndfile.so)
    
    target_link_libraries(test_algo aspmodule
                        ${CMAKE_SOURCE_DIR}/thirdparty/fftw/lib/libfftw3f.so
                        /usr/lib/x86_64-linux-gnu/libsndfile.so)

endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    target_link_libraries(aspmodule /usr/local/lib/libgsl.dylib)
    target_link_libraries(aspmodule pthread)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Android")
    target_link_libraries(aspmodule stdc++)
endif()
