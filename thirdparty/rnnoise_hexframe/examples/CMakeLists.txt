#编译选项
set(RNNOISE_TRAIN "rnnoise_train")
set(RNNOISE_DEMO "rnnoise_demo")
set(RNNOISE_TRAIN_MULTI "rnnoise_train_multi")
set(RNNOISE_TRAIN_CMVN "rnnoise_train_cmvn")
set(WEBRTC_AGC_DEMO "webrtc_agc_demo")
set(RNNOISE_AGC_DEMO "rnnoise_agc_demo")
set(RNNOISE_AGC_MULTI "ring_buffer_test")
set(RNNOISE_AGC_TEST "ragc_exe")
#featur提取目标		
#add_executable(${RNNOISE_TRAIN} rnnoise_train.c)
#target_link_libraries(${RNNOISE_TRAIN} rnnoise_dynamic -lm)

add_executable(${RNNOISE_TRAIN_MULTI} rnnoise_train_multi.c)
target_link_libraries(${RNNOISE_TRAIN_MULTI} rnnoise -lm)

#demo
add_executable(${RNNOISE_DEMO} rnnoise_demo.c)
target_link_libraries(${RNNOISE_DEMO} rnnoise -lm)

#归一化提取feature 控制snr
#add_executable(${RNNOISE_TRAIN_CMVN} rnnoise_train_cmvn.c)
#target_link_libraries(${RNNOISE_TRAIN_CMVN} rnnoise_d -lm)

#webrtc-agc demo
#add_executable(${WEBRTC_AGC_DEMO} webrtc-agc.cpp)
#target_include_directories(${WEBRTC_AGC_DEMO} PUBLIC ../webrtc-native-code)
#target_link_libraries(${WEBRTC_AGC_DEMO} webrtc_apm_d)

#rnnoise+agc demo
#add_executable(${RNNOISE_AGC_DEMO} rnnoise_agc.cpp)
#target_include_directories(${RNNOISE_AGC_DEMO} PUBLIC ../webrtc-native-code)
#target_link_libraries(${RNNOISE_AGC_DEMO} webrtc_apm_d rnnoise_d -lm)

#rnnoise+agc multi thread
#add_executable(${RNNOISE_AGC_MULTI} ring_buffer_test.cpp ../utility/media_buffer.cpp ../utility/media_buffer.h)
#target_include_directories(${RNNOISE_AGC_MULTI} PUBLIC ../webrtc-native-code)
#target_link_libraries(${RNNOISE_AGC_MULTI} webrtc_apm_d rnnoise_d -lm -lpthread)

#ragc lib test
#add_executable(${RNNOISE_AGC_TEST} ragc_test.cpp)
#target_include_directories(${RNNOISE_AGC_TEST} PUBLIC ../wrapper)
#target_link_libraries(${RNNOISE_AGC_TEST} ragc -lm)