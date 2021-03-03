#编译选项
set(RNNOISE_STATIC_LIB "rnnoise")
#set(RNNOISE_STATIC_LIB "rnnoise_static")
#set(RNNOISE_DYNAMIC_LIB "rnnoise_dynamic")

aux_source_directory(. DIR_RNN_SRCS)
#add_library(${RNNOISE_DYNAMIC_LIB} SHARED ${DIR_RNN_SRCS})
add_library(${RNNOISE_STATIC_LIB} STATIC ${DIR_RNN_SRCS})