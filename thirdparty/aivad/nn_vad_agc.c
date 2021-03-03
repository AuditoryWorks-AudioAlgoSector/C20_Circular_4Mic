#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <memory.h>
#include "awi_ai_vad.h"
#include "awi_ai_vad_cfg.h"
#include "agc.h"
#include "nn_vad_agc.h"
//using namespace webrtc;
//using namespace awi;

#define FRAME_SIZE 256
#define BYTES_PER_FRAME_10MS (FRAME_SIZE*2)
#define SAMPLE_RATE 16000
#define SPLIT_BAND_NUM (3)
#define GAIN_NUM (11)
#define AGC_COMPRESSION_GAIN_DB 40
#define AGC_LIMITER_ENABLE 1
#define AGC_TARGET_LEVEL_DBFS (6)
#define FRAME_TIME_LEN 16
#define DELAY_FRAMES 5
#define DELAY_TIME 80 //80ms
#define DELAY_SAMPLES (DELAY_TIME*16)
#define SHORT_MAX (32767.0f)
#define ALPHA_COEFF 1

struct NnvadAgcState {
    awi_ai_vad_t vad;
    awi_ai_vad_cfg_t vad_cfg;
    void *agcInst;
    //MediaBuffer *inBuffer;
    //MediaBuffer *outBuffer;
    int32_t capture_levels;
    int32_t new_capture_level;
    uint8_t saturation_warning;
    int stream_has_echo;
    int16_t split_band_data[SPLIT_BAND_NUM][FRAME_SIZE];
    int16_t *split_bands[SPLIT_BAND_NUM];
    float split_band_data_f[SPLIT_BAND_NUM][FRAME_SIZE];
    float *split_bands_f[SPLIT_BAND_NUM];
    int32_t gains[GAIN_NUM];
};

int awi_nnvad_agc_init(NnvadAgcState *st) {
    printf("nnvad agc version 0927 15:51\n");
    awi_ai_vad_cfg_init(&st->vad_cfg);
   
    awi_ai_vad_init(&st->vad, &st->vad_cfg);

    int status = WebRtcAgc_Init(st->agcInst, 0, 255, 2, 16000);
    if (status != 0) {
        printf("WebRtcAgc_Init fail\n");
        WebRtcAgc_Free(st->agcInst);
        return -1;
    }

    WebRtcAgcConfig agcConfig;
    agcConfig.compressionGaindB = AGC_COMPRESSION_GAIN_DB;//在Fixed模式下，越大声音越大
    agcConfig.limiterEnable = AGC_LIMITER_ENABLE;
    agcConfig.targetLevelDbfs = AGC_TARGET_LEVEL_DBFS;   //dbfs表示相对于full scale的下降值，0表示full scale，越小声音越大
    status = WebRtcAgc_set_config(st->agcInst, agcConfig);
    if (status != 0) {
        printf("WebRtcAgc_set_config fail\n");
        WebRtcAgc_Free(st->agcInst);
        return -1;
    }
    st->capture_levels = 0;
    st->new_capture_level = 0;
    st->saturation_warning = 0;
    st->stream_has_echo = 0;
    for (int i = 0; i < SPLIT_BAND_NUM; ++i) {
        st->split_bands[i] = st->split_band_data[i];
        st->split_bands_f[i] = st->split_band_data_f[i];
    }
    return 0;
}

NnvadAgcState *awi_nnvad_agc_create() {
    NnvadAgcState *st = (NnvadAgcState*)malloc(sizeof(NnvadAgcState));
    memset(st, 0, sizeof(NnvadAgcState));

    st->agcInst = WebRtcAgc_Create();
    return st;
}

int awi_nnvad_agc_free(NnvadAgcState *st) {
    awi_ai_vad_deinit(&st->vad);
    WebRtcAgc_Free(st->agcInst);
    return 0;
}

int awi_nnvad_agc_process(NnvadAgcState *st, float *input, short *output) {  
    //for (int j = 0; j < FRAME_SIZE; ++j) {
    //    st->split_bands[0][j] = input[j] * SHORT_MAX * ALPHA_COEFF;
    //}
    int nAgcRet = WebRtcAgc_Process(st->agcInst,
                                    &st->vad, 
                                    (const int16_t *const *) st->split_bands,
                                    1,
                                    FRAME_SIZE,
                                    (int16_t *const *) st->split_bands, 
                                    output, 
                                    st->capture_levels,
                                    &st->new_capture_level,
                                    st->stream_has_echo,
                                    &st->saturation_warning);
    return nAgcRet;
}

int awi_nnvad_determine(NnvadAgcState* st, float* input){
    float input_tmp[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; ++i) {
        input_tmp[i] = input[i] * ALPHA_COEFF;
    }
    
    return awi_ai_vad_determine(&st->vad, input_tmp);
}

float awi_nnvad_determine_f(NnvadAgcState* st, float* input){
    float input_tmp[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; ++i) {
        input_tmp[i] = input[i] * ALPHA_COEFF;
    }
    
    return awi_ai_vad_determine_f(&st->vad, input_tmp);
}

int awi_nnvad_access(NnvadAgcState* st, float* input, short *output){
    //input from ns out witch belong to [-1,1],so 
    for (int j = 0; j < FRAME_SIZE; ++j) {
       st->split_bands[0][j] = input[j] * SHORT_MAX * ALPHA_COEFF;
    }
    return awi_ai_vad_access(&st->vad, st->split_bands[0], output);
}
