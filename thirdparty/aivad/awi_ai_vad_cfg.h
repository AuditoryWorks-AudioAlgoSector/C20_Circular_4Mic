#ifndef AWI_AI_VAD_CFG_H
#define AWI_AI_VAD_CFG_H

#include "params.h"

typedef struct awi_ai_vad_cfg {
    
    float th;
    short vad_flag;
    int tempFrames;
    int num_total_pred;
    int num_mfcc_frames;
    int num_mfcc_features;
    int features_size;
    int nfilt;
    int numcep;
    int ceplifter;
    int nfft;
    int sr;
    int frame_size;
    int frame_step;
    int appendEnergy;
    float preem_tempdata;
    float Preemphasis;

} awi_ai_vad_cfg_t;

#ifdef __cplusplus
extern "C"{
#endif

void awi_ai_vad_cfg_init(awi_ai_vad_cfg_t *cfg);


#ifdef __cplusplus
};
#endif

#endif //AWI_AI_VAD_CFG_H