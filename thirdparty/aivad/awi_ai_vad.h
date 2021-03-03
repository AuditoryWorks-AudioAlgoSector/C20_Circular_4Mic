#ifndef AWI_AI_VAD_H
#define AWI_AI_VAD_H

#include "awi_ai_vad_cfg.h"
#include "forward.h"
#include "mfcc.h"
#include "ringbuf.h"

typedef struct awi_ai_vad
{
    awi_ai_vad_cfg_t *cfg;
    rb_s16_t temp_inputdata_frames;
    rb_t num_frames_mfcc;
    rb_t one_frames_input;
    rb_t pred_array;
    float *temp_feature_input;
    float *delta;
    float *features;
    float *mfcc_features;
    float *input_frame;
    float *powspectrum;
    float *local_data;
    float *total_mfcc_out;
} awi_ai_vad_t;

#ifdef __cplusplus
extern "C"{
#endif

void awi_ai_vad_init(awi_ai_vad_t *ai_vad, awi_ai_vad_cfg_t *cfg);

short awi_ai_vad_process(awi_ai_vad_t *ai_vad, short *input,short *output);

void input_to_mfcc(awi_ai_vad_t *ptr,float *input);

void get_delta(awi_ai_vad_t *ptr, float *input, float *delta, float *temp_feature_input);

void concatenate(awi_ai_vad_t *ptr,float *mfcc,float *delta,float *features);

void mfcc_to_pred(awi_ai_vad_t *ptr,float *step_pred);

void get_mfcc(awi_ai_vad_t *ptr ,float  *input, float *mfcc_out);

//store data to list and get data from list
short awi_ai_vad_access(awi_ai_vad_t *ai_vad, short *input, short *output);

//awi determin
short awi_ai_vad_determine(awi_ai_vad_t *ai_vad, float *input);

float awi_ai_vad_determine_f(awi_ai_vad_t *ai_vad, float *input);

void awi_ai_vad_deinit(awi_ai_vad_t *ai_vad);

#ifdef __cplusplus
};
#endif
#endif //AWI_AI_VAD_H
