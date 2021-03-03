#ifndef RNNOISE16_RNNOISE_AGC_H
#define RNNOISE16_RNNOISE_AGC_H


 #ifdef __cplusplus
 extern "C"{
 #endif

    typedef struct NnvadAgcState NnvadAgcState;

    NnvadAgcState* awi_nnvad_agc_create();

    int awi_nnvad_agc_init(NnvadAgcState* st);

    int awi_nnvad_agc_process(NnvadAgcState* st, float* input, short *output);

    int awi_nnvad_agc_free(NnvadAgcState* st);

    int awi_nnvad_determine(NnvadAgcState* st, float* input);

    float awi_nnvad_determine_f(NnvadAgcState* st, float* input);

    int awi_nnvad_access(NnvadAgcState* st, float* input, short *output);

 #ifdef __cplusplus
 };
 #endif



#endif //RNNOISE16_RNNOISE_AGC_H

