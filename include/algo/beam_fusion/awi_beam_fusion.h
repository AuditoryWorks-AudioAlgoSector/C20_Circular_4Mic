//
// Created by xzc on 20-4-8.
//

#ifndef ASP_AWI_BEAM_FUSION_H
#define ASP_AWI_BEAM_FUSION_H

#include "awi_beam_fusion_cfg.h"
#include "awi_filterbankparams.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif    /* __cplusplus */

typedef struct awi_bf
{
    awi_bf_cfg_t cfg;
    bool complex_weighting;
    int direction_id;
    float inst_auto_psd[AWI_FRAME_BAND_COUNT*AWI_BEAM_CHANNEL];
    float recur_auto_psd[AWI_FRAME_BAND_COUNT*AWI_BEAM_CHANNEL];
    float recur_energy[AWI_BEAM_CHANNEL];
    float beam_weight[AWI_BEAM_CHANNEL];
    float recur_beam_weight[AWI_BEAM_CHANNEL];
}awi_bf_t;

void awi_bf_init(awi_bf_t *p, bool complex_weight);

void awi_bf_process(awi_bf_t *p, float *ref_sp, float *input_sp, float *output_sp);

#ifdef __cplusplus
}
#endif   /* __cplusplus  */

#endif //ASP_AWI_BEAM_FUSION_H
