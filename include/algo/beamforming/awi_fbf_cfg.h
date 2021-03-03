#ifndef _SRC_MODULE_BEAMFORMING_AWI_FBF_CFG_H_
#define _SRC_MODULE_BEAMFORMING_AWI_FBF_CFG_H_

#include "awi_constant.h"
#include "awi_filterbankparams.h"

#include "fbf_params.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct awi_fbf_cfg awi_fbf_cfg_t;

struct awi_fbf_cfg
{
    int block_angle[AWI_MIC_CHANNEL * AWI_ABF_REFERENCE];
    float *mvdr_fbf_weight;  
};

void awi_fbf_cfg_init(awi_fbf_cfg_t *cfg);

#ifdef __cplusplus
};
#endif

#endif
