#ifndef _SRC_MODULE_BEAMFORMING_AWI_FBF_H_
#define _SRC_MODULE_BEAMFORMING_AWI_FBF_H_

#include "awi_fbf_cfg.h"
#include "awi_constant.h"
#include "awi_filterbankparams.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct awi_fbf awi_fbf_t;

struct awi_fbf
{
    awi_fbf_cfg_t cfg;
};

void awi_fbf_init(awi_fbf_t *fbf);
void awi_fbf_process_single(awi_fbf_t *fbf, float *frame_in, float *frame_out, int beam);


#ifdef __cplusplus
};
#endif
#endif
