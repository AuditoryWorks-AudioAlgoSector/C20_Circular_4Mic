#ifndef _AWI_ABF_CFG_H_
#define _AWI_ABF_CFG_H_

#include "awi_constant.h"
#include "awi_filterbankparams.h"

typedef struct awi_abf_cfg
{
    
    int aicTap;
    int references;
    int AFlen;

    float vss_div_factor;
    float vss_stat_factor;
    float alp;
    float mufb;
    float midPoint;
    float slope;
    float minSNR;
    float maxSNR;
    float leaking_factor;
    float power_reg_factor;

    float divergence_th;
    float alpha_snr;
    
} awi_abf_cfg_t;


void awi_abf_cfg_init(awi_abf_cfg_t *cfg);

#endif
