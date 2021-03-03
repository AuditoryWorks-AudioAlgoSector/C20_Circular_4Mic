//
// Created by xzc on 20-4-8.
//

#include "awi_beam_fusion_cfg.h"

void awi_bf_cfg_init(awi_bf_cfg_t *p)
{
    p->alpha_psd_at = 0.95;
    p->alpha_psd_ar = 0.97;
    p->alpha_weight_at = 0.985;
    p->alpha_weight_ar = 0.991;
}

