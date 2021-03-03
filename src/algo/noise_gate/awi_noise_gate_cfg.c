//
// Created by xzc on 19-7-13.
//

#include "awi_noise_gate_cfg.h"

void awi_noise_gate_cfg_init(awi_noise_gate_cfg_t *p)
{
    p->threshold             = 0.020;
    p->alpha_at              = 0.9983;
    p->alpha_ar              = 0.9662;
    p->scale_down            = 0.7071;
    p->attack_hold_samples   = 512;
    p->release_hold_samples  = 0;
}
