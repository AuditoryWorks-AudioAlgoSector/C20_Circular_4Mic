//
// Created by xzc on 19-7-1.
//

#include "awi_limiter_cfg.h"

void awi_limiter_cfg_init(awi_limiter_cfg_t *p)
{
    p->limiter_am      = 0.9999;
    p->desired_am      = 0.35;
    p->alpha_at        = 0.0;
    p->alpha_ar        = 0.99;
    p->peak_am         = 0;
    p->gain            = 1;
    p->amp_step_inc    = 1;
    p->max_amplifier   = 40;
    p->recur_amplifier = 40;
    p->alpha_amp_dec   = 0.9823;
    p->alpha_amp_inc   = 0.9975;
    p->at_hold_th      = 128;
    p->ar_hold_th      = 0;
    p->amp_low_th      = 2;
    p->amp_high_th     = 50;
    p->max_amp_inc     = 1.01;
    p->max_amp_dec     = 0.99;
}
