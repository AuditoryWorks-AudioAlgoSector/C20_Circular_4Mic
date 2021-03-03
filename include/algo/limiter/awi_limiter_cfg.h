//
// Created by xzc on 19-7-1.
//

#ifndef ASP_AWI_LIMITER_CFG_H
#define ASP_AWI_LIMITER_CFG_H

typedef struct awi_limiter_cfg
{
    float limiter_am;
    float desired_am;
    float alpha_at;
    float alpha_ar;
    float peak_am;
    float gain;
    float amp_step_inc;
    float max_amplifier;
    float recur_amplifier;
    float at_hold_th;
    float ar_hold_th;
    float alpha_amp_dec;
    float alpha_amp_inc;
    float amp_low_th;
    float amp_high_th;
    float max_amp_inc;
    float max_amp_dec;
}awi_limiter_cfg_t;

void awi_limiter_cfg_init(awi_limiter_cfg_t *p);

#endif //ASP_AWI_LIMITER_CFG_H
