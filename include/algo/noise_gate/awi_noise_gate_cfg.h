//
// Created by xzc on 19-7-13.
//

#ifndef ASP_AWI_NOISE_GATE_CFG_H
#define ASP_AWI_NOISE_GATE_CFG_H

typedef struct awi_noise_gate_cfg
{
    float threshold;
    float alpha_at;
    float alpha_ar;
    float scale_down;
    int attack_hold_samples;
    int release_hold_samples;

}awi_noise_gate_cfg_t;


void awi_noise_gate_cfg_init(awi_noise_gate_cfg_t *p);

#endif //ASP_AWI_NOISE_GATE_CFG_H
