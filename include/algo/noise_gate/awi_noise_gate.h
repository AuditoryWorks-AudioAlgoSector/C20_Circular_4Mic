//
// Created by xzc on 19-7-13.
//

#ifndef ASP_AWI_NOISE_GATE_H
#define ASP_AWI_NOISE_GATE_H

#include "awi_noise_gate_cfg.h"

typedef struct awi_noise_gate
{
    awi_noise_gate_cfg_t cfg;
    float recur_gain;
    float recur_gain_prev;
    int attack_hold_counter;
    int release_hold_counter;
}awi_noise_gate_t;

void awi_noise_gate_init(awi_noise_gate_t *p);

void awi_noise_gate_process(awi_noise_gate_t *p, float *audio);

#endif //ASP_AWI_NOISE_GATE_H
