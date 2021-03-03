//
// Created by xzc on 19-7-1.
//

#ifndef ASP_AWI_LIMITER_H
#define ASP_AWI_LIMITER_H

#include "awi_limiter_cfg.h"

typedef struct awi_limiter
{
    awi_limiter_cfg_t cfg;
    float at_hold_counter;
    float ar_hold_counter;
}awi_limiter_t;

void awi_limiter_init(awi_limiter_t *p);

void awi_limiter_process(awi_limiter_t *p, float *audio, int is_speech);

#endif //ASP_AWI_LIMITER_H
