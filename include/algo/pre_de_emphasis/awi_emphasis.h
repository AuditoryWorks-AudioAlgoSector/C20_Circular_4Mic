//
// Created by xzc on 19-7-26.
//

#ifndef ASP_AWI_EMPHASIS_H

#include "awi_emphasis_cfg.h"
#include "awi_constant.h"

typedef struct awi_emphasis
{
    awi_emphasis_cfg_t cfg;
    float de_emp_mem;
    float de_emp_first_out;
    float pre_emp_first_in_spk;
    float pre_emp_first_in_chs[AWI_MIC_CHANNEL];
}awi_emphasis_t;

void awi_emphasis_init(awi_emphasis_t *p);

void awi_pre_emphasis_process(awi_emphasis_t *p, float *audio, float *first_in);

void awi_de_emphasis_process(awi_emphasis_t *p, float *audio);

#define ASP_AWI_EMPHASIS_H

#endif //ASP_AWI_EMPHASIS_H
