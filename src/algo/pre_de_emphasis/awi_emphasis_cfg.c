//
// Created by xzc on 19-7-26.
//

#include "awi_emphasis_cfg.h"

void awi_emphasis_cfg_init(awi_emphasis_cfg_t *p)
{
    p->alpha_pre    = 0.90;
    p->alpha_de     = 0.90;
    p->enable_flag  = 0;
}
