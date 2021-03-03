//
// Created by xzc on 20-7-29.
//

#include "awi_wpe_cfg.h"

void awi_wpe_cfg_init(awi_wpe_cfg_t *cfg, int channels)
{
    cfg->alpha_psd     = 0.9999;
    cfg->taps          = 3;
    cfg->predict_delay = 2;
    cfg->channels      = channels;
}