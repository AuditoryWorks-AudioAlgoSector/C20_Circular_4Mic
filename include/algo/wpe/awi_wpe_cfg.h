//
// Created by xzc on 20-7-29.
//

#ifndef ASP_AWI_WPE_CFG_H
#define ASP_AWI_WPE_CFG_H

typedef struct wpe_cfg
{
    int channels;
    int predict_delay;
    int taps;
    float alpha_psd;
}awi_wpe_cfg_t;



void awi_wpe_cfg_init(awi_wpe_cfg_t *cfg, int channels);

#endif //ASP_AWI_WPE_CFG_H
