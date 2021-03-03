//
// Created by xzc on 20-4-8.
//

#ifndef ASP_AWI_BEAM_FUSION_CFG_H
#define ASP_AWI_BEAM_FUSION_CFG_H

typedef struct awi_bf_cfg
{
    float alpha_psd_at;
    float alpha_psd_ar;
    float alpha_weight_at;
    float alpha_weight_ar;
}awi_bf_cfg_t;


void awi_bf_cfg_init(awi_bf_cfg_t *p);

#endif //ASP_AWI_BEAM_FUSION_CFG_H
