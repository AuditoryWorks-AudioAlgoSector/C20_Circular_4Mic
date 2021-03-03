#ifndef _AWI_CNG_CFG_H_
#define _AWI_CNG_CFG_H_

#include "awi_filterbankparams.h"

extern float Comfort_Noise_Psd_Ref[AWI_FRAME_BAND_COUNT+1];

typedef struct awi_cng_cfg
{
    float noise_sp_scale;
    float fullband_bkg_energy_th_low;
    float fullband_bkg_energy_th_high;
    float alpha_cn_psd_inc;
    float alpha_cn_psd_dec;
    float initial_fullband_cn_energy;
    float cn_energy_ratio;
    float input_energy_th;
    float change_rate_th;
    float alpha_energy;
    float noise_scale;
} awi_cng_cfg_t;

void awi_cng_cfg_init(awi_cng_cfg_t *cfg);

#endif
