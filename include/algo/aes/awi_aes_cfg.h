#ifndef _AWI_AES_CFG_H_
#define _AWI_AES_CFG_H_

#include <stdbool.h>
#include "awi_ns_cfg.h"

extern float SpreadFunc[AWI_FREQ_BANDS-2];


typedef struct awi_aes_cfg
{
    float alpha_psd;
    float alpha_echo_at;
    float alpha_echo_ar;
    int head_length;
    float spk_head_energy_th;
    float spk_head_freq_psd_th;
    float sub_aes_gain_floor;
    float psy_relax_factor;
    int mid_band_id_ey;
    float spk_energy_th_high;
    float full_band_energy_th;
    float alpha_energy;
    int high_band_id_th;
    float high_band_gain_cnt;
    float alpha_fullband;
    float mid_fullband_gain;
    float slope_fullband_gain;
    float weak_aes_gain_th;
    float hl_energy_ratio_change;
    float alpha_sub_erle_low;
    float alpha_sub_erle_high;
    float alpha_sub_erl;

} awi_aes_cfg_t;


void awi_aes_cfg_init(awi_aes_cfg_t *cfg);

#endif
