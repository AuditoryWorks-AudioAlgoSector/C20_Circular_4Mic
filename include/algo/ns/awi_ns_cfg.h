#ifndef AWI_NS_CFG_H
#define AWI_NS_CFG_H

#include "awi_constant.h"
#include "string.h"

#define VK_TABLE1_LEN 100
#define VK_TABLE2_LEN 46
#define VK_TABLE3_LEN 46
#define VK_TABLE4_LEN 46

#define VK_TABLE1_START  1.000000000000e-05f
#define VK_TABLE1_END   0.00496000000000000f
#define VK_TABLE2_START VK_TABLE1_END
#define VK_TABLE2_END   0.0499600000000000f
#define VK_TABLE3_START VK_TABLE2_END 
#define VK_TABLE3_END   0.499960000000000f
#define VK_TABLE4_START VK_TABLE3_END
#define VK_TABLE4_END   4.99996000000000f

extern float EXP_INT_TABLE1[VK_TABLE1_LEN];
extern float EXP_INT_TABLE2[VK_TABLE2_LEN];
extern float EXP_INT_TABLE3[VK_TABLE3_LEN];
extern float EXP_INT_TABLE4[VK_TABLE4_LEN];
extern float mask_band_freq_st[AWI_FREQ_BANDS];

typedef struct awi_ns_cfg
{
    float vk_table_step1;
    float vk_table_step2;
    float vk_table_step3;
    float vk_table_step4;

    float gain_min;                           /* minimum gain in linear scale      */
    float alpha_priori_snr;
    float alpha_recur_priori_snr;
    float snr_min;
    float snr_max;
    float avg_recur_priori_snr_prev;
    float avg_recur_priori_snr_max;
    float avg_recur_priori_snr_min;
    float avg_recur_prior_snr_peak;
    float norm_recur_priori_snr_max;
    float norm_recur_priori_snr_min;
    float max_speech_absence_prob;           /* maximum speech absence probability */
    float high_gain_th;
    float low_counter_ratio;
    float high_counter_ratio;
    int low_freq_id;
    int high_freq_id;
    float alpha_sub_gain;
    float high_band_avg_gain_th;
    float alpha_reverb;
    float alpha_recur_gain;
    float delta_reverb_level;

} awi_ns_cfg_t;

void awi_ns_cfg_init(awi_ns_cfg_t *p);

#endif // AWI_NS_CFG_H
