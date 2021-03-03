#ifndef AWI_AGC_CFG_H
#define AWI_AGC_CFG_H


typedef struct awi_agc_cfg
{
    int low_cutoff_freq_id;
    int low_cutoff_critical_band_id;
    float agc_gain_floor_low;
    float agc_gain_floor_high;
    float sigmoid_mid_point;
    float sigmoid_slope;
    float max_critical_band_snr;
    float max_weighted_factor;
    float m_subband_snrTh;
    float m_fullband_snrTh;
    float sum_power_th;
    float offset_gain_dB;
    float recur_energy_dB_th;
    float alpha_offset_inc;
    float alpha_offset_dec;
    float applied_gain_max_dB;
    float applied_gain_min_dB;
    float applied_gain_dB;
    float offset_gain_dB_inc;
    float offset_gain_dB_dec;
    float alpha_sum;
    float offset_gain_step;
    float desired_sum_power_dBth;
    float alpha_freq_gain_high;
    float alpha_freq_gain_low;
    float alpha_inc;
    float alpha_dec;
    int mid_freq_id;
    float high2full_ratio_th;

    /* for limiter */
    float recur_gain_dB_limit;
    float recur_gain_dB_limit_prev;
    float alpha_gain_at;
    float alpha_gain_ar;
    float recur_gain_gap;
    float recur_gain_gap_prev;

} awi_agc_cfg_t;

void awi_agc_cfg_init(awi_agc_cfg_t *p);


#endif // AWI_AGC_CFG_H

