#include "awi_agc_cfg.h"
#include "awi_aes_cfg.h"


void awi_agc_cfg_init(awi_agc_cfg_t *p)
{
    p->low_cutoff_freq_id                   = 2;
    p->low_cutoff_critical_band_id          = 4;
    p->agc_gain_floor_low                   = 0.01f;
    p->agc_gain_floor_high                  = 0.4467f;
    p->sigmoid_mid_point                    = 1.2589f;
    p->sigmoid_slope                        = 0.1f;
    p->max_critical_band_snr                = 10.f;
    p->max_weighted_factor                  = 1.0f;
    p->m_subband_snrTh                      = 1.5819f;
    p->m_fullband_snrTh                     = 1.2589f;
    p->sum_power_th                         = -70.f;
    p->offset_gain_dB                       = 0.f;
    p->recur_energy_dB_th                   = -108.f;
    p->alpha_offset_inc                     = 0.90f;
    p->alpha_offset_dec                     = 0.98f;
    p->applied_gain_max_dB                  = 10.f;
    p->applied_gain_min_dB                  = -6.f;
    p->applied_gain_dB                      = 0;
    p->offset_gain_dB_inc                   = 16.f;
    p->offset_gain_dB_dec                   = -8.f;
    p->alpha_sum                            = 0.90f;
    p->offset_gain_step                     = 3.0f;
    p->desired_sum_power_dBth               = -80.f;
    p->alpha_freq_gain_high                 = 0.30f;
    p->alpha_freq_gain_low                  = 0.05f;
    p->alpha_inc                            = 0.90f;
    p->alpha_dec                            = 0.10f;
    p->mid_freq_id                          = 64;
    p->high2full_ratio_th                   = 0.3;

    p->recur_gain_dB_limit                  = 0;
    p->recur_gain_dB_limit_prev             = 0;
    p->alpha_gain_at                        = 0.0f;
    p->alpha_gain_ar                        = 0.5f;
    p->recur_gain_gap                       = 0;
    p->recur_gain_gap_prev                  = 0;

}
