#include "awi_simple_vad_cfg.h"


void awi_simple_vad_cfg_init(awi_simple_vad_cfg_t *cfg)
{
    cfg->subband_snr_ratio      = 1.6f;
    cfg->fullband_snr_ratio     = 1.2589f;
    cfg->snr_counter_th         = 5;
    cfg->speech_counter_th      = 3;
    cfg->non_speech_counter_th  = 25;
    cfg->nn_speech_counter_th      = 0;
    cfg->nn_non_speech_counter_th  = 30;

    cfg->alpha_valid_bkg_fast   = 0.00f;
    cfg->alpha_valid_bkg_slow   = 0.9995f;
    cfg->alpha_valid_bkg        = 0.9995f;
    cfg->alpha_global_snr       = 0.95;
    cfg->alpha_local_snr        = 0.90;
    cfg->bkg_energy_high_th     = 1e-6f;
};

