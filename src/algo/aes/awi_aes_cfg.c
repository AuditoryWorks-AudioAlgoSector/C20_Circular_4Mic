#include "awi_aes_cfg.h"
#include "awi_constant.h"
#include <math.h>

void awi_aes_cfg_init(awi_aes_cfg_t *cfg)
{
  cfg->alpha_psd               = 0.95;
  cfg->alpha_echo_at           = 0.0;
  cfg->alpha_echo_ar           = 0.92;
  cfg->head_length             = 300 * AWI_AEC_CHANNELS;
  cfg->spk_head_energy_th      = 1e-7;
  cfg->spk_head_freq_psd_th    = 5e-9;
  cfg->sub_aes_gain_floor      = 0.001f;
  cfg->psy_relax_factor        = 3;     /* 1 ~ 3.1623 */
  cfg->mid_band_id_ey          = 96;
  cfg->spk_energy_th_high      = 1e-4f;
  cfg->full_band_energy_th     = 5e-9;
  cfg->alpha_energy            = 0.92;
  cfg->high_band_id_th         = 120;
  cfg->high_band_gain_cnt      = 0.01;
  cfg->alpha_fullband          = 0.9;
  cfg->mid_fullband_gain       = 0.1;
  cfg->slope_fullband_gain     = 20;
  cfg->weak_aes_gain_th        = 0.1;
  cfg->hl_energy_ratio_change  = 20;
  cfg->alpha_sub_erle_low      = 0.10;
  cfg->alpha_sub_erle_high     = 0.99;
  cfg->alpha_sub_erl           = 0.10;
}


float SpreadFunc[AWI_FREQ_BANDS-2] =
{ 0.999680209586331,
  0.371021214256294,
  0.058438075989367,
  0.007246709261353,
  0.000819990912477,
  0.000088810696868,
  0.000009389461625,
  0.000000978243099,
  0.000000100948526,
  0.000000010348986,
  0.000000001055964,
  0.000000000107370,
  0.000000000010888,
  0.000000000001102,
  0.000000000000111,
  0.000000000000011,
  0.000000000000001,
  0.000000000000000,
  0.000000000000000,
  0.000000000000000,
  0.000000000000000,
  0.000000000000000
  };




