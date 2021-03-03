#ifndef AWI_AGC_H
#define AWI_AGC_H

#include <math.h>
#include "awi_agc_cfg.h"
#include "awi_aes.h"
#include "awi_filterbankparams.h"


typedef struct awi_agc
{
    awi_agc_cfg_t cfg;
    float recur_energy_dB;
    float high_band_energy_ratio;
    float inst_agc_in_psd[AWI_FRAME_BAND_COUNT];
    float agc_gain[AWI_FRAME_BAND_COUNT];

} awi_agc_t;


void awi_agc_init(awi_agc_t *agc);

void awi_agc_process(awi_agc_t *agc, awi_aes_t *aes, float *input_sp, float *out_sp, float fullband_snr, float *critical_band_snr,
        int is_speech, int agc_enable_status, float *weighted_factor);

float awi_agc_pre_process(awi_agc_t *p, float *input_sp);

float compute_inst_agc_gain_in_dB(awi_agc_t *p, float inst_input_db_energy, int is_speech);

float compute_inst_average_agc_gain_in_linear_scale(awi_agc_t *p, float fullband_snr,
        float *critical_band_snr, float fullband_gain, int is_speech);

void apply_agc_gain(float *input_sp, float *out_sp, float *applied_agc_gain, int low_freq_id_eof, float lower_gain_th, float higher_gain_th);

float dl_sigmoid(float in, float mid_point, float slope);

#endif // AWI_AGC_H
