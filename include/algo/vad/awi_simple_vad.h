#ifndef _AWI_SIMPLE_VAD_H_
#define _AWI_SIMPLE_VAD_H_
#include <stdbool.h>
#include "awi_filterbankparams.h"
#include "awi_simple_vad_cfg.h"


typedef struct awi_simple_vad
{
    awi_simple_vad_cfg_t cfg;
    bool is_speech_triggered;

    int speech_counter;
    int non_speech_counter;

    int nn_speech_counter;
    int nn_non_speech_counter;
    int nn_post_vad_flag;

    float fullband_snr;
    float recur_fullband_bkg_energy;

    float critical_band_snr[AWI_FREQ_BANDS];
    float bkg_est_psd_floor[AWI_FRAME_BAND_COUNT];

} awi_simple_vad_t;


void awi_simple_vad_init(awi_simple_vad_t *p);

void awi_simple_vad_process(awi_simple_vad_t *p, float *bkg_est_psd, float *recur_block_psd, float *inst_noisy_psd);

void awi_simple_vad_compute_bkg_floor(awi_simple_vad_t *p, float *bkg_est_psd,
        float *inst_noisy_psd, float fullband_bkg_energy, float fullband_sig_energy);



#endif

