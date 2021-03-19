#ifndef AWI_NS_H
#define AWI_NS_H

#include <string.h>
#include <math.h>
#include "awi_ns_cfg.h"
#include "awi_filterbankparams.h"

typedef struct awi_ns
{
    awi_ns_cfg_t cfg;
    float inst_post_snr_cng[AWI_FRAME_BAND_COUNT];
    float priori_snr[AWI_FRAME_BAND_COUNT];
    float gain_speech[AWI_FRAME_BAND_COUNT];
    float inst_post_snr_prev[AWI_FRAME_BAND_COUNT];
    float recur_priori_snr[AWI_FRAME_BAND_COUNT];
    float local_recur_priori_snr[AWI_FRAME_BAND_COUNT];
    float global_recur_priori_snr[AWI_FRAME_BAND_COUNT];
    float local_speech_freq_prob[AWI_FRAME_BAND_COUNT];
    float global_speech_freq_prob[AWI_FRAME_BAND_COUNT];
    float ns_gain_seq[AWI_FRAME_BAND_COUNT];
    float recur_ns_gain_seq[AWI_FRAME_BAND_COUNT];
    float low_avg_gain[AWI_FREQ_BANDS];
    float high_avg_gain[AWI_FREQ_BANDS];
    float high_gain_stat_counter[AWI_FREQ_BANDS];
    float full_avg_gain[AWI_FREQ_BANDS];
    float reverb_psd_est[AWI_FRAME_BAND_COUNT];
    
} awi_ns_t;

void awi_ns_init(awi_ns_t *p);

void awi_ns_process(awi_ns_t *p, float *input_sp, float *out_sp, float *input_psd, float *recur_block_psd, float* bkg_est_psd, int speech_flag);

#endif // AWI_NS_H
