#ifndef _AWI_AES_H_
#define _AWI_AES_H_

#include "awi_constant.h"
#include "awi_filterbankparams.h"
#include "awi_aes_cfg.h"
#include "awi_subband_aec.h"

typedef struct awi_aes
{
    awi_aes_cfg_t cfg;

    int head_counter;
    float min_fullband_aes_gain;
    float min_aes_gain[AWI_FRAME_BAND_COUNT];
    float aes_gain[AWI_FRAME_BAND_COUNT];
    float recur_spk_auto_psd[AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT];
    float recur_fbf_auto_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float recur_aec_eo_auto_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float recur_aec_yo_auto_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float recur_fbf_aec_eo_cross_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float recur_fbf_aec_eo_cross_imag_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float residue_echo_psd[AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT];
    float psy_near_auto_psd[AWI_FRAME_BAND_COUNT];
    float recur_low_band_energy[AWI_AEC_CHANNELS];
    float recur_high_band_energy[AWI_AEC_CHANNELS];
    float recur_gsc_low_band_energy[AWI_AEC_CHANNELS];
    float recur_gsc_high_band_energy[AWI_AEC_CHANNELS];

} awi_aes_t;

void awi_aes_init(awi_aes_t *p);


void awi_aes_process(awi_aes_t *p, awi_subband_aec_part1_t *aec, float *fbf_sp, float *aec_sp, float *gsc_sp,
        float *aes_sp, float *bkg_psd_floor);

void apply_psycho_acoustic_postfilter(float *desired_sig_psd_est, float *residue_echo_psd_est,
                                      float *aes_gain_beam, float psy_relax_factor);

#endif
