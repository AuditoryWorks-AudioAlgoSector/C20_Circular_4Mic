#ifndef _AWI_SUBBAND_AEC_H_
#define _AWI_SUBBAND_AEC_H_

#include "awi_subband_aec_cfg.h"
#include "awi_filterbankframe.h"
#include "awi_constant.h"


typedef struct awi_subband_aec_part1
{
    awi_subband_aec_cfg_t cfg;

    int is_aec_on;
    int is_aec_update;
    int divergence_flag_fast;

    float sum_inst_Far_autoPsd;
    float sum_recur_Far_autoPsd;
    float sum_aec_Far_autoPsd;
    float max_spk_psd;
    int header_counter;
    int is_head;
    int init_aec_flag;

    int init_global_far_psd;
    int init_enable_aec_psd;

    float recur_full_band_erle_chs[AWI_BEAM_CHANNEL];
    float fullband_spk_mic_erl_chs[AWI_BEAM_CHANNEL];
    float recur_aec_eo_fast_energy[AWI_BEAM_CHANNEL];
    float recur_aec_eo_slow_energy[AWI_BEAM_CHANNEL];

    float xfm[AWI_SPK_CHANNEL * AWI_AEC_TAP * AWI_FRAME_BAND_COUNT * 2];
    float wfb[AWI_BEAM_CHANNEL * TOTAL_TAP_NUM * 2];
    float yfm[AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * 2];
    float nlms_mEk[AWI_FRAME_BAND_COUNT * 2];
    float sye[AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * 2];
    float sy[AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT];
    float se[AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT];
    float aec_Far_autoPsd[AWI_FRAME_BAND_COUNT];
    float recur_Far_autoPsd[AWI_FRAME_BAND_COUNT];
    float inst_Far_autoPsd[AWI_FRAME_BAND_COUNT];
    float aec_Near_autoPsd[AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT];
    float tmp_yo_slow[AWI_FRAME_BAND_COUNT * 2];
    float max_spk_freq_psd[AWI_FRAME_BAND_COUNT];

} awi_subband_aec_part1_t;


typedef struct awi_subband_aec_part2
{
    /* for dual filter */
    float wfb_fast[AWI_BEAM_CHANNEL * TOTAL_TAP_NUM * 2];
}awi_subband_aec_part2_t;

void awi_subband_aec_init(awi_subband_aec_part1_t *p1, awi_subband_aec_part2_t *p2);

void awi_subband_aec_process(awi_subband_aec_part1_t *p1, awi_subband_aec_part2_t *p2, float *fbFrameIn, float *refFrameIn, float *fbFrameOut);

void awi_subband_aec_filter(float *input, float *weight, float *output);

void awi_subband_aec_update_weight(float *input, float *weight, float *nlms_eo);

void awi_subband_aec_compute_step_init(awi_subband_aec_part1_t *p, float *aec_eo, float *tmp_yo, float *tmp_mic_sp, float *tmp_mic_psd,
        float *tmp_ch_se, float *tmp_ch_sy, float *tmp_ch_ye, int ch_id);

void awi_subband_aec_compute_step(awi_subband_aec_part1_t *p, float *aec_eo, float *tmp_yo, float *tmp_mic_sp, float *tmp_mic_psd,
        float *tmp_ch_se, float *tmp_ch_sy, float *tmp_ch_ye, int ch_id);
#endif

