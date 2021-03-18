#ifndef _AWI_NSE_H_
#define _AWI_NSE_H_

#include "awi_filterbankparams.h"
#include "awi_nse_cfg.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct awi_nse
{
    awi_nse_cfg_t cfg;

    int init_counter;
    int init_complete;
    float avg_inst_noisy_psd[AWI_FRAME_BAND_COUNT];
    float bkg_est_psd[AWI_FRAME_BAND_COUNT];
    float recur_block_psd[AWI_FRAME_BAND_COUNT];
    float recur_block_psd_prev[AWI_FRAME_BAND_COUNT * REVERB_DELAYS];

    float first_moment_recur_psd[AWI_FRAME_BAND_COUNT];
    float second_moment_recur_psd[AWI_FRAME_BAND_COUNT];
    float inv_Qeq[AWI_FRAME_BAND_COUNT];
    float Qeq_inv_global[AWI_FRAME_BAND_COUNT];
    float global_bias_origin[AWI_FRAME_BAND_COUNT];
    float Qeq_inv_local[AWI_FRAME_BAND_COUNT];
    float local_bias_origin[AWI_FRAME_BAND_COUNT];

    float local_ms_win[AWI_FRAME_BAND_COUNT * AWI_BKG_V];
    float global_ms_win[AWI_FRAME_BAND_COUNT * AWI_BKG_U];
    float local_ms_ns_est[AWI_FRAME_BAND_COUNT];
    float global_ms_ns_est[AWI_FRAME_BAND_COUNT];

} awi_nse_t;

void awi_nse_init(awi_nse_t *p);

void awi_nse_process_multichannel(awi_nse_t *p, float *input_sp, int auxi_vad_flag, int enable_vad_auxi);

void awi_nse_process_mono(awi_nse_t *p, float *input_sp, int auxi_vad_flag, int enable_vad_auxi);

void awi_nse_process_local_func(awi_nse_t *p, int auxi_vad_flag, int enable_vad_auxi);


#ifdef __cplusplus
};
#endif

#endif
