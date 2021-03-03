#ifndef _AWI_CNG_H_
#define _AWI_CNG_H_
#include <stdbool.h>
#include "awi_cng_cfg.h"
#include "awi_filterbankframe.h"

#define MATH_PI 3.14159f


typedef struct awi_cng
{
    awi_cng_cfg_t cfg;

    int rand_counts;

    bool is_full;

    bool is_bkg_psd_valid;

    bool init_flag;

    int pure_noise_flag;

    float fullband_bkg_energy_prev;

    float recur_input_energy;

    float bkg_psd_est_mat[AWI_FRAME_BAND_COUNT * CNG_RAN_NUM];

    float rand_seq_sp[AWI_FB_FRAME_LENGTH * 2];

    float recur_cn_psd_used[AWI_FRAME_BAND_COUNT];

} awi_cng_t;


void awi_cng_init(awi_cng_t *p);

void awi_cng_process(awi_cng_t *p, float *subband_aes_gain, float *input_sp, float *out_sp, 
                     float *bkg_est_psd_floor, float bkg_energy);

void awi_cng_process_mod(awi_cng_t *p, float *input_sp, float *out_sp, int vad_flag, 
                     float *bkg_est_psd_floor, float bkg_energy);

float gaussrand(void);


#endif

