#include "awi_cng.h"
#include "awi_fft.h"
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "awi_utils.h"


void awi_cng_init(awi_cng_t *p)
{
    p->rand_counts              = 0;
    p->is_full                  = 0;
    p->is_bkg_psd_valid         = 0;
    p->init_flag                = 1;
    p->pure_noise_flag               = 1;
    p->fullband_bkg_energy_prev = -10;
    p->recur_input_energy       = -10;

    memset(p->bkg_psd_est_mat, 0, CNG_RAN_NUM * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->rand_seq_sp, 0, AWI_FB_FRAME_LENGTH * 2 * sizeof(float));
    memset(p->recur_cn_psd_used, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
}

void awi_cng_process(awi_cng_t *p, float *subband_aes_gain, float *input_sp, float *out_sp, 
                     float *bkg_est_psd_floor, float bkg_energy)
{
    float fullband_bkg_energy = 0, inst_fullband_sig_energy = 0;
    int rand_seq;
    float *current_cn_psd_used;
    float scale;
    float sp_re, sp_im;
    float fullband_ref_cn_energy = Comfort_Noise_Psd_Ref[AWI_FRAME_BAND_COUNT];
    float added_fullband_cn_energy = 0;
    float decision_factor = 0;
    float change_rate = 0;
    int idx_re, idx_im;
    float alpha_cn_psd, alpha_cn_psd_1;

    for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        fullband_bkg_energy += bkg_est_psd_floor[band_id];
        sp_re = input_sp[idx_re];
        sp_im = input_sp[idx_im];
        inst_fullband_sig_energy += ( sp_re * sp_re + sp_im * sp_im );
    }
    if ( p->fullband_bkg_energy_prev < 0 )
        p->fullband_bkg_energy_prev = fullband_bkg_energy;

    if ( inst_fullband_sig_energy > p->recur_input_energy )
        p->recur_input_energy = inst_fullband_sig_energy;
    else
    {
        p->recur_input_energy = p->cfg.alpha_energy * p->recur_input_energy +
                ( 1 - p->cfg.alpha_energy ) * inst_fullband_sig_energy;
    }


    if ( !p->is_bkg_psd_valid )
    {
        if ( fullband_bkg_energy > p->cfg.fullband_bkg_energy_th_low && fullband_bkg_energy < p->cfg.fullband_bkg_energy_th_high )
        {
            p->is_bkg_psd_valid = 1;
            memcpy(p->bkg_psd_est_mat, bkg_est_psd_floor, sizeof(float) * AWI_FRAME_BAND_COUNT);
            p->rand_counts++;
        }
    }
    else
    {
        if ( fullband_bkg_energy > p->cfg.fullband_bkg_energy_th_low && fullband_bkg_energy < p->cfg.fullband_bkg_energy_th_high )
        {
            change_rate = fullband_bkg_energy / ( p->fullband_bkg_energy_prev + AWI_EPS );
            if ( change_rate < p->cfg.change_rate_th )
            {
                float *bkg_psd_est_start = p->bkg_psd_est_mat + p->rand_counts * AWI_FRAME_BAND_COUNT;
                memcpy(bkg_psd_est_start, bkg_est_psd_floor, sizeof(float) * AWI_FRAME_BAND_COUNT);
                p->rand_counts++;
            }
            p->is_bkg_psd_valid = 1;
        }
    }
    p->fullband_bkg_energy_prev = fullband_bkg_energy;


    if (p->rand_counts == CNG_RAN_NUM)
    {
        p->rand_counts = 0;
        p->is_full = 1;
    }


    for (int band_id = 0; band_id < AWI_FB_FRAME_LENGTH; band_id++)
    {
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        p->rand_seq_sp[idx_re] = gaussrand();
        p->rand_seq_sp[idx_im] = 0;
    }
#ifdef AWI_ARM
    awi_fft_256(p->rand_seq_sp);
#else
    awi_fft_1d(p->rand_seq_sp, p->rand_seq_sp, AWI_FB_FRAME_LENGTH);
#endif


    if ( p->is_bkg_psd_valid )
    {
        if( p->is_full != 1 )
        {
            rand_seq = rand() % p->rand_counts;
        }
        else
        {
            rand_seq = rand() % CNG_RAN_NUM;
        }
        current_cn_psd_used = p->bkg_psd_est_mat + rand_seq * AWI_FRAME_BAND_COUNT;

        if ( p->init_flag )
        {
            memcpy(p->recur_cn_psd_used, current_cn_psd_used, AWI_FRAME_BAND_COUNT * sizeof(float));
            p->init_flag = 0;
        }
        else
        {
            for ( int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                if (current_cn_psd_used[band_id] > p->cfg.change_rate_th * p->recur_cn_psd_used[band_id])
                    alpha_cn_psd = p->cfg.alpha_cn_psd_inc;
                else
                    alpha_cn_psd = p->cfg.alpha_cn_psd_dec;
                alpha_cn_psd_1 = 1 - alpha_cn_psd;

                p->recur_cn_psd_used[band_id] = alpha_cn_psd * p->recur_cn_psd_used[band_id] +
                        alpha_cn_psd_1 * current_cn_psd_used[band_id];
            }
        }

        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            idx_re = band_id + band_id;
            idx_im = idx_re + 1;
            scale = ( 1 - subband_aes_gain[band_id] ) *
                    p->cfg.noise_sp_scale * sqrtf(p->recur_cn_psd_used[band_id]);
            sp_re = scale * p->rand_seq_sp[idx_re];
            sp_im = scale * p->rand_seq_sp[idx_im];
            added_fullband_cn_energy += sp_re * sp_re + sp_im * sp_im;

            out_sp[idx_re]  = input_sp[idx_re] + sp_re;
            out_sp[idx_im]  = input_sp[idx_im] + sp_im;
        }
    }
    else
    {
        added_fullband_cn_energy = AWI_MIN(p->cfg.initial_fullband_cn_energy, bkg_energy);

        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            idx_re = band_id + band_id;
            idx_im = idx_re + 1;
            scale = ( 1 - subband_aes_gain[band_id] ) * sqrtf(added_fullband_cn_energy
                    / fullband_ref_cn_energy * Comfort_Noise_Psd_Ref[band_id]) * p->cfg.noise_sp_scale;
            out_sp[idx_re] = input_sp[idx_re] + scale * p->rand_seq_sp[idx_re];
            out_sp[idx_im] = input_sp[idx_re] + scale * p->rand_seq_sp[idx_im];
        }
    }

    decision_factor = ( p->recur_input_energy + 1e-20f ) / ( added_fullband_cn_energy + 1e-20f );

    if ( decision_factor > p->cfg.cn_energy_ratio && p->recur_input_energy > p->cfg.input_energy_th )
        p->pure_noise_flag = 1;
    else
        p->pure_noise_flag = 0;

}

void awi_cng_process_mod(awi_cng_t *p, float *input_sp, float *out_sp, int vad_flag, 
                     float *bkg_est_psd_floor, float bkg_energy)
{
    int idx_re, idx_im;
    float fullband_ref_cn_energy = Comfort_Noise_Psd_Ref[AWI_FRAME_BAND_COUNT];
    float scale;
    
    memcpy(out_sp, input_sp, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    

    for (int band_id = 0; band_id < AWI_FB_FRAME_LENGTH; band_id++)
    {
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        p->rand_seq_sp[idx_re] = gaussrand();
        p->rand_seq_sp[idx_im] = 0;
    }

#ifdef AWI_ARM
    awi_fft_256(p->rand_seq_sp);
#else
    awi_fft_1d(p->rand_seq_sp, p->rand_seq_sp, AWI_FB_FRAME_LENGTH);
#endif

    if ( p->is_bkg_psd_valid )
    {
        if ( !vad_flag )
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;
                scale = p->cfg.noise_sp_scale * sqrtf(p->recur_cn_psd_used[band_id]) / 2;
                out_sp[idx_re] = p->cfg.noise_scale * input_sp[idx_re] + scale * p->rand_seq_sp[idx_re];
                out_sp[idx_im] = p->cfg.noise_scale * input_sp[idx_im] + scale * p->rand_seq_sp[idx_im];
            }
        }
    }
    else
    {
        float applied_cn_energy = AWI_MIN(p->cfg.initial_fullband_cn_energy, bkg_energy);
        applied_cn_energy = AWI_MAX(applied_cn_energy, p->cfg.fullband_bkg_energy_th_low);
        if ( !vad_flag )
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;
                scale = sqrtf(applied_cn_energy / fullband_ref_cn_energy * 
                        Comfort_Noise_Psd_Ref[band_id]) * p->cfg.noise_sp_scale;
                out_sp[idx_re] = p->cfg.noise_scale * input_sp[idx_re] + scale * p->rand_seq_sp[idx_re];
                out_sp[idx_im] = p->cfg.noise_scale * input_sp[idx_im] + scale * p->rand_seq_sp[idx_im];
            }
        }
    }
}



float gaussrand()
{
    static float U = 0.f, V = 0.f, m1 = 0.f, m2 = 0.f;
    static int phase = 0;
    float Z = 0.f;

    if( phase == 0 )
    {
        U = rand() / (RAND_MAX + 1.0f);

        V = rand() / (RAND_MAX + 1.0f);

        if (U == 0)
            U = 1.0f / RAND_MAX;

        m1 = sqrtf(-2.0f * logf(U));
        m2 = 2.0f * MATH_PI * V;
        Z = m1 * sinf(m2);
    }
    else
    {
        Z = m1 * cosf(m2);
    }

    phase = 1 - phase;

    return Z;
}

