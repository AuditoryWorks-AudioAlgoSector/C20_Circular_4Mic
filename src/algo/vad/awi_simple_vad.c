#include "awi_simple_vad.h"
#include "awi_utils.h"
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void awi_simple_vad_init(awi_simple_vad_t *p)
{
    p->speech_counter        = 0;
    p->is_speech_triggered   = 0;
    p->non_speech_counter    = 0;
    p->fullband_snr          = 0;
    p->nn_speech_counter     = 0;
    p->nn_non_speech_counter = 0;
    p->nn_post_vad_flag      = 0;

    p->recur_fullband_bkg_energy = 0;

    memset(p->critical_band_snr, 0, AWI_FREQ_BANDS * sizeof(float));
    memset(p->bkg_est_psd_floor, 0, AWI_FRAME_BAND_COUNT * sizeof(float));

}


void awi_simple_vad_process(awi_simple_vad_t *p, float *bkg_est_psd, float *recur_block_psd, float *inst_noisy_psd)
{
    float subband_bkg_energy[AWI_FREQ_BANDS], subband_sig_energy[AWI_FREQ_BANDS];
    int snr_counter = 0;
    int start_freq_id = 0, end_freq_id = 0;
    float fullband_bkg_energy, fullband_sig_energy;
    float inst_critical_band_snr, inst_fullband_snr;
    float alpha_global_snr_1 = 1 - p->cfg.alpha_local_snr;
    float alpha_local_snr_1 = 1 - p->cfg.alpha_local_snr;


    fullband_bkg_energy = 0, fullband_sig_energy = 0;
    for (int band_id = 0; band_id < AWI_FREQ_BANDS; band_id++)
    {
        start_freq_id = AWI_FREQ_LOW[band_id];
        end_freq_id = AWI_MIN(AWI_FREQ_HIGH[band_id]-1, AWI_FRAME_BAND_COUNT);

        if( AWI_FREQ_LOW[band_id] > AWI_FRAME_BAND_COUNT )
            break;

        subband_bkg_energy[band_id] = 0;
        subband_sig_energy[band_id] = 0;
        for (int freq_id = start_freq_id - 1; freq_id < end_freq_id; freq_id++)
        {
            subband_bkg_energy[band_id] += bkg_est_psd[freq_id];
            subband_sig_energy[band_id] += recur_block_psd[freq_id];
        }
        fullband_bkg_energy += subband_bkg_energy[band_id];
        fullband_sig_energy += subband_sig_energy[band_id];

        inst_critical_band_snr = ( subband_sig_energy[band_id] + AWI_EPS ) / ( subband_bkg_energy[band_id] + AWI_EPS );

        p->critical_band_snr[band_id] = p->cfg.alpha_local_snr * p->critical_band_snr[band_id] + alpha_local_snr_1 * inst_critical_band_snr;

        p->critical_band_snr[band_id] = AWI_MIN(p->critical_band_snr[band_id], 10);
        p->critical_band_snr[band_id] = AWI_MAX(p->critical_band_snr[band_id], 1);

        if ( band_id >= 3 && p->critical_band_snr[band_id] >= p->cfg.subband_snr_ratio )
            snr_counter++;
    }

    inst_fullband_snr = ( fullband_sig_energy + AWI_EPS ) / ( fullband_bkg_energy + AWI_EPS );
    p->fullband_snr = p->cfg.alpha_global_snr * p->fullband_snr + alpha_global_snr_1 * inst_fullband_snr;
    p->fullband_snr = AWI_MIN(p->fullband_snr, 10);
    p->fullband_snr = AWI_MAX(p->fullband_snr, 1);
    p->recur_fullband_bkg_energy = fullband_bkg_energy;

    if ( snr_counter >= p->cfg.snr_counter_th && p->fullband_snr > p->cfg.fullband_snr_ratio )
    {
        p->speech_counter++;
        if ( p->speech_counter >= p->cfg.speech_counter_th )
        {
            p->is_speech_triggered = 1;
            p->non_speech_counter = 0;
        }
    }
    else
    {
        p->non_speech_counter++;
        if ( p->non_speech_counter >= p->cfg.non_speech_counter_th )
        {
            p->is_speech_triggered = 0;
            p->speech_counter = 0;
        }
    }

    awi_simple_vad_compute_bkg_floor(p, bkg_est_psd, inst_noisy_psd, fullband_bkg_energy, fullband_sig_energy);

}


void awi_simple_vad_compute_bkg_floor(awi_simple_vad_t *p, float *bkg_est_psd,
                                      float *inst_noisy_psd, float fullband_bkg_energy, float fullband_sig_energy)
{
    float alpha_valid_bkg, alpha_valid_bkg_1;
    float tmp_freq_bkg_psd;
    float bkg_floor_energy = 0;
    
    if ( p->recur_fullband_bkg_energy < p->cfg.bkg_energy_high_th )
    {
        if ( !p->is_speech_triggered )
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                tmp_freq_bkg_psd = bkg_est_psd[band_id];
                alpha_valid_bkg = p->cfg.alpha_valid_bkg;
                alpha_valid_bkg_1 = 1 - alpha_valid_bkg;
                p->bkg_est_psd_floor[band_id] = alpha_valid_bkg * p->bkg_est_psd_floor[band_id] +
                                                alpha_valid_bkg_1 * tmp_freq_bkg_psd;
            }
        }
        else
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                tmp_freq_bkg_psd = bkg_est_psd[band_id];
                alpha_valid_bkg = p->cfg.alpha_valid_bkg_fast;
                alpha_valid_bkg_1 = 1 - alpha_valid_bkg;
                if ( tmp_freq_bkg_psd > p->bkg_est_psd_floor[band_id])
                {
                    alpha_valid_bkg = p->cfg.alpha_valid_bkg_slow;
                    alpha_valid_bkg_1 = 1 - p->cfg.alpha_valid_bkg_slow;
                }
                p->bkg_est_psd_floor[band_id] = alpha_valid_bkg * p->bkg_est_psd_floor[band_id] +
                                                alpha_valid_bkg_1 * tmp_freq_bkg_psd;
            }
        }     
    }

    for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        bkg_floor_energy += p->bkg_est_psd_floor[band_id];
    }

}



