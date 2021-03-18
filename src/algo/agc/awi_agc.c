#include "awi_agc.h"
#include "awi_filterbankparams.h"
#include "awi_utils.h"
#include "awi_aes.h"
#include "awi_ns_cfg.h"
#include <string.h>


void awi_agc_init(awi_agc_t *agc)
{
    agc->recur_energy_dB    = -115.f;

    agc->high_band_energy_ratio = 0.1;
    
    memset(agc->inst_agc_in_psd, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);

    memset(agc->agc_gain, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);

}



float dl_sigmoid(float in, float mid_point, float slope)
{
    float exp_val;
    exp_val = expf( -( in - mid_point ) / slope );
    return 1.f / ( 1.f + exp_val );
}


void awi_agc_process(awi_agc_t *p, awi_aes_t *aes, float *input_sp, float *out_sp, float fullband_snr, float *critical_band_snr,
        int is_speech, int agc_enable_status, float *weighted_factor)
{
    float summed_inst_agc_in_psd_dB;
    float adj_factor;
    float inst_applied_gain_dB, mOffsetGain;
    float avg_agc_gain;
    float alpha_freq_gain, alpha_freq_gain_1;
    float weighted_adj, multi_aes_gain, aes_gain_adj;
    float high_band_crt_factor;

    summed_inst_agc_in_psd_dB = awi_agc_pre_process(p, input_sp);

    adj_factor = dl_sigmoid(fullband_snr, p->cfg.sigmoid_mid_point, p->cfg.sigmoid_slope);
    if ( is_speech )
    {
        alpha_freq_gain = p->cfg.alpha_freq_gain_high;
    }
    else
    {
        alpha_freq_gain = p->cfg.alpha_freq_gain_low;
    }

    inst_applied_gain_dB = compute_inst_agc_gain_in_dB(p, summed_inst_agc_in_psd_dB, is_speech);
    inst_applied_gain_dB *= agc_enable_status;
    inst_applied_gain_dB = AWI_MAX(AWI_MIN(inst_applied_gain_dB, p->cfg.applied_gain_max_dB),
                p->cfg.agc_gain_floor_high);

    if ( inst_applied_gain_dB > p->cfg.applied_gain_dB )
    {
        p->cfg.applied_gain_dB = p->cfg.alpha_inc * p->cfg.applied_gain_dB +
                                  ( 1 - p->cfg.alpha_inc ) * inst_applied_gain_dB;
    }
    else
    {
        p->cfg.applied_gain_dB = p->cfg.alpha_dec * p->cfg.applied_gain_dB +
                                  ( 1 - p->cfg.alpha_dec ) * inst_applied_gain_dB;
    }

    mOffsetGain = powf(10, p->cfg.applied_gain_dB / 20 ) * adj_factor * agc_enable_status;

    avg_agc_gain =  compute_inst_average_agc_gain_in_linear_scale(p, fullband_snr, critical_band_snr, mOffsetGain, is_speech);

    alpha_freq_gain_1 = 1 - alpha_freq_gain;
    for(int band = p->cfg.low_cutoff_freq_id; band < AWI_FRAME_BAND_COUNT; band++)
    {
        p->agc_gain[band] = alpha_freq_gain * p->agc_gain[band] + alpha_freq_gain_1 * avg_agc_gain;
    }

    for (int band = p->cfg.low_cutoff_freq_id; band < p->cfg.mid_freq_id; band++)
    {
        weighted_adj = AWI_MIN(1, weighted_factor[band] / p->cfg.max_weighted_factor);
        multi_aes_gain = aes->min_fullband_aes_gain * aes->min_aes_gain[band];
        aes_gain_adj = AWI_MIN(1, 2 * multi_aes_gain / ( multi_aes_gain + 0.5f ));
        p->agc_gain[band] = p->agc_gain[band] * weighted_adj * aes_gain_adj;
    }

    high_band_crt_factor = AWI_MIN(p->cfg.high2full_ratio_th / p->high_band_energy_ratio, 1);
    for (int band = p->cfg.mid_freq_id; band < AWI_FRAME_BAND_COUNT; band++)
    {
        weighted_adj = high_band_crt_factor * AWI_MIN(1, weighted_factor[band] / p->cfg.max_weighted_factor);
        multi_aes_gain = aes->min_fullband_aes_gain * aes->min_aes_gain[band];
        aes_gain_adj = AWI_MIN(1, 2 * multi_aes_gain / ( multi_aes_gain + 0.75f ));
        p->agc_gain[band] = p->agc_gain[band] * weighted_adj * aes_gain_adj;
    }

    // for(int band = p->cfg.low_cutoff_freq_id; band < AWI_FRAME_BAND_COUNT; band++)
    // {
    //     if ( p->inst_agc_in_psd[band] < 1e-9f )
    //     {
    //         p->agc_gain[band] = p->cfg.agc_gain_floor_high;
    //     }
    // }

    apply_agc_gain(input_sp, out_sp, p->agc_gain, p->cfg.low_cutoff_freq_id, p->cfg.agc_gain_floor_low,
                   p->cfg.agc_gain_floor_high);
    
}


float awi_agc_pre_process(awi_agc_t *p, float *input_sp)
{
    float sp_re, sp_im;
    int idx_re, idx_im;
    float summed_inst_agc_in_psd, summed_inst_agc_in_psd_dB;
    float low_band_energy, high_band_energy;

    // compute instant auto psd of input spectrum
    for(int band = p->cfg.low_cutoff_freq_id; band < AWI_FRAME_BAND_COUNT; band++)
    {
        idx_re = band + band;
        idx_im = idx_re + 1;
        sp_re = input_sp[idx_re];
        sp_im = input_sp[idx_im];
        p->inst_agc_in_psd[band] = sp_re * sp_re + sp_im * sp_im;
    }

    // compute input energy in frequency domain using dB scale corresponding to that in time domain
    low_band_energy = 0;
    high_band_energy = 0;
    for(int band = p->cfg.low_cutoff_freq_id; band < p->cfg.mid_freq_id; band++)
    {
        low_band_energy += p->inst_agc_in_psd[band];
    }

    for(int band = p->cfg.mid_freq_id; band < AWI_FRAME_BAND_COUNT - 1; band++)
    {
        high_band_energy += p->inst_agc_in_psd[band];
    }

    summed_inst_agc_in_psd = low_band_energy + high_band_energy;

    p->high_band_energy_ratio = high_band_energy / (summed_inst_agc_in_psd + AWI_EPS);

    summed_inst_agc_in_psd = summed_inst_agc_in_psd * 2 / AWI_FB_FRAME_LENGTH;

    summed_inst_agc_in_psd_dB = 10 * log10f( summed_inst_agc_in_psd + AWI_EPS );

    // smooth dB energy
    p->recur_energy_dB = p->cfg.alpha_sum * p->recur_energy_dB +
            ( 1 - p->cfg.alpha_sum ) * summed_inst_agc_in_psd_dB;

    return summed_inst_agc_in_psd_dB;
}


float compute_inst_agc_gain_in_dB(awi_agc_t *p, float inst_input_db_energy, int is_speech)
{
    float inst_gaindB_limit;
    float inst_alpha, inst_gain_gap;
    float inst_applied_gain_dB;

    /* compute upper gain limit */
    inst_gaindB_limit = p->cfg.desired_sum_power_dBth - p->recur_energy_dB;
    if ( inst_gaindB_limit > p->cfg.recur_gain_dB_limit_prev )
    {
        inst_alpha = p->cfg.alpha_gain_ar;
    }
    else
    {
        inst_alpha = p->cfg.alpha_gain_at;
    }
    p->cfg.recur_gain_dB_limit = inst_alpha * p->cfg.recur_gain_dB_limit + ( 1 - inst_alpha ) * inst_gaindB_limit;
    p->cfg.recur_gain_dB_limit_prev = p->cfg.recur_gain_dB_limit;

    if ( p->recur_energy_dB < p->cfg.recur_energy_dB_th )
    {
        p->cfg.offset_gain_dB = p->cfg.alpha_offset_dec * p->cfg.offset_gain_dB;
    }
    else if( p->recur_energy_dB + p->cfg.offset_gain_dB < p->cfg.desired_sum_power_dBth && is_speech)
    {
        p->cfg.offset_gain_dB = p->cfg.alpha_offset_inc * p->cfg.offset_gain_dB + \
                                 ( 1 - p->cfg.alpha_offset_inc ) * AWI_MAX(p->cfg.offset_gain_dB_inc, p->cfg.offset_gain_dB + p->cfg.offset_gain_step);
    } else
    {
        p->cfg.offset_gain_dB = p->cfg.alpha_offset_dec * p->cfg.offset_gain_dB + \
                                 ( 1 - p->cfg.alpha_offset_dec ) * p->cfg.offset_gain_dB_dec;
    }

    /* compute gap between gradual gain and its recursive limit in dB */
    inst_gain_gap = p->cfg.offset_gain_dB - p->cfg.recur_gain_dB_limit;
    inst_gain_gap = AWI_MAX(inst_gain_gap, 0);
    if ( inst_gain_gap < p->cfg.recur_gain_gap_prev )
    {
        inst_alpha = p->cfg.alpha_gain_at;
    }
    else
    {
        inst_alpha = p->cfg.alpha_gain_ar;
    }
    p->cfg.recur_gain_gap = inst_alpha * p->cfg.recur_gain_gap + ( 1 - inst_alpha ) * inst_gain_gap;
    p->cfg.recur_gain_gap_prev = p->cfg.recur_gain_gap;

    /* compute instant agc gain in dB */
    inst_applied_gain_dB = p->cfg.offset_gain_dB - p->cfg.recur_gain_gap;

    /* impose constraint for total output energy of agc */
    if ( inst_input_db_energy + inst_applied_gain_dB > p->cfg.sum_power_th )
        inst_applied_gain_dB = AWI_MIN(AWI_MAX(p->cfg.sum_power_th - inst_input_db_energy - 1, 0),
                                       inst_applied_gain_dB - 1);

    return inst_applied_gain_dB;
}


float compute_inst_average_agc_gain_in_linear_scale(awi_agc_t *p, float fullband_snr,
        float *critical_band_snr, float fullband_gain, int is_speech)
{
    float critical_band_gain, subband_scale;
    float avg_agc_gain = 0;
    for(int band_id = 0; band_id < AWI_FREQ_BANDS; band_id++)
    {
        if ( mask_band_freq_st[band_id] < 8000 )
        {
            int upper = AWI_MIN(AWI_FREQ_HIGH[band_id]-1, AWI_FRAME_BAND_COUNT);


            if ( critical_band_snr[band_id] > p->cfg.m_subband_snrTh && fullband_snr > p->cfg.m_fullband_snrTh && is_speech)
            {
                critical_band_gain = fullband_gain;
            }
            else
            {

                subband_scale = critical_band_snr[band_id] / p->cfg.max_critical_band_snr;
                subband_scale = AWI_MIN(subband_scale, 1);

                critical_band_gain = subband_scale * fullband_gain;

            }

            if ( band_id >= p->cfg.low_cutoff_critical_band_id )
                avg_agc_gain += critical_band_gain;

            for(int k = AWI_FREQ_LOW[band_id] - 1; k < upper; k++)
            {
                p->agc_gain[k] = critical_band_gain;
            }
        }
    }

    avg_agc_gain /= ( AWI_FREQ_BANDS_VALID - p->cfg.low_cutoff_critical_band_id );

    return avg_agc_gain;

}

void apply_agc_gain(float *input_sp, float *out_sp, float *applied_agc_gain, int low_freq_id_eof, float lower_gain_th, float higher_gain_th)
{
    float agc_gain_var;
    float sp_re, sp_im;
    int idx_re, idx_im;
    for(int band = 0; band < low_freq_id_eof * 2; band++)
        out_sp[band] = input_sp[band] * lower_gain_th;

    for(int band = low_freq_id_eof; band < AWI_FRAME_BAND_COUNT; band++)
    {
        idx_re = band + band;
        idx_im = idx_re + 1;
        agc_gain_var = AWI_MAX(applied_agc_gain[band], higher_gain_th);
        sp_re = input_sp[idx_re] * agc_gain_var;
        sp_im = input_sp[idx_im] * agc_gain_var;
        out_sp[idx_re] = sp_re;
        out_sp[idx_im] = sp_im;
    }
}

