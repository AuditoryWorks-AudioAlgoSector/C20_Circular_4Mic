#include "awi_filterbankframe.h"
#include "awi_constant.h"
#include "awi_utils.h"
#include "awi_aes.h"
#include <math.h>


void awi_aes_init(awi_aes_t *p)
{
    p->head_counter = 0;
    p->min_fullband_aes_gain = 1000;
    memset(p->min_aes_gain, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->aes_gain, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_spk_auto_psd, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->psy_near_auto_psd, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_fbf_auto_psd, 0, AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_aec_eo_auto_psd, 0, AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_aec_yo_auto_psd, 0, AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_fbf_aec_eo_cross_psd, 0, AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->residue_echo_psd, 0, AWI_AEC_CHANNELS * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p->recur_low_band_energy, 0, AWI_AEC_CHANNELS * sizeof(float));
    memset(p->recur_high_band_energy, 0, AWI_AEC_CHANNELS * sizeof(float));
    memset(p->recur_gsc_low_band_energy, 0, AWI_AEC_CHANNELS * sizeof(float));
    memset(p->recur_gsc_high_band_energy, 0, AWI_AEC_CHANNELS * sizeof(float));
}


float Sigmoid_mod(float val, float midPoint, float slope)
{
    float step  = (val - midPoint) * slope;
    float coeff = 1 - 0.5f * ( 1 + step / (0.1f + fabsf(step)));
    return fmaxf_local(fminf_local(coeff, 0.99), 0.01);
}


void awi_aes_process(awi_aes_t *p, awi_subband_aec_part1_t *aec, float *fbf_sp, float *aec_sp, float *gsc_sp,
                     float *aes_sp, float *bkg_psd_floor)
{
    float *fbf_sp_beam, *aec_eo_sp_beam, *aec_yo_sp_beam, *gsc_sp_beam, *aes_sp_beam;
    float *fbf_auto_psd_beam, *aec_eo_auto_psd_beam, *aec_yo_auto_psd_beam, *fbf_aec_eo_cross_psd_beam;
    float *residue_echo_psd_beam, *aes_gain_beam;
    float sp_re, sp_im, sp_re2, sp_im2;
    float inst_low_band_energy, inst_high_band_energy;
    float max_spk_energy, near_energy, gsc_energy;
    int idx_re, idx_im;
    float gsc_sp_re, gsc_sp_im;
    float tmp_gain, tmp_gain_minus;
    float freq_spk_psd;
    float tmp_residue_echo_psd;
    float alpha_psd_1, alpha_echo_at_1, alpha_echo_ar_1, alpha_energy_1;
    float global_erle_weight, global_erle_weight_inv;
    float fullband_aes_gain, alpha_fullband_1;
    float tmp_psd_fbf, tmp_psd_eo, tmp_psd_yo, tmp_cpsd, tmp_psd_gsc;
    float sub_erle, combined_erle;
    float sub_erl, sub_erl_inv;
    float initial_est_echo_psd;
    float alpha_sub_erle, alpha_full_erle;
    float near_hl_energy_ratio, gsc_hl_energy_ratio, hl_energy_ratio_change;
    float inst_gsc_high_energy, inst_gsc_low_energy;
    float spk_energy_th_high;
    float fullband_erl_inv;
    float inst_post_ser;
    float alpha_global_erl;
    const int sp_offset = AWI_FRAME_BAND_COUNT * 2;
    const int psd_offset = AWI_FRAME_BAND_COUNT;

    max_spk_energy = 0;
    if (aec->is_aec_on)
    {
        max_spk_energy = AWI_MAX(aec->sum_inst_Far_autoPsd, aec->sum_aec_Far_autoPsd);
    }

    p->min_fullband_aes_gain = 1000;
    for(int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        p->min_aes_gain[band_id] = 1000;
    }

    alpha_echo_at_1 = 1 - p->cfg.alpha_echo_at;
    alpha_echo_ar_1 = 1 - p->cfg.alpha_echo_ar;
    alpha_energy_1  = 1 - p->cfg.alpha_energy;
    alpha_global_erl = 1 - p->cfg.alpha_sub_erl;
    alpha_psd_1 = 1 - p->cfg.alpha_psd;

    fbf_sp_beam = fbf_sp;
    aec_eo_sp_beam = aec_sp;
    aec_yo_sp_beam = aec->yfm;
    fbf_auto_psd_beam = p->recur_fbf_auto_psd;
    aec_eo_auto_psd_beam = p->recur_aec_eo_auto_psd;
    aec_yo_auto_psd_beam = p->recur_aec_yo_auto_psd;
    fbf_aec_eo_cross_psd_beam = p->recur_fbf_aec_eo_cross_psd;
    gsc_sp_beam = gsc_sp;
    residue_echo_psd_beam = p->residue_echo_psd;
    aes_sp_beam = aes_sp;
    aes_gain_beam = p->aes_gain;

    for (int beam_id = 0; beam_id < AWI_AEC_CHANNELS; beam_id++)
    {
        // fullband erle
        if (aec->is_aec_on)
        {
            global_erle_weight = AWI_MAX(aec->recur_full_band_erle_chs[beam_id], 1.0f);
            global_erle_weight_inv = 1.f / global_erle_weight;
        }
        else
        {
            global_erle_weight = 0;
            global_erle_weight_inv = 0;
        }
        
        inst_low_band_energy = 0;
        inst_high_band_energy = 0;
        inst_gsc_high_energy = 0;
        inst_gsc_low_energy = 0;

        spk_energy_th_high = p->cfg.spk_energy_th_high;

        // fullband erl (spk ---> fbf)
        if (aec->is_aec_on)
        {
            aec->fullband_spk_mic_erl_chs[beam_id] = AWI_MAX(aec->fullband_spk_mic_erl_chs[beam_id], 1);
            fullband_erl_inv = 1 / ( aec->fullband_spk_mic_erl_chs[beam_id] + AWI_EPS );
        }
        else
        {
            fullband_erl_inv = 0;
        }

        // compute recursive psd for spk reference signal 
        if ( beam_id == 0 )
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                p->recur_spk_auto_psd[band_id] = p->cfg.alpha_psd * p->recur_spk_auto_psd[band_id] +
                                        alpha_psd_1 * aec->inst_Far_autoPsd[band_id];
            }
        }
        
        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            idx_re = band_id + band_id;
            idx_im = idx_re + 1;

            // compute psd for gsc output signal 
            gsc_sp_re = gsc_sp_beam[idx_re];
            gsc_sp_im = gsc_sp_beam[idx_im];
            tmp_psd_gsc = gsc_sp_re * gsc_sp_re + gsc_sp_im * gsc_sp_im + AWI_EPS;

            if ( aec->max_spk_freq_psd[band_id] > 1e-12f && aec->sum_aec_Far_autoPsd > aec->cfg.spk_energy_th )
            {

                // compute psd for fbf output signal 
                sp_re = fbf_sp_beam[idx_re];
                sp_im = fbf_sp_beam[idx_im];
                tmp_psd_fbf = sp_re * sp_re + sp_im * sp_im + AWI_EPS;
                fbf_auto_psd_beam[band_id] = p->cfg.alpha_psd * fbf_auto_psd_beam[band_id] +
                                            alpha_psd_1 * tmp_psd_fbf;

                // compute psd for aec output signal 
                sp_re2 = aec_eo_sp_beam[idx_re];
                sp_im2 = aec_eo_sp_beam[idx_im];
                tmp_psd_eo = sp_re2 * sp_re2 + sp_im2 * sp_im2 + AWI_EPS;
                aec_eo_auto_psd_beam[band_id] = p->cfg.alpha_psd * aec_eo_auto_psd_beam[band_id] +
                                                alpha_psd_1 * tmp_psd_eo;

                // compute real part of cross psd between fbf and aec output 
                tmp_cpsd = sp_re * sp_re2 + sp_im * sp_im2;
                tmp_cpsd = tmp_cpsd > 0 ? tmp_cpsd : -tmp_cpsd;
                fbf_aec_eo_cross_psd_beam[band_id] = p->cfg.alpha_psd * fbf_aec_eo_cross_psd_beam[band_id] +
                                                    alpha_psd_1 * tmp_cpsd;

                // compute psd for estimated linear echo 
                sp_re = aec_yo_sp_beam[idx_re];
                sp_im = aec_yo_sp_beam[idx_im];
                tmp_psd_yo = sp_re * sp_re + sp_im * sp_im + AWI_EPS;
                aec_yo_auto_psd_beam[band_id] = p->cfg.alpha_psd * aec_yo_auto_psd_beam[band_id] +
                                                alpha_psd_1 * tmp_psd_yo;
                // compute subband erle
                sub_erle = AWI_MAX(fbf_auto_psd_beam[band_id] / (aec_eo_auto_psd_beam[band_id] + AWI_EPS), 0.1);    

                // compute subband erl
                sub_erl = aec->aec_Far_autoPsd[band_id] / (fbf_auto_psd_beam[band_id] + AWI_EPS);      
                
                sub_erl = AWI_MAX(sub_erl, aec->fullband_spk_mic_erl_chs[beam_id]);
                
                sub_erl_inv = 1.f / sub_erl;
                // estimate residue echo psd 
                // freq_spk_psd = p->recur_spk_auto_psd[band_id];
                freq_spk_psd = aec->max_spk_freq_psd[band_id];
                initial_est_echo_psd = sub_erl_inv * freq_spk_psd / sub_erle;  
                initial_est_echo_psd = p->cfg.alpha_sub_erl * initial_est_echo_psd + alpha_global_erl * fullband_erl_inv * global_erle_weight_inv * freq_spk_psd;
                initial_est_echo_psd = AWI_MAX(initial_est_echo_psd, aec_yo_auto_psd_beam[band_id]);
               
                // estimate signal-to-residue_echo ratio 
                inst_post_ser = tmp_psd_gsc / initial_est_echo_psd;
                alpha_sub_erle = p->cfg.alpha_sub_erle_low;

                // apply different smoothing factor based on ser
                if (inst_post_ser < 3) 
                {
                    alpha_sub_erle = p->cfg.alpha_sub_erle_high; 
                }
                alpha_full_erle = 1 - alpha_sub_erle;
                combined_erle = alpha_sub_erle * sub_erle + alpha_full_erle * global_erle_weight;
                
                // compute initial aes gain based on wiener filter method
                aes_gain_beam[band_id] = (fbf_aec_eo_cross_psd_beam[band_id] + AWI_EPS) /
                                        (fbf_auto_psd_beam[band_id] + 
                                        combined_erle * initial_est_echo_psd + AWI_EPS);

                aes_gain_beam[band_id] = tmp_gain = AWI_MIN(aes_gain_beam[band_id], 1);
                
                p->min_aes_gain[band_id] = AWI_MIN(p->min_aes_gain[band_id], tmp_gain);

                // compute psd of estimated near end signal 
                tmp_gain_minus = 1 - tmp_gain;
                tmp_residue_echo_psd = tmp_gain_minus * tmp_gain_minus * tmp_psd_gsc;
                p->psy_near_auto_psd[band_id] = tmp_gain * tmp_gain * tmp_psd_gsc;
                
                // smooth psd of residue echo 
                if( tmp_residue_echo_psd > residue_echo_psd_beam[band_id] )
                {
                    residue_echo_psd_beam[band_id] = p->cfg.alpha_echo_at * residue_echo_psd_beam[band_id] +
                                                    alpha_echo_at_1 * tmp_residue_echo_psd;
                }
                else
                {
                    residue_echo_psd_beam[band_id] = p->cfg.alpha_echo_ar * residue_echo_psd_beam[band_id] +
                                                    alpha_echo_ar_1 * tmp_residue_echo_psd;
                }
            }
            else
            {
                p->psy_near_auto_psd[band_id] = tmp_psd_gsc;

                residue_echo_psd_beam[band_id] = 0;

                aes_gain_beam[band_id] = 1;
            }

            // compute low and high band energy 
            if (band_id > 2 && band_id < p->cfg.mid_band_id_ey)
            {
                inst_low_band_energy += p->psy_near_auto_psd[band_id];
                inst_gsc_low_energy += tmp_psd_gsc;
            }
            else if ( band_id < p->cfg.high_band_id_th )
            {
                inst_high_band_energy += p->psy_near_auto_psd[band_id];
                inst_gsc_high_energy += tmp_psd_gsc;
            }

        }

        // psychological postfilter
        apply_psycho_acoustic_postfilter(p->psy_near_auto_psd, residue_echo_psd_beam, aes_gain_beam, p->cfg.psy_relax_factor);
        
        // smooth low and high band energy, accumulate low and high band energy to get fullband energy, and compute high-to-low energy ratio for estimated near end signal
        p->recur_low_band_energy[beam_id] = p->cfg.alpha_energy * p->recur_low_band_energy[beam_id] +
                                            alpha_energy_1 * inst_low_band_energy;
        p->recur_high_band_energy[beam_id] = p->cfg.alpha_energy * p->recur_high_band_energy[beam_id] +
                                             alpha_energy_1 * inst_high_band_energy;
        near_energy = p->recur_low_band_energy[beam_id] + p->recur_high_band_energy[beam_id];

        near_hl_energy_ratio = p->recur_high_band_energy[beam_id] / (p->recur_low_band_energy[beam_id] + AWI_EPS);

        // smooth low and high band energy, accumulate low and high band energy to get fullband energy, and compute high-to-low energy ratio for gsc output signal 
        p->recur_gsc_low_band_energy[beam_id] = p->cfg.alpha_energy * p->recur_gsc_low_band_energy[beam_id] +
                                                alpha_energy_1 * inst_gsc_low_energy;
        p->recur_gsc_high_band_energy[beam_id] = p->cfg.alpha_energy * p->recur_gsc_high_band_energy[beam_id] +
                                                alpha_energy_1 * inst_gsc_high_energy;
        gsc_energy = p->recur_gsc_low_band_energy[beam_id] + p->recur_gsc_high_band_energy[beam_id];

        gsc_hl_energy_ratio = p->recur_gsc_high_band_energy[beam_id] / (p->recur_gsc_low_band_energy[beam_id] + AWI_EPS);

        // compute change of high-to-low energy ratio after aes 
        hl_energy_ratio_change = near_hl_energy_ratio / (gsc_hl_energy_ratio + AWI_EPS);

        // compute fullband aes gain 
        fullband_aes_gain = ( near_energy + AWI_EPS ) / ( gsc_energy + AWI_EPS );
        
        // special processing to eliminate initial head echo 
        if ( p->head_counter < p->cfg.head_length && max_spk_energy > p->cfg.spk_head_energy_th)
        {
            p->head_counter++;
            float echo_energy = 0;
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                if ( residue_echo_psd_beam[band_id] > p->cfg.spk_head_freq_psd_th )
                {
                    aes_gain_beam[band_id] *= AWI_MIN(p->cfg.spk_head_freq_psd_th / (residue_echo_psd_beam[band_id] + AWI_EPS), 1);
                }
                echo_energy += residue_echo_psd_beam[band_id];
            }
            fullband_aes_gain = AWI_MAX(( near_energy - echo_energy + AWI_EPS ) / ( gsc_energy + AWI_EPS ), AWI_EPS);
        }

        p->cfg.alpha_fullband = Sigmoid_mod(fullband_aes_gain, p->cfg.mid_fullband_gain, p->cfg.slope_fullband_gain);
        alpha_fullband_1 = 1 - p->cfg.alpha_fullband;

        p->min_fullband_aes_gain = AWI_MIN(p->min_fullband_aes_gain, fullband_aes_gain);

        // apply fullband and subband aes gain 
        if (max_spk_energy > spk_energy_th_high &&
        ( hl_energy_ratio_change > p->cfg.hl_energy_ratio_change
        || near_energy < p->cfg.full_band_energy_th 
        || ( fullband_aes_gain < p->cfg.weak_aes_gain_th ) ) )
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;
                tmp_gain = aes_gain_beam[band_id] = p->cfg.sub_aes_gain_floor;
                p->min_aes_gain[band_id] = AWI_MIN(p->min_aes_gain[band_id], tmp_gain);
                aes_sp_beam[idx_re] = tmp_gain * gsc_sp_beam[idx_re];
                aes_sp_beam[idx_im] = tmp_gain * gsc_sp_beam[idx_im];
            }
        }
        else
        {
            for (int band_id = 0; band_id < p->cfg.high_band_id_th; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;
                tmp_gain = aes_gain_beam[band_id];
                tmp_gain = p->cfg.alpha_fullband * fullband_aes_gain + alpha_fullband_1 * tmp_gain;
                aes_gain_beam[band_id] = tmp_gain = AWI_MAX(tmp_gain, p->cfg.sub_aes_gain_floor);
                p->min_aes_gain[band_id] = AWI_MIN(p->min_aes_gain[band_id], tmp_gain);
                aes_sp_beam[idx_re] = tmp_gain * gsc_sp_beam[idx_re];
                aes_sp_beam[idx_im] = tmp_gain * gsc_sp_beam[idx_im];
            }

            for (int band_id = p->cfg.high_band_id_th; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;
                aes_gain_beam[band_id] = tmp_gain = p->cfg.high_band_gain_cnt;
                p->min_aes_gain[band_id] = AWI_MIN(p->min_aes_gain[band_id], tmp_gain);
                aes_sp_beam[idx_re] = tmp_gain * gsc_sp_beam[idx_re];
                aes_sp_beam[idx_im] = tmp_gain * gsc_sp_beam[idx_im];
            }
        }

        fbf_sp_beam += sp_offset;
        aec_eo_sp_beam += sp_offset;
        aec_yo_sp_beam += sp_offset;
        fbf_auto_psd_beam += psd_offset;
        aec_yo_auto_psd_beam += psd_offset;
        fbf_aec_eo_cross_psd_beam += psd_offset;
        gsc_sp_beam += sp_offset;
        residue_echo_psd_beam += psd_offset;
        aes_sp_beam += sp_offset;
    }
}



void apply_psycho_acoustic_postfilter(float *desired_sig_psd_est, float *residue_echo_psd_est,
                                      float *aes_gain_beam, float psy_relax_factor)
{
    int lower_id, upper_id, diff_band_id;
    float bark_band_energy[AWI_FREQ_BANDS_VALID], spread_bark_band_energy[AWI_FREQ_BANDS_VALID];
    float algebraic_avg, geometric_avg, geometric_avg_low, geometric_avg_high, alg2geo_ratio;
    const float max_alg2geo_ratio = -60;
    float alpha, offset_dB, spread_gain;
    float max_spread_bark_band_energy = -1000;
    float mask_th[AWI_FREQ_BANDS_VALID];

    /* step1: compute energy in Bark scale */
    for(int band = 0; band < AWI_FREQ_BANDS_VALID; band++)
    {
        bark_band_energy[band] = 0;
        if ( mask_band_freq_st[band] < 8000 )
        {
            lower_id = AWI_FREQ_LOW[band] - 1;
            upper_id = AWI_MIN(AWI_FREQ_HIGH[band]-1, AWI_FRAME_BAND_COUNT);
            for(int k = lower_id; k < upper_id; k++)
                bark_band_energy[band] += desired_sig_psd_est[k];
        }
    }
    /* step2: compute masking between different critical bark bands */
    for(int band_i = 0; band_i < AWI_FREQ_BANDS_VALID; band_i++)
    {
        spread_bark_band_energy[band_i] = 0;
        for(int  band_j = 0; band_j < AWI_FREQ_BANDS_VALID; band_j++)
        {
            diff_band_id = AWI_MAX(band_i - band_j, band_j - band_i);
            spread_bark_band_energy[band_i] += bark_band_energy[band_j] * SpreadFunc[diff_band_id];
        }
    }

    /* step3: compute masking threshold */
    /* compute algebraic and geometric average of spreading bark energy */
    algebraic_avg = 0, geometric_avg_low = 1, geometric_avg_high = 1;
    for(int band_id = 0; band_id < AWI_FREQ_BANDS_VALID; band_id++)
    {
        algebraic_avg += spread_bark_band_energy[band_id];
        max_spread_bark_band_energy = AWI_MAX(max_spread_bark_band_energy, spread_bark_band_energy[band_id]);
    }

    for(int band_id = 0; band_id < 17; band_id++)
    {
        geometric_avg_low *= (spread_bark_band_energy[band_id] + 1e-20f) / (max_spread_bark_band_energy + 1e-20f) * 10;
    }

    for(int band_id = 17; band_id < AWI_FREQ_BANDS_VALID; band_id++)
    {
        geometric_avg_high *= (spread_bark_band_energy[band_id] + 1e-20f) / (max_spread_bark_band_energy + 1e-20f) * 1000;
    }
    geometric_avg_high *= 10;

    algebraic_avg /= AWI_FREQ_BANDS_VALID;

    geometric_avg_low = powf(geometric_avg_low, 1.f / ( AWI_FREQ_BANDS - 2.f ) );
    geometric_avg_high = powf(geometric_avg_high, 1.f / ( AWI_FREQ_BANDS - 2.f ) );
    geometric_avg = geometric_avg_low * geometric_avg_high * max_spread_bark_band_energy / 10.f / sqrtf(10);

    alg2geo_ratio = 10 * log10f((geometric_avg + 1e-20f) / (algebraic_avg + 1e-20f));
    alg2geo_ratio = AWI_MIN(alg2geo_ratio, 0);
    alpha = AWI_MIN(alg2geo_ratio / max_alg2geo_ratio, 1);
    for(int band_id = 0; band_id < AWI_FREQ_BANDS_VALID; band_id++)
    {
        offset_dB = alpha * ( 14.5f + band_id + 1 ) + ( 1 - alpha ) * 5.5f;
        mask_th[band_id] = spread_bark_band_energy[band_id] / powf(2, 0.33219281f * offset_dB);
        /* step4: renormalization */
        spread_gain = (spread_bark_band_energy[band_id] + 1e-20f) / (bark_band_energy[band_id] + 1e-20f);
        mask_th[band_id] /= spread_gain;
    }

    /* step5: set frequency bins within the same critical bark band with the same masking threshold */
    /* step6: re-compute gain considering masking threshold and estimated residue echo PSD */
    for(int band_id = 0; band_id < AWI_FREQ_BANDS_VALID; band_id++)
    {
        if ( mask_band_freq_st[band_id] < 8000 )
        {
            lower_id = AWI_FREQ_LOW[band_id] - 1;
            upper_id = AWI_MIN(AWI_FREQ_HIGH[band_id]-1, AWI_FRAME_BAND_COUNT);
            for(int k = lower_id; k < upper_id; k++)
            {
                if ( desired_sig_psd_est[band_id] > 1e-9f && k < 80 )
                {
                    aes_gain_beam[k] = psy_relax_factor * sqrtf((mask_th[band_id] + 1e-20f)
                            / (residue_echo_psd_est[k] + 1e-20f));
                }
                aes_gain_beam[k] = AWI_MIN(aes_gain_beam[k], 1);
            }
        }
    }

}

