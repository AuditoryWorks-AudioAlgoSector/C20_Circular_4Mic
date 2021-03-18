#include "awi_ns.h"
#include "awi_utils.h"

void awi_ns_init(awi_ns_t *p)
{
    for(int i = 0; i < AWI_FRAME_BAND_COUNT; i++)
    {
        p->gain_speech[i]        = 0.1f;
        p->inst_post_snr_prev[i] = 0.001f;
        p->recur_priori_snr[i]   = 0.001f;
    }

    memset(p->local_recur_priori_snr, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->global_recur_priori_snr, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->priori_snr, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->inst_post_snr_cng, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->local_speech_freq_prob, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->global_speech_freq_prob, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->ns_gain_seq, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->low_avg_gain, 0, sizeof(float) * AWI_FREQ_BANDS);
    memset(p->high_avg_gain, 0, sizeof(float) * AWI_FREQ_BANDS);
    memset(p->full_avg_gain, 0, sizeof(float) * AWI_FREQ_BANDS);
    memset(p->high_gain_stat_counter, 0, sizeof(float) * AWI_FREQ_BANDS);
}


void awi_ns_process(awi_ns_t *p, float *input_sp, float *out_sp, float *input_psd, float *recur_block_psd, float* bkg_est_psd, int speech_flag)
{
    float instant_noisy_psd = 0;
    float ns_gain = 0;
    float speech_absence_freq = 0;
    float vk = 0;
    float priori_speech_presence = 0;
    float avg_recur_priori_snr = 0;
    float tmp_snr = 0;
    float speech_frame_prob = 0;
    float speech_absence_prob = 0;
    int lower = 0, upper = 0, vk_id = 0;
    float gain_speech_band = 0, priori_snr_band = 0;
    float norm_recur_priori_snr = 0;
    int idx_re, idx_im;
    float high_band_avg_gain;

    float cnt1 = VK_TABLE1_START - 0.5f * p->cfg.vk_table_step1;
    float cnt2 = VK_TABLE2_START - 0.5f * p->cfg.vk_table_step2;
    float cnt3 = VK_TABLE3_START - 0.5f * p->cfg.vk_table_step3;
    float cnt4 = VK_TABLE4_START - 0.5f * p->cfg.vk_table_step4;
    float step1_inv = 1.f / p->cfg.vk_table_step1;
    float step2_inv = 1.f / p->cfg.vk_table_step2;
    float step3_inv = 1.f / p->cfg.vk_table_step3;
    float step4_inv = 1.f / p->cfg.vk_table_step4;
    int len1_minus1 = VK_TABLE1_LEN - 1;
    int len2_minus1 = VK_TABLE2_LEN - 1;
    int len3_minus1 = VK_TABLE3_LEN - 1;
    int len4_minus1 = VK_TABLE4_LEN - 1;

    float log_cnt1 = 1.f / log10f( p->cfg.norm_recur_priori_snr_max / p->cfg.norm_recur_priori_snr_min );
    float multi_cnt1 = p->cfg.norm_recur_priori_snr_min * p->cfg.avg_recur_prior_snr_peak;
    float multi_cnt2 = p->cfg.norm_recur_priori_snr_max * p->cfg.avg_recur_prior_snr_peak;
    float multi_cnt3 = 1.f / p->cfg.avg_recur_prior_snr_peak / p->cfg.norm_recur_priori_snr_min;

    float alpha_priori_snr_1, alpha_recur_priori_snr_1;
    alpha_priori_snr_1 = 1 - p->cfg.alpha_priori_snr;
    alpha_recur_priori_snr_1 = 1 - p->cfg.alpha_recur_priori_snr;

    float alpha_sub_gain_1 = 1 - p->cfg.alpha_sub_gain;

    if ( !speech_flag )
    {
        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            gain_speech_band = p->gain_speech[band];
            instant_noisy_psd = input_psd[band];
            p->inst_post_snr_cng[band] = 1;
            tmp_snr = AWI_EPS;
            p->priori_snr[band] = p->cfg.alpha_priori_snr * (gain_speech_band * gain_speech_band) * \
                    p->inst_post_snr_prev[band] + alpha_priori_snr_1 * tmp_snr;

            p->recur_priori_snr[band] = p->cfg.alpha_recur_priori_snr * p->recur_priori_snr[band] + \
                    alpha_recur_priori_snr_1 *  p->priori_snr[band];

            avg_recur_priori_snr += p->recur_priori_snr[band];
        }
    }
    else
    {
        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            gain_speech_band = p->gain_speech[band];
            instant_noisy_psd = input_psd[band];
            tmp_snr = (instant_noisy_psd + AWI_EPS) / (bkg_est_psd[band] + AWI_EPS);
            p->inst_post_snr_cng[band] = tmp_snr;
            tmp_snr = AWI_MAX(tmp_snr-1, AWI_EPS);
            p->priori_snr[band] = p->cfg.alpha_priori_snr * (gain_speech_band * gain_speech_band) * \
                    p->inst_post_snr_prev[band] + alpha_priori_snr_1 * tmp_snr;

            p->recur_priori_snr[band] = p->cfg.alpha_recur_priori_snr * p->recur_priori_snr[band] + \
                    alpha_recur_priori_snr_1 *  p->priori_snr[band];

            p->recur_priori_snr[band] = AWI_MAX(p->recur_priori_snr[band], p->cfg.snr_min);
            p->recur_priori_snr[band] = AWI_MIN(p->recur_priori_snr[band], p->cfg.snr_max);

            avg_recur_priori_snr += p->recur_priori_snr[band];

        }
    }

    avg_recur_priori_snr /= AWI_FRAME_BAND_COUNT;
    
    if( avg_recur_priori_snr < p->cfg.avg_recur_priori_snr_min )
    {
        speech_frame_prob = 0;
    }
    else if ( avg_recur_priori_snr > p->cfg.avg_recur_priori_snr_prev )
    {
        speech_frame_prob = 1;
        p->cfg.avg_recur_prior_snr_peak = AWI_MIN(AWI_MAX(avg_recur_priori_snr,p->cfg.avg_recur_priori_snr_min),
                p->cfg.avg_recur_priori_snr_max);
    }
    else
    {
        if ( avg_recur_priori_snr <= multi_cnt1 )
        {
            speech_frame_prob = 0;
        }
        else if ( avg_recur_priori_snr >= multi_cnt2 )
        {
            speech_frame_prob = 1;
        }
        else
        {
            speech_frame_prob = log10f(avg_recur_priori_snr * multi_cnt3) * log_cnt1;
        }
    }

    p->cfg.avg_recur_priori_snr_prev = avg_recur_priori_snr;


    awi_smoothing_1D_vector(p->local_recur_priori_snr, p->recur_priori_snr, LOCAL_SMOOTHING_WINDOW, AWI_FRAME_BAND_COUNT, NS_LOCAL_WIN_LEN);
    awi_smoothing_1D_vector(p->global_recur_priori_snr, p->recur_priori_snr, GLOBAL_SMOOTHING_WINDOW, AWI_FRAME_BAND_COUNT, NS_GLOBAL_WIN_LEN);
    

    for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {
        norm_recur_priori_snr = p->local_recur_priori_snr[band] / p->cfg.snr_max;
        if ( norm_recur_priori_snr <= p->cfg.norm_recur_priori_snr_min )
            p->local_speech_freq_prob[band] = 0;
        else if ( norm_recur_priori_snr >= p->cfg.norm_recur_priori_snr_max )
            p->local_speech_freq_prob[band] = 1;
        else
        {
            p->local_speech_freq_prob[band] = log10f( norm_recur_priori_snr / p->cfg.norm_recur_priori_snr_min ) \
                        * log_cnt1;
        }

        norm_recur_priori_snr = p->global_recur_priori_snr[band] / p->cfg.snr_max;
        if ( norm_recur_priori_snr <= p->cfg.norm_recur_priori_snr_min )
            p->global_speech_freq_prob[band] = 0;
        else if ( norm_recur_priori_snr >= p->cfg.norm_recur_priori_snr_max )
            p->global_speech_freq_prob[band] = 1;
        else
        {
            p->global_speech_freq_prob[band] = log10f( norm_recur_priori_snr / p->cfg.norm_recur_priori_snr_min ) \
                        * log_cnt1;
        }
    }


    for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {

        priori_snr_band = p->priori_snr[band];

        speech_absence_prob = 1 - p->local_speech_freq_prob[band] * speech_frame_prob * p->global_speech_freq_prob[band];

        speech_absence_freq = AWI_MIN(speech_absence_prob, p->cfg.max_speech_absence_prob);

        vk =( priori_snr_band * p->inst_post_snr_cng[band]) / (1 + priori_snr_band);

        priori_speech_presence = 1 / ( 1 + speech_absence_freq / ( 1 - speech_absence_freq) * ( 1 + priori_snr_band )
                * expf(-vk));

        if(vk < VK_TABLE1_START)
            vk = VK_TABLE1_START;

        if(vk > VK_TABLE4_END)
            vk = VK_TABLE4_END;

        if ( vk >= VK_TABLE1_START && vk < VK_TABLE1_END )
        {
            vk_id = (int) ( ( vk - cnt1 ) * step1_inv );
            vk_id = AWI_MIN(vk_id, len1_minus1);
            gain_speech_band = priori_snr_band / ( 1 + priori_snr_band ) * EXP_INT_TABLE1[vk_id];
        }
        else if( vk >= VK_TABLE2_START && vk < VK_TABLE2_END)
        {
            vk_id = (int) ( ( vk - cnt2 ) * step2_inv );
            vk_id = AWI_MIN(vk_id, len2_minus1);
            gain_speech_band = priori_snr_band / ( 1 + priori_snr_band ) * EXP_INT_TABLE2[vk_id];
        }
        else if ( vk >= VK_TABLE3_START && vk < VK_TABLE3_END)
        {
            vk_id = (int) ( ( vk - cnt3 ) * step3_inv );
            vk_id = AWI_MIN(vk_id, len3_minus1);
            gain_speech_band = priori_snr_band / ( 1 + priori_snr_band ) * EXP_INT_TABLE3[vk_id];
        }
        else
        {
            vk_id = (int) ( ( vk - cnt4 ) * step4_inv );
            vk_id = AWI_MIN(vk_id, len4_minus1);
            gain_speech_band = priori_snr_band / ( 1 + priori_snr_band ) * EXP_INT_TABLE4[vk_id];
        }

        gain_speech_band = AWI_MIN(AWI_MAX(gain_speech_band, p->cfg.gain_min), 1);

        p->gain_speech[band] = gain_speech_band;

        p->ns_gain_seq[band] = powf(gain_speech_band, priori_speech_presence) *
                powf(p->cfg.gain_min, ( 1- priori_speech_presence));

    }


    for(int j = 0; j < AWI_FREQ_BANDS_VALID; j++)
    {
        p->high_gain_stat_counter[j] = 0;
        p->high_avg_gain[j] = 0;
        p->low_avg_gain[j] = 0;
        if ( mask_band_freq_st[j] < 8000 )
        {
            lower = AWI_FREQ_LOW[j] - 1;
            upper = AWI_MIN(AWI_FREQ_HIGH[j], AWI_FRAME_BAND_COUNT);
            for(int k = lower; k < upper; k++)
            {
                if ( p->ns_gain_seq[k] > p->cfg.high_gain_th )
                {
                    p->high_gain_stat_counter[j]++;
                    p->high_avg_gain[j] += p->ns_gain_seq[k];
                }
                else
                    p->low_avg_gain[j] += p->ns_gain_seq[k];
            }
            float total_counter = upper - lower;
            if ( p->high_gain_stat_counter[j] == total_counter )
            {
                p->high_avg_gain[j] /= total_counter;
                p->low_avg_gain[j] = 0;
            }
            else if ( p->high_gain_stat_counter[j] == 0 )
            {
                p->low_avg_gain[j] /= total_counter;
                p->high_avg_gain[j] = 0;
            } else
            {
                p->high_avg_gain[j] /= p->high_gain_stat_counter[j];
                p->low_avg_gain[j] /= ( total_counter - p->high_gain_stat_counter[j]);
            }

            if ( p->high_gain_stat_counter[j] > p->cfg.high_counter_ratio * total_counter )
                p->full_avg_gain[j] = p->high_avg_gain[j];
            else if ( p->high_gain_stat_counter[j] < p->cfg.low_counter_ratio * total_counter )
                p->full_avg_gain[j] = p->low_avg_gain[j];
            else
            {
                float ratio = p->high_gain_stat_counter[j] / total_counter;
                p->full_avg_gain[j] = ratio * p->high_avg_gain[j] + ( 1 - ratio ) * p->low_avg_gain[j];
            }

            if ( j <= 6 )
            {
                for(int k = lower; k < upper; k++)
                {
                    if ( recur_block_psd[k] < 5e-8f ) p->ns_gain_seq[k] = p->full_avg_gain[j];
                }
            }
            else
            {
                for(int k = lower; k < upper; k++)
                {
                    p->ns_gain_seq[k] = p->full_avg_gain[j];
                }
            }   
        }
    }

    high_band_avg_gain = 0;
    for(int band = p->cfg.low_freq_id; band < p->cfg.high_freq_id; band++) 
    {
        high_band_avg_gain += p->ns_gain_seq[band];
    }
    high_band_avg_gain /= (p->cfg.high_freq_id - p->cfg.low_freq_id);

    if (high_band_avg_gain < p->cfg.high_band_avg_gain_th)
    {
        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            p->ns_gain_seq[band] = p->cfg.alpha_sub_gain * p->ns_gain_seq[band] + alpha_sub_gain_1 * high_band_avg_gain;
        }
    }


    for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {
        idx_re = band + band;
        idx_im = idx_re + 1;
        ns_gain = AWI_MAX(p->ns_gain_seq[band], p->cfg.gain_min);
        out_sp[idx_re] = EQ_COMPEN[band] * ns_gain * input_sp[idx_re];
        out_sp[idx_im] = EQ_COMPEN[band] * ns_gain * input_sp[idx_im];
    }

    memcpy(p->inst_post_snr_prev, p->inst_post_snr_cng, sizeof(float) * AWI_FRAME_BAND_COUNT);
}

