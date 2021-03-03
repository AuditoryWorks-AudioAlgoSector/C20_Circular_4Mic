#include <string.h>
#include <math.h>
#include "awi_constant.h"
#include "awi_nse.h"
#include "awi_utils.h"

void awi_nse_init(awi_nse_t *p)
{
    p->init_counter   = 0;
    p->init_complete  = 0;
    memset(p->avg_inst_noisy_psd, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->bkg_est_psd, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->recur_block_psd, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->inv_Qeq, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->Qeq_inv_global, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->Qeq_inv_local, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->global_bias_origin, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->local_bias_origin, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->local_ms_ns_est, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);
    memset(p->global_ms_ns_est, 0, sizeof(float) * AWI_FRAME_BAND_COUNT);

    for(int i = 0; i < AWI_FRAME_BAND_COUNT * p->cfg.V; i++)
        p->local_ms_win[i] = 1000.f;

    for(int i = 0; i < AWI_FRAME_BAND_COUNT * p->cfg.U; i++)
        p->global_ms_win[i] = 1000.f;
}


void awi_nse_process_multichannel(awi_nse_t *p, float *input_sp, int auxi_vad_flag, int enable_vad_auxi)
{
    float *input_sp_beam;
    float sp_re, sp_im;
    int idx_re, idx_im;
    float *avg_inst_aec_eo_psd = p->avg_inst_noisy_psd;
    const int offset_sp = AWI_FRAME_BAND_COUNT * 2;

    input_sp_beam = input_sp;
    for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {
        idx_re = band + band;
        idx_im = idx_re + 1;
        sp_re = input_sp_beam[idx_re];
        sp_im = input_sp_beam[idx_im];
        avg_inst_aec_eo_psd[band] = sp_re * sp_re + sp_im * sp_im;
    }

    input_sp_beam = input_sp + offset_sp;
    for(int channel = 1; channel < AWI_MIC_CHANNEL; channel++)
    {
        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            idx_re = band + band;
            idx_im = idx_re + 1;
            sp_re = input_sp_beam[idx_re];
            sp_im = input_sp_beam[idx_im];
            avg_inst_aec_eo_psd[band] = AWI_MAX(sp_re * sp_re + sp_im * sp_im, avg_inst_aec_eo_psd[band]);
        }
        input_sp_beam += offset_sp;
    }

    awi_nse_process_local_func(p, auxi_vad_flag, enable_vad_auxi);
}

void awi_nse_process_mono(awi_nse_t *p, float *input_sp, int auxi_vad_flag, int enable_vad_auxi)
{
    float re, im;
    int idx_re, idx_im;
    float *avg_inst_aec_eo_psd = p->avg_inst_noisy_psd;

    for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {
        idx_re = band + band;
        idx_im = idx_re + 1;
        re = input_sp[idx_re];
        im = input_sp[idx_im];
        avg_inst_aec_eo_psd[band] = re * re + im * im;
    }

    awi_nse_process_local_func(p, auxi_vad_flag, enable_vad_auxi);
}

void awi_nse_process_local_func(awi_nse_t *p, int auxi_vad_flag, int enable_vad_auxi)
{
    float tmp_recur, power;

    float *avg_inst_aec_eo_psd = p->avg_inst_noisy_psd;
    float *bkg_est_psd = p->bkg_est_psd;
    float *recur_block_psd = p->recur_block_psd;
    float *first_moment_recur_psd = p->first_moment_recur_psd;
    float *second_moment_recur_psd = p->second_moment_recur_psd;
    float *inv_Qeq = p->inv_Qeq;
    float *Qeq_inv_global = p->Qeq_inv_global ;
    float *global_bias_origin = p->global_bias_origin;
    float *local_bias_origin = p->local_bias_origin;
    float *Qeq_inv_local = p->Qeq_inv_local;
    float *local_ms_win = p->local_ms_win;
    float *global_ms_win = p->global_ms_win;
    float *local_ms_ns_est = p->local_ms_ns_est;
    float *global_ms_ns_est = p->global_ms_ns_est;

    float mvby2 = 2 * p->cfg.MV;
    float mdby2 = 2 * p->cfg.MD;
    float md_minus1_inv = 1.f / (1.f - p->cfg.MD);
    float mv_minus1_inv = 1.f / (1.f - p->cfg.MV);
    float dminus1by2 = (p->cfg.D - 1) * 2;
    float vminus1by2 = (p->cfg.V - 1) * 2;

    float alpha_bkg_inc_1 = 1 - p->cfg.alpha_bkg_inc;
    float alpha_bkg_dec_1 = 1 - p->cfg.alpha_bkg_dec;
    float alpha_init_speech_1 = 1 - p->cfg.alpha_init_speech;
    float alpha_init_bkg = 1 - p->cfg.alpha_init_bkg;

    if (!p->init_complete)
    {
        p->init_counter++;
        if(p->init_counter <= 2)
        {
            for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
            {
                power = avg_inst_aec_eo_psd[band];
                recur_block_psd[band] = power;
                bkg_est_psd[band] = power;
            }
            return;
        }
        else if (p->init_counter < p->cfg.nis_blocks)
        {
            for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
            {
                power = avg_inst_aec_eo_psd[band];
                recur_block_psd[band] = p->cfg.alpha_init_speech * recur_block_psd[band] + alpha_init_speech_1 * power;
                bkg_est_psd[band] = p->cfg.alpha_init_bkg * bkg_est_psd[band] + alpha_init_bkg * power;
            }

            return;
        }

        if(p->init_counter == p->cfg.nis_blocks)
        {
            for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
            {

                power = avg_inst_aec_eo_psd[band];

                recur_block_psd[band] = tmp_recur = p->cfg.alpha_init_speech * recur_block_psd[band] + alpha_init_speech_1 * power;
                bkg_est_psd[band] = p->cfg.alpha_init_bkg * bkg_est_psd[band] + alpha_init_bkg * power;
                first_moment_recur_psd[band] = tmp_recur;
                second_moment_recur_psd[band] = tmp_recur * tmp_recur;
            }
            p->init_complete = 1;
        }
    }

    if(p->init_complete)
    {

        float summed_noisy_psd = 0;
        float summed_inst_psd = 0;
        float summed_bkg_psd = 0;
        float overall_SNR, alpha_opt_min;
        float psd_cmp_ratio_1, inst_alpha_crt, alpha_multi;
        float alpha_opt_, beta_, beta_1, recur_;
        float first_moment, second_moment;
        float bkg_psd;
        float iv_eq, eq;
        float avg_iv_eqs = 0;
        float inv_qeq_g, inv_qeq_l;
        float min_psd;
        float noise_slope_max;
        float *local_ms_win_jth;
        float *global_ms_win_jth;
        int offset = 0;

        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            summed_noisy_psd += recur_block_psd[band];
            summed_inst_psd += avg_inst_aec_eo_psd[band];
            summed_bkg_psd += bkg_est_psd[band];
        }

        overall_SNR = summed_inst_psd / (summed_bkg_psd + p->cfg.epsth);
        overall_SNR = AWI_MAX(overall_SNR - 1, 0.001f);
        alpha_opt_min = powf(overall_SNR, p->cfg.snr_power);
        alpha_opt_min = AWI_MIN(alpha_opt_min, p->cfg.alpha_opt_min_max);

        psd_cmp_ratio_1 = ( summed_noisy_psd / ( summed_inst_psd + p->cfg.epsth ) - 1 );
        inst_alpha_crt = AWI_MAX(1 / ( 1 +  psd_cmp_ratio_1 * psd_cmp_ratio_1), 0.7f);

        p->cfg.alpha_correction_factor = 0.7f * p->cfg.alpha_correction_factor + 0.3f * inst_alpha_crt;
        alpha_multi = p->cfg.alpha_correction_factor * p->cfg.alpha_opt_max;

        for(int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            psd_cmp_ratio_1 = ( recur_block_psd[band] / ( bkg_est_psd[band] + p->cfg.epsth ) - 1 );
            alpha_opt_ = alpha_multi / ( 1 + psd_cmp_ratio_1 * psd_cmp_ratio_1 );
            alpha_opt_ = AWI_MIN(AWI_MAX(alpha_opt_, alpha_opt_min), p->cfg.alpha_opt_max);
            recur_block_psd[band] = recur_ = alpha_opt_ * recur_block_psd[band] +
                    (1 - alpha_opt_) * avg_inst_aec_eo_psd[band];
            beta_ = alpha_opt_ * alpha_opt_;
            beta_ = AWI_MIN(beta_, p->cfg.beta_max);
            beta_1 = 1 - beta_;
            first_moment_recur_psd[band]  = first_moment = beta_ * first_moment_recur_psd[band] +
                    beta_1 * recur_;
            second_moment_recur_psd[band] = second_moment = beta_ * second_moment_recur_psd[band] +
                    beta_1 * recur_ * recur_;

            bkg_psd = p->bkg_est_psd[band];
            iv_eq = (second_moment - first_moment * first_moment) /
                    (2 * bkg_psd  * bkg_psd + p->cfg.epsth);
            iv_eq = AWI_MAX(iv_eq, -iv_eq);
            iv_eq = AWI_MIN(iv_eq, 0.5f);
            inv_Qeq[band] = iv_eq;
            avg_iv_eqs += iv_eq;
            eq = 1.f / iv_eq;
            Qeq_inv_global[band] = inv_qeq_g  = (eq - mdby2) * md_minus1_inv;
            global_bias_origin[band] = 1 + dminus1by2  / inv_qeq_g;
            Qeq_inv_local[band] = inv_qeq_l = (eq - mvby2) * mv_minus1_inv;
            local_bias_origin[band] = 1 + vminus1by2  /  inv_qeq_l;
        }

        avg_iv_eqs = avg_iv_eqs / AWI_FRAME_BAND_COUNT;

        offset = p->cfg.local_ms_counter * AWI_FRAME_BAND_COUNT;
        for(int i = 0; i < AWI_FRAME_BAND_COUNT; i++)
        {
            local_ms_win[i + offset] = recur_block_psd[i];
            min_psd = 10000.f;
            local_ms_win_jth = local_ms_win;
            for(int j = 0; j < p->cfg.V; j++)
            {
                min_psd = AWI_MIN(local_ms_win_jth[i], min_psd);
                local_ms_win_jth += AWI_FRAME_BAND_COUNT;
            }
            local_ms_ns_est[i] = min_psd;
        }

        p->cfg.local_ms_counter++;

        if( p->cfg.local_ms_counter == p->cfg.V)
        {

            p->cfg.local_ms_counter = 0;

            offset = p->cfg.sub_win_counter * AWI_FRAME_BAND_COUNT;
            for(int i =  0; i < AWI_FRAME_BAND_COUNT; i++)
            {
                global_ms_win[offset + i] = local_ms_ns_est[i];
            }

            p->cfg.sub_win_counter++;

            if( p->cfg.sub_win_counter == p->cfg.U)
                p->cfg.sub_win_counter = 0;

            for(int i = 0; i < AWI_FRAME_BAND_COUNT; i++)
            {
                min_psd = 10000.f;
                global_ms_win_jth = global_ms_win;
                for(int j = 0; j < p->cfg.U; j++)
                {
                    min_psd = AWI_MIN(min_psd, global_ms_win_jth[i]);
                    global_ms_win_jth += AWI_FRAME_BAND_COUNT;
                }
                global_ms_ns_est[i] = min_psd;
            }

            if(avg_iv_eqs < 0.03f)
                noise_slope_max = 1.5849f;
            else if(avg_iv_eqs < 0.05f)
                noise_slope_max = 1.4125f;
            else if(avg_iv_eqs < 0.06f)
                noise_slope_max = 1.2589f;
            else
                noise_slope_max = 1.2f;

            for(int i = 0; i < AWI_FRAME_BAND_COUNT; i++)
            {
                if( local_ms_ns_est[i] < noise_slope_max * global_ms_ns_est[i] &&
                    local_ms_ns_est[i] > global_ms_ns_est[i] )
                {
                    global_ms_ns_est[i] = local_ms_ns_est[i];

                    global_ms_win_jth = global_ms_win;
                    for(int j = 0; j < p->cfg.U; j++)
                    {
                        global_ms_win_jth[i] = local_ms_ns_est[i];
                        global_ms_win_jth += AWI_FRAME_BAND_COUNT;
                    }
                }
            }
        }
        else
        {
            p->cfg.alpha_bkg_inc = 0.98;
            if (overall_SNR > 1.2589)
                p->cfg.alpha_bkg_inc = 0.992;
            alpha_bkg_inc_1 = 1 - p->cfg.alpha_bkg_inc;

            if (!enable_vad_auxi || ( enable_vad_auxi && auxi_vad_flag) )
            {
                for(int i = 0; i < p->cfg.high_freq_id_th; i++)
                {
                    min_psd = AWI_MIN(global_ms_ns_est[i], local_ms_ns_est[i]);

                    if ( bkg_est_psd[i] * p->cfg.bkg_update_th < min_psd )
                    {
                        bkg_est_psd[i] = p->cfg.alpha_bkg_inc * bkg_est_psd[i] +
                                        alpha_bkg_inc_1 * min_psd;
                    }
                    else
                    {
                        bkg_est_psd[i] = p->cfg.alpha_bkg_dec * bkg_est_psd[i] +
                                        alpha_bkg_dec_1 * min_psd;
                    }
                    bkg_est_psd[i] = AWI_NSE_ALPHA_UPDATE[i] * bkg_est_psd[i] + ( 1 - AWI_NSE_ALPHA_UPDATE[i] ) * local_ms_ns_est[i];
                }

                for(int i = p->cfg.high_freq_id_th; i < AWI_FRAME_BAND_COUNT; i++)
                {
                    bkg_est_psd[i] = local_ms_ns_est[i];
                }
            }
            else
            {
                for(int i = 0; i < AWI_FRAME_BAND_COUNT; i++)
                {
                    bkg_est_psd[i] = local_ms_ns_est[i];
                }
            }   
        }
    }
}
