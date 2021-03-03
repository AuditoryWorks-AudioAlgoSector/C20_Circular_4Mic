#include <math.h>
#include "awi_utils.h"
#include "awi_subband_aec.h"

void awi_subband_aec_init(awi_subband_aec_part1_t *p1, awi_subband_aec_part2_t *p2)
{
    p1->is_aec_on = 0;
    p1->is_aec_update = 0;
    p1->divergence_flag_fast = 0;
    p1->sum_inst_Far_autoPsd = 0.f;
    p1->sum_recur_Far_autoPsd = 0.f;
    p1->sum_aec_Far_autoPsd = 0.f;
    p1->max_spk_psd = -1000;
    p1->header_counter = 0;
    p1->is_head = 1;
    p1->init_aec_flag = 1;
    p1->init_global_far_psd = 1;
    p1->init_enable_aec_psd = 1;

    memset(p1->xfm,  0,  AWI_FRAME_BAND_COUNT * 2 * p1->cfg.taps * sizeof(float));
    memset(p1->wfb,  0,  AWI_BEAM_CHANNEL * TOTAL_TAP_NUM * 2 * sizeof(float));
    memset(p1->yfm,  0,  AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    memset(p1->nlms_mEk, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    memset(p1->sye,  0,  AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    memset(p1->sy,   0,  AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p1->se,   0,  AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * sizeof(float));

    memset(p1->aec_Near_autoPsd,  0, AWI_BEAM_CHANNEL * AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p1->aec_Far_autoPsd,   0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p1->recur_Far_autoPsd, 0, AWI_FRAME_BAND_COUNT * sizeof(float));
    memset(p1->inst_Far_autoPsd,  0, AWI_FRAME_BAND_COUNT * sizeof(float));

    memset(p1->fullband_spk_mic_erl_chs, 0, AWI_BEAM_CHANNEL * sizeof(float));
    memset(p1->recur_full_band_erle_chs, 0, AWI_BEAM_CHANNEL * sizeof(float));
    memset(p1->recur_aec_eo_fast_energy, 0, AWI_BEAM_CHANNEL * sizeof(float));
    memset(p1->recur_aec_eo_slow_energy, 0, AWI_BEAM_CHANNEL * sizeof(float));
    memset(p1->tmp_yo_slow, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    memset(p1->max_spk_freq_psd, 0, AWI_FRAME_BAND_COUNT * sizeof(float));

    /* for dual filter */
    memset(p2->wfb_fast,  0,  AWI_BEAM_CHANNEL * TOTAL_TAP_NUM * 2 * sizeof(float));

}


void awi_subband_aec_process(awi_subband_aec_part1_t *p, awi_subband_aec_part2_t *p2, float *fbFrameIn,
        float *refFrameIn, float *fbFrameOut)
{
    float applied_alpha_far;
    float tmp_re, tmp_ie;
    float *aec_eo;
    float *tmp_mic_sp, *tmp_mic_psd;
    float alpha_far_1, gamma_1, alpha_eo_energy_1;

    float *tmp_ch_se, *tmp_ch_sy, *tmp_ch_ye, *tmp_yo, *tmp_ch_wfb, *tmp_ch_wfb_fast;
    int idx_re, idx_im;
    float tmp_re_fast, tmp_ie_fast, tmp_re_slow, tmp_ie_slow;
    float tmp_psd_fast, tmp_psd_slow;
    float eo_energy_slow, eo_energy_fast;
    float *tmp_eo_slow, *tmp_yfm_slow;

    int freqbinsby2 = AWI_FRAME_BAND_COUNT * 2;
    size_t single_sp_bytes = freqbinsby2 * sizeof(float);
    size_t single_psd_bytes = AWI_FRAME_BAND_COUNT * sizeof(float);
    size_t weight_bytes = TOTAL_TAP_NUM * 2 * sizeof(float);
    size_t xfm_bytes = AWI_AEC_TAP * AWI_FRAME_BAND_COUNT * 2 * sizeof(float);
    int weight_size = TOTAL_TAP_NUM * 2;

    float *spk_sp_curr = refFrameIn + 0 * AWI_FRAME_BAND_COUNT * 2;
    p->sum_inst_Far_autoPsd = 0;
    /* compute PSD and corresponding sum of spk signal within a certain frequency band */
    for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        /* compute instant PSD of spk signal */
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        tmp_re = spk_sp_curr[idx_re];
        tmp_ie = spk_sp_curr[idx_im];

        p->inst_Far_autoPsd[band_id] = tmp_re * tmp_re + tmp_ie * tmp_ie + AWI_EPS;

        /* compute sum of instant PSD of spk signal*/
        p->sum_inst_Far_autoPsd += p->inst_Far_autoPsd[band_id];
    }


    if ( p->sum_inst_Far_autoPsd > p->sum_recur_Far_autoPsd )
    {
        applied_alpha_far = p->cfg.alpha_far_at;
    }
    else
    {
        applied_alpha_far = p->cfg.alpha_far_ar;
    }

    alpha_eo_energy_1 = 1 - p->cfg.alpha_eo_energy;
    gamma_1 = 1 - p->cfg.gamma;
    p->sum_recur_Far_autoPsd = 0;
    p->max_spk_psd = -1000;
    alpha_far_1 = 1 - applied_alpha_far;
    /* compute recursive PSD for spk signal */
    if (p->init_global_far_psd)
    {
        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            p->recur_Far_autoPsd[band_id] = p->inst_Far_autoPsd[band_id];
            p->max_spk_freq_psd[band_id] = p->inst_Far_autoPsd[band_id];
            p->max_spk_psd = AWI_MAX(p->max_spk_psd, p->max_spk_freq_psd[band_id]);
            /* update sum of recursive PSD of spk signal */
            p->sum_recur_Far_autoPsd += p->recur_Far_autoPsd[band_id];
        }
        p->init_global_far_psd = 0;
    }
    else
    {
        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            p->recur_Far_autoPsd[band_id] = applied_alpha_far * p->recur_Far_autoPsd[band_id] +
                                            alpha_far_1 * p->inst_Far_autoPsd[band_id];
            p->max_spk_freq_psd[band_id] = AWI_MAX(p->inst_Far_autoPsd[band_id], p->recur_Far_autoPsd[band_id]);
            p->max_spk_psd = AWI_MAX(p->max_spk_psd, p->max_spk_freq_psd[band_id]);
            /* update sum of recursive PSD of spk signal */
            p->sum_recur_Far_autoPsd += p->recur_Far_autoPsd[band_id];
        }
    }

    /* aec controlling logic */
    if ( p->sum_recur_Far_autoPsd > p->cfg.spk_energy_th )
    {
        p->is_aec_on = 1;
        if ( p->sum_inst_Far_autoPsd > p->cfg.spk_energy_th )
            p->is_aec_update = 1;
        else
            p->is_aec_update = 0;

        if (p->is_head)
        {
            if ( p->sum_recur_Far_autoPsd > p->cfg.spk_energy_head_th &&
                 p->header_counter < p->cfg.header_counter_th )
            {
                p->header_counter++;
            }
            if (p->header_counter == p->cfg.header_counter_th)
                p->is_head = 0;
        }
    }
    else
    {
        p->is_aec_on = 0;
        p->is_aec_update = 0;
    }

    if ( p->is_aec_on )
    {

        p->sum_aec_Far_autoPsd = 0;
        if (p->init_enable_aec_psd)
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                p->aec_Far_autoPsd[band_id] = p->inst_Far_autoPsd[band_id];
                p->sum_aec_Far_autoPsd += p->aec_Far_autoPsd[band_id];
                p->max_spk_psd = AWI_MAX(p->max_spk_psd, p->aec_Far_autoPsd[band_id]);
            }
        }
        else
        {
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                p->cfg.gamma = 1 - 1.f / FREQ_TAPS[band_id];
                gamma_1 = 1-p->cfg.gamma;
                p->aec_Far_autoPsd[band_id] = p->cfg.gamma * p->aec_Far_autoPsd[band_id] +
                                              gamma_1 * p->inst_Far_autoPsd[band_id];
                p->sum_aec_Far_autoPsd += p->aec_Far_autoPsd[band_id];
                p->max_spk_psd = AWI_MAX(p->max_spk_psd, p->aec_Far_autoPsd[band_id]);
            }
        }


        /* construct spk reference matrix , new spk reference is placed in the foremost */
        int addressOffset = 0;
        if (p->init_aec_flag)
        {
            memcpy(p->xfm + addressOffset, refFrameIn, 2 * AWI_SPK_CHANNEL * freqbinsby2 * sizeof(float));
            p->init_aec_flag = 0;
        }
        else
        {
            memcpy(p->xfm + addressOffset, refFrameIn, AWI_SPK_CHANNEL * freqbinsby2 * sizeof(float));
        }
        /* initialization */
        aec_eo = fbFrameOut;
        tmp_yo = p->yfm;
        tmp_mic_sp = fbFrameIn;
        tmp_mic_psd = p->aec_Near_autoPsd;
        tmp_ch_se = p->se;
        tmp_ch_sy = p->sy;
        tmp_ch_ye = p->sye;
        tmp_ch_wfb = p->wfb;
        tmp_ch_wfb_fast = p2->wfb_fast;
        for (int ch_id = 0; ch_id < AWI_BEAM_CHANNEL; ch_id++)
        {
            tmp_eo_slow = p->nlms_mEk;
            tmp_yfm_slow = p->tmp_yo_slow;

            /* start filtering */
            awi_subband_aec_filter(p->xfm, tmp_ch_wfb, tmp_yfm_slow);
            awi_subband_aec_filter(p->xfm, tmp_ch_wfb_fast, tmp_yo);

            eo_energy_fast = 0, eo_energy_slow = 0;
            for(int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                idx_re = band_id + band_id;
                idx_im = idx_re + 1;

                aec_eo[idx_re] = tmp_re_fast = tmp_mic_sp[idx_re] - tmp_yo[idx_re];
                aec_eo[idx_im] = tmp_ie_fast = tmp_mic_sp[idx_im] - tmp_yo[idx_im];

                tmp_eo_slow[idx_re] = tmp_re_slow = tmp_mic_sp[idx_re] - tmp_yfm_slow[idx_re];
                tmp_eo_slow[idx_im] = tmp_ie_slow = tmp_mic_sp[idx_im] - tmp_yfm_slow[idx_im];

                tmp_psd_fast = tmp_re_fast * tmp_re_fast + tmp_ie_fast * tmp_ie_fast;
                tmp_psd_slow = tmp_re_slow * tmp_re_slow + tmp_ie_slow * tmp_ie_slow;
                eo_energy_fast += tmp_psd_fast;
                eo_energy_slow += tmp_psd_slow;
            }

            p->recur_aec_eo_fast_energy[ch_id] = p->cfg.alpha_eo_energy * p->recur_aec_eo_fast_energy[ch_id] +
                    alpha_eo_energy_1 * eo_energy_fast;
            p->recur_aec_eo_slow_energy[ch_id] = p->cfg.alpha_eo_energy * p->recur_aec_eo_slow_energy[ch_id] +
                    alpha_eo_energy_1 * eo_energy_slow;

            /* updating logic for copying weights */
            if ( p->recur_aec_eo_slow_energy[ch_id] > p->recur_aec_eo_fast_energy[ch_id] * p->cfg.copy_fast_energy_th )
            {
                memcpy(tmp_ch_wfb, tmp_ch_wfb_fast, weight_bytes);
            }
            else if ( p->recur_aec_eo_fast_energy[ch_id] > p->recur_aec_eo_slow_energy[ch_id] * p->cfg.copy_slow_energy_th )
            {
                memcpy(tmp_ch_wfb_fast, tmp_ch_wfb, weight_bytes);
                memcpy(aec_eo, tmp_eo_slow, single_sp_bytes);
                memcpy(tmp_yo, tmp_yfm_slow, single_sp_bytes);
            }

            /* compute frequency-wise variable adaptive step via dtd method */
            if (p->init_enable_aec_psd)
            {
                awi_subband_aec_compute_step_init(p, aec_eo, tmp_yo, tmp_mic_sp, tmp_mic_psd, tmp_ch_se, tmp_ch_sy,
                                             tmp_ch_ye, ch_id);
            }
            else
            {
                awi_subband_aec_compute_step(p, aec_eo, tmp_yo, tmp_mic_sp, tmp_mic_psd, tmp_ch_se, tmp_ch_sy,
                                             tmp_ch_ye, ch_id);
            }

            /* update adaptive filter coefficients */
            if (p->is_aec_update)
            {
                /*
             * check for divergence
             * if divergence happens, restart adaptive filter
             * */
                if (p->divergence_flag_fast)
                {
                    memset(tmp_ch_wfb_fast, 0, weight_bytes);
                    memset(p->xfm, 0, xfm_bytes);
                    memset(p->aec_Far_autoPsd, 0, single_psd_bytes);
                    memset(tmp_mic_psd, 0, single_psd_bytes);
                    memset(tmp_ch_sy, 0, single_psd_bytes);
                    memset(tmp_ch_se, 0, single_psd_bytes);
                    memset(tmp_ch_ye, 0, single_sp_bytes);
                    memcpy(aec_eo, tmp_mic_sp, single_sp_bytes);
                }
                else
                {
                    awi_subband_aec_update_weight(p->xfm, tmp_ch_wfb_fast, p->nlms_mEk);
                }

            }
            /* update */
            aec_eo += freqbinsby2;
            tmp_yo += freqbinsby2;
            tmp_mic_sp += freqbinsby2;
            tmp_mic_psd += AWI_FRAME_BAND_COUNT;
            tmp_ch_se += AWI_FRAME_BAND_COUNT;
            tmp_ch_sy += AWI_FRAME_BAND_COUNT;
            tmp_ch_ye += freqbinsby2;
            tmp_ch_wfb += weight_size;
            tmp_ch_wfb_fast += weight_size;
        }

        if (p->init_enable_aec_psd)
            p->init_enable_aec_psd = 0;

        /* update spk reference matrix */
        addressOffset = p->cfg.spkChannels * freqbinsby2;
        memmove(p->xfm + addressOffset, p->xfm, ( p->cfg.taps - 1 ) * p->cfg.spkChannels * freqbinsby2 * sizeof(float));
    }
    else
    {
        size_t total_beams_bytes = AWI_BEAM_CHANNEL * freqbinsby2 * sizeof(float);
        memcpy(fbFrameOut, fbFrameIn, total_beams_bytes);
        memset(p->yfm,  0,  total_beams_bytes);
        memset(p->fullband_spk_mic_erl_chs, 0, AWI_BEAM_CHANNEL * sizeof(float));
        memset(p->recur_full_band_erle_chs, 0, AWI_BEAM_CHANNEL * sizeof(float));
    }

    // for (int ch = 0; ch < AWI_BEAM_CHANNEL; ch++)
    // {
    //     printf("%.15f %.15f\n", p->fullband_spk_mic_erl_chs[ch], p->recur_full_band_erle_chs[ch]);
    // }

    // printf("%.15f %.15f\n", p->sum_inst_Far_autoPsd, p->sum_recur_Far_autoPsd);

}



void awi_subband_aec_filter(float *input, float *weight, float *output)
{
    float *pSrcA, *pSrcB, *pDst, *pSrcA_Prev;
    unsigned int blkCnt;
    float a1, b1, c1, d1;
    float a2, b2, c2, d2;
    float acc1, acc2, acc3, acc4;
    /* process first tap: initialization */
    /* with full frequency bins */
    pSrcA = input; /* filter input */
    pSrcB = weight; /* filter weight */
    pDst = output; /* to store output of adaptive filter */
    for(int tap = 0; tap < 1; tap++)
    {
        blkCnt = 32;

        while (blkCnt > 0U)
        {
            a1 = *pSrcA;
            c1 = *pSrcB;

            b1 = *(pSrcA + 1);
            acc1 = a1 * c1;

            a2 = *(pSrcA + 2);
            acc2 = (b1 * c1);

            d1 = *(pSrcB + 1);
            c2 = *(pSrcB + 2);
            acc1 -= b1 * d1;

            d2 = *(pSrcB + 3);
            acc3 = a2 * c2;

            b2 = *(pSrcA + 3);
            acc2 += (a1 * d1);


            a1 = *(pSrcA + 4);
            acc4 = (a2 * d2);

            c1 = *(pSrcB + 4);
            acc3 -= (b2 * d2);
            *pDst = acc1;

            b1 = *(pSrcA + 5);
            acc4 += b2 * c2;

            *(pDst + 1) = acc2;
            acc1 = (a1 * c1);

            d1 = *(pSrcB + 5);
            acc2 = (b1 * c1);

            *(pDst + 2) = acc3;
            *(pDst + 3) = acc4;

            a2 = *(pSrcA + 6);
            acc1 -= (b1 * d1);

            c2 = *(pSrcB + 6);
            acc2 += (a1 * d1);

            b2 = *(pSrcA + 7);
            acc3 = (a2 * c2);

            d2 = *(pSrcB + 7);
            acc4 = (b2 * c2);

            *(pDst + 4) = acc1;
            pSrcA += 8U;

            acc3 -= (b2 * d2);
            acc4 += (a2 * d2);

            *(pDst + 5) = acc2;
            pSrcB += 8U;

            *(pDst + 6) = acc3;
            *(pDst + 7) = acc4;

            pDst += 8U;
            blkCnt--;
        }

        blkCnt = 1;

        while(blkCnt > 0U)
        {
            a1 = *pSrcA;
            b1 = *(pSrcA+1);
            c1 = *(pSrcB);
            d1 = *(pSrcB+1);
            *pDst = a1 * c1 - b1 * d1;
            *(pDst+1) = c1 * b1 + a1 * d1;
            blkCnt--;
            pSrcA += 2U;
            pSrcB += 2U;
            pDst += 2U;
        }
    }

    /* process taps with full frequency bins */
    for(int tap = 1; tap < AWI_AEC_TAP_LOWEST; tap++)
    {
        pDst = output; /* to store output of adaptive filter */
        blkCnt = 32;

        while (blkCnt > 0U)
        {
            a1 = *pSrcA;
            c1 = *pSrcB;

            b1 = *(pSrcA + 1);
            acc1 = a1 * c1;

            a2 = *(pSrcA + 2);
            acc2 = (b1 * c1);

            d1 = *(pSrcB + 1);
            c2 = *(pSrcB + 2);
            acc1 -= b1 * d1;

            d2 = *(pSrcB + 3);
            acc3 = a2 * c2;

            b2 = *(pSrcA + 3);
            acc2 += (a1 * d1);

            a1 = *(pSrcA + 4);
            acc4 = (a2 * d2);

            c1 = *(pSrcB + 4);
            acc3 -= (b2 * d2);
            *pDst += acc1;

            b1 = *(pSrcA + 5);
            acc4 += b2 * c2;

            *(pDst + 1) += acc2;
            acc1 = (a1 * c1);

            d1 = *(pSrcB + 5);
            acc2 = (b1 * c1);

            *(pDst + 2) += acc3;
            *(pDst + 3) += acc4;

            a2 = *(pSrcA + 6);
            acc1 -= (b1 * d1);

            c2 = *(pSrcB + 6);
            acc2 += (a1 * d1);

            b2 = *(pSrcA + 7);
            acc3 = (a2 * c2);

            d2 = *(pSrcB + 7);
            acc4 = (b2 * c2);

            *(pDst + 4) += acc1;
            pSrcA += 8U;

            acc3 -= (b2 * d2);
            acc4 += (a2 * d2);

            *(pDst + 5) += acc2;
            pSrcB += 8U;

            *(pDst + 6) += acc3;
            *(pDst + 7) += acc4;

            pDst += 8U;
            blkCnt--;
        }

        blkCnt = 1;
        while(blkCnt > 0U)
        {
            a1 = *pSrcA;
            b1 = *(pSrcA+1);
            c1 = *pSrcB;
            d1 = *(pSrcB+1);
            *pDst += a1 * c1 - b1 * d1;
            *(pDst+1) += c1 * b1 + a1 * d1;
            blkCnt--;
            pSrcA += 2U;
            pSrcB += 2U;
            pDst += 2U;
        }
    }

    /* process residue taps with partial frequency bins */
    pSrcA_Prev = pSrcA;
    for (int tap = AWI_AEC_TAP_LOWEST; tap < AWI_AEC_TAP; tap++)
    {
        pDst = output; /* to store output of adaptive filter */
        pSrcA = pSrcA_Prev;
        for (int band_id = 0; band_id < TAP_FREQS[tap]; band_id++)
        {
            a1 = *pSrcA;
            b1 = *(pSrcA+1);
            c1 = *pSrcB;
            d1 = *(pSrcB+1);
            *pDst += a1 * c1 - b1 * d1;
            *(pDst+1) += c1 * b1 + a1 * d1;
            pSrcA += 2U;
            pSrcB += 2U;
            pDst += 2U;
        }
        pSrcA_Prev += AWI_FRAME_BAND_COUNT * 2;
    }
}



void awi_subband_aec_update_weight(float *input, float *weight, float *nlms_eo)
{
    float *pSrcA, *pSrcB, *pDst, *pSrcA_Prev;
    unsigned int blkCnt;
    float a1, b1, c1, d1;
    float a2, b2, c2, d2;
    float acc1, acc2, acc3, acc4;
    pSrcA = input;
    pDst = weight;
    /* process taps with full frequency bins */
    for (int tap = 0; tap < AWI_AEC_TAP_LOWEST; tap++)
    {
        pSrcB = nlms_eo;

        blkCnt = 32;

        while (blkCnt > 0U) {
            a1 = *pSrcA;
            c1 = *pSrcB;

            b1 = -*(pSrcA + 1);
            acc1 = a1 * c1;

            a2 = *(pSrcA + 2);
            acc2 = (b1 * c1);

            d1 = *(pSrcB + 1);
            c2 = *(pSrcB + 2);
            acc1 -= b1 * d1;

            d2 = *(pSrcB + 3);
            acc3 = a2 * c2;

            b2 = -*(pSrcA + 3);
            acc2 += (a1 * d1);

            a1 = *(pSrcA + 4);
            acc4 = (a2 * d2);

            c1 = *(pSrcB + 4);
            acc3 -= (b2 * d2);
            *pDst += acc1;

            b1 = -*(pSrcA + 5);
            acc4 += b2 * c2;

            *(pDst + 1) += acc2;
            acc1 = (a1 * c1);

            d1 = *(pSrcB + 5);
            acc2 = (b1 * c1);

            *(pDst + 2) += acc3;
            *(pDst + 3) += acc4;

            a2 = *(pSrcA + 6);
            acc1 -= (b1 * d1);

            c2 = *(pSrcB + 6);
            acc2 += (a1 * d1);

            b2 = -*(pSrcA + 7);
            acc3 = (a2 * c2);

            d2 = *(pSrcB + 7);
            acc4 = (b2 * c2);

            *(pDst + 4) += acc1;
            pSrcA += 8U;

            acc3 -= (b2 * d2);
            acc4 += (a2 * d2);

            *(pDst + 5) += acc2;
            pSrcB += 8U;

            *(pDst + 6) += acc3;
            *(pDst + 7) += acc4;

            pDst += 8U;

            blkCnt--;
        }

        blkCnt = 1;

        while (blkCnt > 0U) {
            a1 = *pSrcA;
            b1 = -*(pSrcA+1);
            c1 = *pSrcB;
            d1 = *(pSrcB+1);

            *pDst += (a1 * c1) - (b1 * d1);
            *(pDst+1) += (b1 * c1) + (a1 * d1);

            blkCnt--;
            pSrcA += 2U;
            pSrcB += 2U;
            pDst += 2U;
        }

    }
    /* process residue taps with partial frequency bins */
    pSrcA_Prev = pSrcA;
    for (int tap = AWI_AEC_TAP_LOWEST; tap < AWI_AEC_TAP; tap++)
    {
        pSrcB = nlms_eo;
        pSrcA = pSrcA_Prev;
        for (int band_id = 0; band_id < TAP_FREQS[tap]; band_id++)
        {
            a1 = *pSrcA;
            b1 = -*(pSrcA+1);
            c1 = *pSrcB;
            d1 = *(pSrcB+1);

            *pDst += (a1 * c1) - (b1 * d1);
            *(pDst+1) += (b1 * c1) + (a1 * d1);

            pSrcA += 2U;
            pSrcB += 2U;
            pDst += 2U;
        }
        pSrcA_Prev += AWI_FRAME_BAND_COUNT * 2;
    }
}


void awi_subband_aec_compute_step_init(awi_subband_aec_part1_t *p, float *aec_eo, float *tmp_yo, float *tmp_mic_sp,
                                  float *tmp_mic_psd, float *tmp_ch_se, float *tmp_ch_sy, float *tmp_ch_ye,
                                  int ch_id)
{
    int idx_re, idx_im;
    float tmp_yo_re, tmp_yo_ie, tmp_mic_re, tmp_mic_im, tmp_eo_re, tmp_eo_ie;
    float tmp_psd_mic, tmp_psd_eo, tmp_psd;
    float inst_aec_Near_autoPsd_Sum, aec_Near_autoPsd_Sum, inst_aec_Eo_autoPsd_sum, aec_Eo_autoPsd_sum;
    float ch_ye_re, ch_ye_im, ch_sy_ri, ch_se_ri;
    float coh_yo_eo;
    float dtd_fast, norm_factor_fast;
    float factor_fast_1 = p->cfg.factor_fast - 1;
    float err_step_size_scale;
    float reg_factor, dtd_reg_factor;
    float inst_fullband_erle, inst_fullband_spk_mic_erl;

    reg_factor = p->cfg.reg_factor_basis;
    dtd_reg_factor = p->cfg.dtd_reg_factor_basis;
    err_step_size_scale = p->cfg.err_step_size_scale;

    inst_aec_Near_autoPsd_Sum = 0;
    aec_Near_autoPsd_Sum = 0;
    inst_aec_Eo_autoPsd_sum = 0;
    aec_Eo_autoPsd_sum = 0;
    for(int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        if (p->is_head)
        {
            reg_factor = p->cfg.reg_factor_basis_head;
            dtd_reg_factor = p->cfg.dtd_reg_factor_basis_head;
        }

        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        tmp_yo_re = tmp_yo[idx_re];
        tmp_yo_ie = tmp_yo[idx_im];
        tmp_mic_re = tmp_mic_sp[idx_re];
        tmp_mic_im = tmp_mic_sp[idx_im];

        /* output after linear aec */
        tmp_eo_re = aec_eo[idx_re];
        tmp_eo_ie = aec_eo[idx_im];

        /* compute PSD for mic signal and its sum */
        tmp_psd_mic = tmp_mic_re * tmp_mic_re + tmp_mic_im * tmp_mic_im;
        inst_aec_Near_autoPsd_Sum += tmp_psd_mic;
        tmp_mic_psd[band_id] = tmp_psd_mic;
        aec_Near_autoPsd_Sum += tmp_mic_psd[band_id];

        /* compute PSD for output signal and its sum */
        tmp_psd_eo = tmp_eo_re * tmp_eo_re + tmp_eo_ie * tmp_eo_ie;
        inst_aec_Eo_autoPsd_sum += tmp_psd_eo;
        tmp_ch_se[band_id] = tmp_psd_eo;
        aec_Eo_autoPsd_sum += tmp_ch_se[band_id];

        /* compute PSD for output of adaptive filter */
        tmp_psd = tmp_yo_re * tmp_yo_re + tmp_yo_ie * tmp_yo_ie;
        tmp_ch_sy[band_id] = tmp_psd;


        /* compute cross PSD of aec output and adaptive filter output */
        if ( p->aec_Far_autoPsd[band_id] > p->cfg.coh_th )
        {
            tmp_ch_ye[idx_re] = ( tmp_eo_re * tmp_yo_re  + tmp_eo_ie * tmp_yo_ie );
            tmp_ch_ye[idx_im] = ( tmp_eo_ie * tmp_yo_re - tmp_eo_re * tmp_yo_ie );
        }


        ch_ye_re = tmp_ch_ye[idx_re];
        ch_ye_im = tmp_ch_ye[idx_im];
        ch_sy_ri = tmp_ch_sy[band_id];
        ch_se_ri = tmp_ch_se[band_id];

        coh_yo_eo = (ch_ye_re * ch_ye_re + ch_ye_im * ch_ye_im) / \
                                      ( ch_sy_ri * ch_se_ri + p->cfg.coh_reg_factor );
        if(coh_yo_eo < 0)
            coh_yo_eo = 0;
        else if(coh_yo_eo > 0.999f)
            coh_yo_eo = 0.999f;

        /* compute adaptive step size for fast adaptive filter */
        dtd_fast  = ( ( 1 + factor_fast_1 * coh_yo_eo ) *
                      (ch_sy_ri + dtd_reg_factor ) ) / \
                                 ( err_step_size_scale * ( 1 - coh_yo_eo ) * \
                                  ch_se_ri + ch_sy_ri + dtd_reg_factor );
        if ( p->aec_Far_autoPsd[band_id] < p->cfg.sig_th )
            dtd_fast = 0;
        norm_factor_fast = dtd_fast  * p->cfg.mufb_fast / ( FREQ_TAPS[band_id] * p->aec_Far_autoPsd[band_id] + reg_factor );

        p->nlms_mEk[idx_re] = norm_factor_fast * tmp_eo_re;
        p->nlms_mEk[idx_im] = norm_factor_fast * tmp_eo_ie;
    }
    
    p->divergence_flag_fast = 0;
    if (aec_Eo_autoPsd_sum > aec_Near_autoPsd_Sum * p->cfg.divergT)
    {
        p->divergence_flag_fast = 1;
    }

    inst_fullband_erle = ( aec_Near_autoPsd_Sum + 1e-20f ) / ( aec_Eo_autoPsd_sum + 1e-20f );
    inst_fullband_spk_mic_erl = ( p->sum_aec_Far_autoPsd + 1e-20f ) / ( aec_Near_autoPsd_Sum + 1e-20f );
    p->recur_full_band_erle_chs[ch_id] = inst_fullband_erle;
    p->fullband_spk_mic_erl_chs[ch_id] = inst_fullband_spk_mic_erl;

}



void awi_subband_aec_compute_step(awi_subband_aec_part1_t *p, float *aec_eo, float *tmp_yo, float *tmp_mic_sp,
                                  float *tmp_mic_psd, float *tmp_ch_se, float *tmp_ch_sy, float *tmp_ch_ye,
                                  int ch_id)
{
    int idx_re, idx_im;
    float tmp_yo_re, tmp_yo_ie, tmp_mic_re, tmp_mic_im, tmp_eo_re, tmp_eo_ie;
    float tmp_psd_mic, tmp_psd_eo, tmp_psd;
    float inst_aec_Near_autoPsd_Sum, aec_Near_autoPsd_Sum, inst_aec_Eo_autoPsd_sum, aec_Eo_autoPsd_sum;
    float ch_ye_re, ch_ye_im, ch_sy_ri, ch_se_ri;
    float coh_yo_eo;
    float dtd_fast, norm_factor_fast;
    float stat_ff_1 = 1-p->cfg.stat_ff;
    float factor_fast_1 = p->cfg.factor_fast - 1;
    float err_step_size_scale;
    float reg_factor, dtd_reg_factor;
    float inst_fullband_erle, inst_fullband_spk_mic_erl;

    reg_factor = p->cfg.reg_factor_basis;
    dtd_reg_factor = p->cfg.dtd_reg_factor_basis;
    err_step_size_scale = p->cfg.err_step_size_scale;

    inst_aec_Near_autoPsd_Sum = 0;
    aec_Near_autoPsd_Sum = 0;
    inst_aec_Eo_autoPsd_sum = 0;
    aec_Eo_autoPsd_sum = 0;
    for(int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        p->cfg.stat_ff = 1 - 1.f / 2.f / FREQ_TAPS[band_id];
        stat_ff_1 = 1-p->cfg.stat_ff;

        if (p->is_head)
        {
            reg_factor = p->cfg.reg_factor_basis_head;
            dtd_reg_factor = p->cfg.dtd_reg_factor_basis_head;
        }

        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        tmp_yo_re = tmp_yo[idx_re];
        tmp_yo_ie = tmp_yo[idx_im];
        tmp_mic_re = tmp_mic_sp[idx_re];
        tmp_mic_im = tmp_mic_sp[idx_im];

        /* output after linear aec */
        tmp_eo_re = aec_eo[idx_re];
        tmp_eo_ie = aec_eo[idx_im];

        tmp_psd_eo = tmp_eo_re * tmp_eo_re + tmp_eo_ie * tmp_eo_ie;
        tmp_psd_mic = tmp_mic_re * tmp_mic_re + tmp_mic_im * tmp_mic_im;

        if (tmp_psd_eo > tmp_psd_mic)
        {
            aec_eo[idx_re] = tmp_mic_re;
            aec_eo[idx_im] = tmp_mic_im;
            tmp_yo_re = 0;
            tmp_yo_ie = 0;
        }

        /* compute PSD for mic signal and its sum */
        tmp_psd_mic = tmp_mic_re * tmp_mic_re + tmp_mic_im * tmp_mic_im;
        inst_aec_Near_autoPsd_Sum += tmp_psd_mic;
        tmp_mic_psd[band_id] = p->cfg.stat_ff * tmp_mic_psd[band_id] + stat_ff_1 * tmp_psd_mic;
        aec_Near_autoPsd_Sum += tmp_mic_psd[band_id];

        /* compute PSD for output signal and its sum */
        tmp_psd_eo = tmp_eo_re * tmp_eo_re + tmp_eo_ie * tmp_eo_ie;
        inst_aec_Eo_autoPsd_sum += tmp_psd_eo;
        tmp_ch_se[band_id] = p->cfg.stat_ff * tmp_ch_se[band_id] + stat_ff_1 * tmp_psd_eo;
        aec_Eo_autoPsd_sum += tmp_ch_se[band_id];

        /* compute PSD for output of adaptive filter */
        tmp_psd = tmp_yo_re * tmp_yo_re + tmp_yo_ie * tmp_yo_ie;
        tmp_ch_sy[band_id] = p->cfg.stat_ff * tmp_ch_sy[band_id] + stat_ff_1 * tmp_psd;


        /* compute cross PSD of aec output and adaptive filter output */
        if ( p->aec_Far_autoPsd[band_id] > p->cfg.coh_th )
        {
            tmp_ch_ye[idx_re] = p->cfg.stat_ff * tmp_ch_ye[idx_re] +
                                stat_ff_1 * ( tmp_eo_re * tmp_yo_re  + tmp_eo_ie * tmp_yo_ie );
            tmp_ch_ye[idx_im] = p->cfg.stat_ff * tmp_ch_ye[idx_im] +
                                stat_ff_1 * ( tmp_eo_ie * tmp_yo_re - tmp_eo_re * tmp_yo_ie );
        }


        ch_ye_re = tmp_ch_ye[idx_re];
        ch_ye_im = tmp_ch_ye[idx_im];
        ch_sy_ri = tmp_ch_sy[band_id];
        ch_se_ri = tmp_ch_se[band_id];

        coh_yo_eo = (ch_ye_re * ch_ye_re + ch_ye_im * ch_ye_im) / \
                                      ( ch_sy_ri * ch_se_ri + p->cfg.coh_reg_factor );
        if(coh_yo_eo < 0)
            coh_yo_eo = 0;
        else if(coh_yo_eo > 0.999f)
            coh_yo_eo = 0.999f;

        /* compute adaptive step size for fast adaptive filter */
        dtd_fast  = ( ( 1 + factor_fast_1 * coh_yo_eo ) *
                      (ch_sy_ri + dtd_reg_factor ) ) / \
                                 ( err_step_size_scale * ( 1 - coh_yo_eo ) * \
                                  ch_se_ri + ch_sy_ri + dtd_reg_factor );
        if ( p->aec_Far_autoPsd[band_id] < p->cfg.sig_th)
            dtd_fast = 0;
        norm_factor_fast = dtd_fast  * p->cfg.mufb_fast / ( FREQ_TAPS[band_id] * p->aec_Far_autoPsd[band_id] + reg_factor );

        p->nlms_mEk[idx_re] = norm_factor_fast * tmp_eo_re;
        p->nlms_mEk[idx_im] = norm_factor_fast * tmp_eo_ie;
    }

    p->divergence_flag_fast = 0;
    if (aec_Eo_autoPsd_sum > aec_Near_autoPsd_Sum * p->cfg.divergT)
    {
        p->divergence_flag_fast = 1;
    }

    inst_fullband_erle = ( aec_Near_autoPsd_Sum + 1e-20f ) / ( aec_Eo_autoPsd_sum + 1e-20f );
    inst_fullband_spk_mic_erl = ( p->sum_aec_Far_autoPsd + 1e-20f ) / ( aec_Near_autoPsd_Sum + 1e-20f );
    p->recur_full_band_erle_chs[ch_id] = inst_fullband_erle;
    p->fullband_spk_mic_erl_chs[ch_id] = inst_fullband_spk_mic_erl;

}
