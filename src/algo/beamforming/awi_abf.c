#include <math.h>
#include <string.h>
#include "awi_utils.h"
#include "awi_abf.h"

void awi_abf_init(awi_abf_part1_t *p1, awi_abf_part2_t *p2)
{
    awi_abf_reset(p1, p2);
}

void awi_abf_reset(awi_abf_part1_t *p1, awi_abf_part2_t *p2)
{
    memset(p1->wfb, 0, sizeof(float) * AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * AWI_DOA_BIN_NUM);
    memset(p2->xfm, 0, sizeof(float) * AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * AWI_DOA_BIN_NUM);
    memset(p1->sd,  0, AWI_FRAME_BAND_COUNT * AWI_DOA_BIN_NUM * sizeof(float));
    memset(p1->se,  0, AWI_FRAME_BAND_COUNT * AWI_DOA_BIN_NUM * sizeof(float));
    memset(p1->sx,  0, AWI_FRAME_BAND_COUNT * AWI_ABF_REFERENCE * AWI_DOA_BIN_NUM * sizeof(float));
    memset(p2->mEk, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));

    memset(p2->yfm, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));

    float *snr, *global_snr;
    for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
    {
        for (int ref_id = 0; ref_id < AWI_ABF_REFERENCE; ref_id++)
        {
            snr = p1->recur_snr + ref_id * AWI_FRAME_BAND_COUNT + beam_id * AWI_ABF_REFERENCE * AWI_FRAME_BAND_COUNT;
            global_snr = p1->recur_global_snr + beam_id * AWI_ABF_REFERENCE;
            global_snr[ref_id] = p1->cfg.midPoint;
            for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
            {
                snr[band_id] = p1->cfg.midPoint;
            }
        }
    }

}


float Sigmoid(float val, float midPoint, float slope)
{
    float step  = (val - midPoint) * slope;
    float coeff = 1 - 0.5f * ( 1 + step / (1 + fabsf(step)));
    return fmaxf_local(fminf_local(coeff, 1), 0);
}


void awi_abf_process(awi_abf_part1_t *p1, awi_abf_part2_t *p2, float *frameIn, float *ref1_sp, float *ref2_sp,
        float *ref3_sp, float *frameOut, int beamId)
{
    /* construct reference , place newly coming references into the foremost location*/
    int offset = beamId * AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2;
    int psd_offset = beamId * AWI_FRAME_BAND_COUNT;
    int psd_offset_sx = beamId * AWI_FRAME_BAND_COUNT * AWI_ABF_REFERENCE;
    int freqbinsby2 = AWI_FRAME_BAND_COUNT * 2;
    memcpy(p2->xfm + offset, ref1_sp, freqbinsby2 * sizeof(float));
    memcpy(p2->xfm + offset + freqbinsby2, ref2_sp, freqbinsby2 * sizeof(float));
    memcpy(p2->xfm + offset + 2 * freqbinsby2, ref3_sp, freqbinsby2 * sizeof(float));


    /* filtering */
    float *tmp_xfm, *tmp_wfb;
    float xfm_re, xfm_ie;
    float  *yy = p2->yfm;

    float *pSrcA = p2->xfm + offset;
    float *pSrcB = p1->wfb + offset;
    float *pDst = yy;

    unsigned int blkCnt = AWI_FRAME_BAND_COUNT >> 2U;

    float a1, b1, c1, d1;
    float a2, b2, c2, d2;
    float acc1, acc2, acc3, acc4;

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

    blkCnt = AWI_FRAME_BAND_COUNT % 0x4U;

    while(blkCnt > 0U)
    {
        a1 = *pSrcA++;
        b1 = *pSrcA++;
        c1 = *pSrcB++;
        d1 = *pSrcB++;
        *pDst++ = a1 * c1 - b1 * d1;
        *pDst++ = c1 * b1 + a1 * d1;
        blkCnt--;
    }


    for(int tap = 1; tap < p1->cfg.AFlen; tap++)
    {
        pSrcA = p2->xfm + offset + tap * freqbinsby2;
        pSrcB = p1->wfb + offset + tap * freqbinsby2;
        pDst = yy;

        blkCnt = AWI_FRAME_BAND_COUNT >> 2U;

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

        blkCnt = AWI_FRAME_BAND_COUNT % 0x4U;
        while(blkCnt > 0U)
        {
            a1 = *pSrcA++;
            b1 = *pSrcA++;
            c1 = *pSrcB++;
            d1 = *pSrcB++;
            *pDst++ += a1 * c1 - b1 * d1;
            *pDst++ += c1 * b1 + a1 * d1;
            blkCnt--;
        }
    }


    float tmp_psd;
    float *tmp_sx1 = p1->sx + psd_offset_sx;
    float *tmp_sx2 = p1->sx + psd_offset_sx + AWI_FRAME_BAND_COUNT;
    float *tmp_sx3 = p1->sx + psd_offset_sx + 2 * AWI_FRAME_BAND_COUNT;

    float *snr0 = p1->recur_snr + beamId * AWI_ABF_REFERENCE * AWI_FRAME_BAND_COUNT;
    float *snr1 = p1->recur_snr + AWI_FRAME_BAND_COUNT + beamId * AWI_ABF_REFERENCE * AWI_FRAME_BAND_COUNT;
    float *snr2 = p1->recur_snr + 2 * AWI_FRAME_BAND_COUNT + beamId * AWI_ABF_REFERENCE * AWI_FRAME_BAND_COUNT;
    float *global_snr = p1->recur_global_snr + beamId * AWI_ABF_REFERENCE;


    float *mEk = p2->mEk;

    float frameIn_re, frameIn_im;
    float frameOut_re, frameOut_im;
    float yy_re, yy_im, sp_re, sp_im;
    float sd, se, sx1, sx2, sx3;
    float ref1, ref2, ref3;
    float sum_sd = 0, sum_se = 0, sum_sx1 = 0, sum_sx2 = 0, sum_sx3 = 0;
    float geo_avg_snr = 0;
    float tmp_vss;
    int idx_re, idx_im;


    for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;
        frameIn_re = frameIn[idx_re];
        frameIn_im = frameIn[idx_im];
        yy_re = yy[idx_re];
        yy_im = yy[idx_im];

        frameOut[idx_re] = frameOut_re = frameIn_re - yy_re;
        frameOut[idx_im] = frameOut_im = frameIn_im - yy_im;

        tmp_psd = frameOut_re * frameOut_re + frameOut_im * frameOut_im;
        p1->se[band_id + psd_offset] = se = ( 1 - p1->cfg.alp ) * p1->se[band_id + psd_offset] + p1->cfg.alp * tmp_psd;
        sum_se += se;

        tmp_psd = frameIn_re * frameIn_re + frameIn_im * frameIn_im;
        p1->sd[band_id + psd_offset] = sd = ( 1 - p1->cfg.alp ) * p1->sd[band_id + psd_offset]
                + p1->cfg.alp * tmp_psd;

        sum_sd += sd;

        sp_re = ref1_sp[idx_re];
        sp_im = ref1_sp[idx_im];
        tmp_psd = sp_re * sp_re + sp_im * sp_im;
        tmp_sx1[band_id] = sx1 = ( 1 - p1->cfg.alp ) * tmp_sx1[band_id] + p1->cfg.alp * tmp_psd;

        sum_sx1 += sx1;

        sp_re = ref2_sp[idx_re];
        sp_im = ref2_sp[idx_im];
        tmp_psd = sp_re * sp_re + sp_im * sp_im;
        tmp_sx2[band_id] = sx2 = ( 1 - p1->cfg.alp ) * tmp_sx2[band_id] + p1->cfg.alp * tmp_psd;

        sum_sx2 += sx2;

        sp_re = ref3_sp[idx_re];
        sp_im = ref3_sp[idx_im];
        tmp_psd = sp_re * sp_re + sp_im * sp_im;
        tmp_sx3[band_id] = sx3 = ( 1 - p1->cfg.alp ) * tmp_sx3[band_id] + p1->cfg.alp * tmp_psd;

        sum_sx3 += sx3;

        ref1 = ( sd + p1->cfg.vss_div_factor ) / ( sx1 + p1->cfg.vss_div_factor );

        if(ref1 > p1->cfg.maxSNR)
            ref1 = p1->cfg.maxSNR;
        else if(ref1 < p1->cfg.minSNR)
            ref1 = p1->cfg.minSNR;

        snr0[band_id] = ( 1 - p1->cfg.vss_stat_factor ) * snr0[band_id] + p1->cfg.vss_stat_factor * ref1;

        ref2 = ( sd + p1->cfg.vss_div_factor ) / ( sx2 + p1->cfg.vss_div_factor );

        if(ref2 > p1->cfg.maxSNR)
            ref2 = p1->cfg.maxSNR;
        else if(ref2 < p1->cfg.minSNR)
            ref2 = p1->cfg.minSNR;

        snr1[band_id] = ( 1 - p1->cfg.vss_stat_factor ) * snr1[band_id] + p1->cfg.vss_stat_factor * ref2;

        ref3 = ( sd + p1->cfg.vss_div_factor ) / ( sx3 + p1->cfg.vss_div_factor );

        if(ref3 > p1->cfg.maxSNR)
            ref3 = p1->cfg.maxSNR;
        else if(ref3 < p1->cfg.minSNR)
            ref3 = p1->cfg.minSNR;

        snr2[band_id] = ( 1 - p1->cfg.vss_stat_factor ) * snr2[band_id] + p1->cfg.vss_stat_factor * ref3;

        mEk[idx_re] = p1->cfg.mufb * frameOut_re;
        mEk[idx_im] = p1->cfg.mufb * frameOut_im;
    }


    global_snr[0] = ( 1 - p1->cfg.vss_stat_factor ) * global_snr[0] + p1->cfg.vss_stat_factor * ( sum_sd + p1->cfg.vss_div_factor ) / ( sum_sx1 + p1->cfg.vss_div_factor );
    if(global_snr[0] > p1->cfg.maxSNR)
        global_snr[0] = p1->cfg.maxSNR;
    else if(global_snr[0] < p1->cfg.minSNR)
        global_snr[0] = p1->cfg.minSNR;

    global_snr[1] = ( 1 - p1->cfg.vss_stat_factor ) * global_snr[1] + p1->cfg.vss_stat_factor * ( sum_sd + p1->cfg.vss_div_factor ) / ( sum_sx2 + p1->cfg.vss_div_factor );
    if(global_snr[1] > p1->cfg.maxSNR)
        global_snr[1] = p1->cfg.maxSNR;
    else if(global_snr[1] < p1->cfg.minSNR)
        global_snr[1] = p1->cfg.minSNR;

    global_snr[2] = ( 1 - p1->cfg.vss_stat_factor ) * global_snr[2] + p1->cfg.vss_stat_factor * ( sum_sd + p1->cfg.vss_div_factor ) / ( sum_sx3 + p1->cfg.vss_div_factor );
    if(global_snr[2] > p1->cfg.maxSNR)
        global_snr[2] = p1->cfg.maxSNR;
    else if(global_snr[2] < p1->cfg.minSNR)
        global_snr[2] = p1->cfg.minSNR;



    // decide if divergence happens
    if ( sum_se > p1->cfg.divergence_th * sum_sd )
    {
        memcpy(frameOut, frameIn, freqbinsby2 * sizeof(float));
        memset(p2->xfm + offset, 0, AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * sizeof(float));
        memset(p1->wfb + offset, 0, AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * sizeof(float));
        memset(p2->mEk, 0, freqbinsby2 * sizeof(float));
    }


    tmp_sx1 = p1->sx + psd_offset_sx;
    tmp_sx2 = p1->sx + psd_offset_sx + AWI_FRAME_BAND_COUNT;
    tmp_sx3 = p1->sx + psd_offset_sx + 2 * AWI_FRAME_BAND_COUNT;

    float mEk_re, mEk_ie;

    float vssf;

    for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
    {
        idx_re = band_id + band_id;
        idx_im = idx_re + 1;

        mEk_re = mEk[idx_re];
        mEk_ie = mEk[idx_im];

        sx1 = tmp_sx1[band_id];
        geo_avg_snr = p1->cfg.alpha_snr * snr0[band_id] + ( 1 - p1->cfg.alpha_snr ) * global_snr[0];
        tmp_vss = Sigmoid(geo_avg_snr, p1->cfg.midPoint, p1->cfg.slope);
        vssf = tmp_vss / ( p1->cfg.AFlen * sx1 + p1->cfg.power_reg_factor );

        for (int tap_id = 0; tap_id < p1->cfg.AFlen; tap_id += AWI_ABF_REFERENCE)
        {
            tmp_wfb = p1->wfb + offset + tap_id * freqbinsby2;
            tmp_xfm = p2->xfm + offset + tap_id * freqbinsby2;
            xfm_re = tmp_xfm[idx_re];
            xfm_ie = tmp_xfm[idx_im];

            tmp_wfb[idx_re] = p1->cfg.leaking_factor * tmp_wfb[idx_re] + ( mEk_re * xfm_re + mEk_ie * xfm_ie) * vssf;
            tmp_wfb[idx_im] = p1->cfg.leaking_factor * tmp_wfb[idx_im] + ( -mEk_re * xfm_ie + mEk_ie * xfm_re) * vssf;

        }

        sx2 = tmp_sx2[band_id];
        geo_avg_snr = p1->cfg.alpha_snr * snr1[band_id] + ( 1 - p1->cfg.alpha_snr ) * global_snr[1];
        tmp_vss = Sigmoid(geo_avg_snr, p1->cfg.midPoint, p1->cfg.slope);
        vssf = tmp_vss / ( p1->cfg.AFlen * sx2 + p1->cfg.power_reg_factor );
        for (int tap_id = 1; tap_id < p1->cfg.AFlen; tap_id += AWI_ABF_REFERENCE)
        {
            tmp_wfb = p1->wfb + offset + tap_id * freqbinsby2;
            tmp_xfm = p2->xfm + offset + tap_id * freqbinsby2;
            xfm_re = tmp_xfm[idx_re];
            xfm_ie = tmp_xfm[idx_im];

            tmp_wfb[idx_re] = p1->cfg.leaking_factor * tmp_wfb[idx_re] +  ( mEk_re * xfm_re + mEk_ie * xfm_ie) * vssf;
            tmp_wfb[idx_im] = p1->cfg.leaking_factor * tmp_wfb[idx_im] + ( -mEk_re * xfm_ie + mEk_ie * xfm_re) * vssf;

        }

        sx3 = tmp_sx3[band_id];
        geo_avg_snr = p1->cfg.alpha_snr * snr2[band_id] + ( 1 - p1->cfg.alpha_snr ) * global_snr[2];
        tmp_vss = Sigmoid(geo_avg_snr, p1->cfg.midPoint, p1->cfg.slope);
        vssf = tmp_vss / ( p1->cfg.AFlen * sx3 + p1->cfg.power_reg_factor );
        for (int tap_id = 2; tap_id < p1->cfg.AFlen; tap_id += AWI_ABF_REFERENCE)
        {
            tmp_wfb = p1->wfb + offset + tap_id * freqbinsby2;
            tmp_xfm = p2->xfm + offset + tap_id * freqbinsby2;
            xfm_re = tmp_xfm[idx_re];
            xfm_ie = tmp_xfm[idx_im];

            tmp_wfb[idx_re] = p1->cfg.leaking_factor * tmp_wfb[idx_re] +  ( mEk_re * xfm_re + mEk_ie * xfm_ie) * vssf;
            tmp_wfb[idx_im] = p1->cfg.leaking_factor * tmp_wfb[idx_im] + ( -mEk_re * xfm_ie + mEk_ie * xfm_re) * vssf;

        }

    }

    /* update reference */
    memmove(p2->xfm + offset + AWI_ABF_REFERENCE * freqbinsby2, p2->xfm + offset,
            ( p1->cfg.AFlen - AWI_ABF_REFERENCE ) * freqbinsby2 * sizeof(float));

}

