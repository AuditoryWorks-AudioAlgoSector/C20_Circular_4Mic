//
// Created by xzc on 20-4-8.
//

#include "awi_beam_fusion.h"
#include <string.h>
#include <math.h>

void awi_bf_init(awi_bf_t *p, bool complex_weight)
{
    p->complex_weighting = complex_weight;
    p->direction_id = 0u;
    memset(p->inst_auto_psd, 0, AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * sizeof(float));
    memset(p->recur_auto_psd, 0, AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * sizeof(float));
    memset(p->recur_energy, 0, AWI_BEAM_CHANNEL * sizeof(float));
    memset(p->beam_weight, 0, AWI_BEAM_CHANNEL * sizeof(float));
    for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
        p->recur_beam_weight[beam_id] = 1.f / ( AWI_BEAM_CHANNEL * 1.f );
}


void awi_bf_process(awi_bf_t *p, float *ref_sp, float *input_sp, float *output_sp)
{
    float sp_re, sp_im, tmp_psd, total_energy, weight, norm_weight, max_weight;
    float alpha_psd_at_1, alpha_psd_ar_1, alpha_weight_at_1, alpha_weight_ar_1;
    int id_re, id_im;
    float *bf_auto_psd_beam, *bf_sp_beam, *inst_psd_beam;
    const int sp_offset = 2 * AWI_FRAME_BAND_COUNT;
    const int psd_offset = AWI_FRAME_BAND_COUNT;
    bf_auto_psd_beam = p->recur_auto_psd;
    inst_psd_beam = p->inst_auto_psd;
    bf_sp_beam = input_sp;
    alpha_psd_at_1 = 1 - p->cfg.alpha_psd_at;
    alpha_psd_ar_1 = 1 - p->cfg.alpha_psd_ar;
    alpha_weight_at_1 = 1 - p->cfg.alpha_weight_at;
    alpha_weight_ar_1 = 1 - p->cfg.alpha_weight_ar;
    total_energy = 0;
    for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
    {
        p->recur_energy[beam_id] = 0;
        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            id_re = band_id + band_id;
            id_im = id_re + 1;
            sp_re = bf_sp_beam[id_re];
            sp_im = bf_sp_beam[id_im];
            tmp_psd = sp_re * sp_re + sp_im * sp_im;
            inst_psd_beam[band_id] = tmp_psd;
            if ( tmp_psd > bf_auto_psd_beam[band_id] )
            {
                bf_auto_psd_beam[band_id] = p->cfg.alpha_psd_at * bf_auto_psd_beam[band_id] +
                        alpha_psd_at_1 * tmp_psd;
            }
            else
            {
                bf_auto_psd_beam[band_id] = p->cfg.alpha_psd_ar * bf_auto_psd_beam[band_id] +
                        alpha_psd_ar_1 * tmp_psd;
            }
            p->recur_energy[beam_id] += bf_auto_psd_beam[band_id];
        }
        total_energy += p->recur_energy[beam_id];
        bf_auto_psd_beam += psd_offset;
        inst_psd_beam += psd_offset;
        bf_sp_beam += sp_offset;
    }

    memset(output_sp, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));

    norm_weight = 0;
    max_weight = -1;
    p->direction_id = 0u;
    for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
    {
        weight = p->recur_energy[beam_id] / ( total_energy + AWI_EPS ) - 0.125f;
        weight = weight > 0 ? weight : 0;
        weight *= weight * weight;
        norm_weight += weight;
        p->beam_weight[beam_id] = weight;

        if ( max_weight < weight )
        {
            max_weight = weight;
            p->direction_id = beam_id;
        }
    }

    if (max_weight > 0)
    {
        for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
        {
            p->beam_weight[beam_id] /= norm_weight;
            if (p->beam_weight[beam_id] < p->recur_beam_weight[beam_id])
            {
                p->recur_beam_weight[beam_id] = p->cfg.alpha_weight_ar * p->recur_beam_weight[beam_id] +
                                                alpha_weight_ar_1 * p->beam_weight[beam_id];
            }
            else
            {
                p->recur_beam_weight[beam_id] = p->cfg.alpha_weight_at * p->recur_beam_weight[beam_id] +
                                                alpha_weight_at_1 * p->beam_weight[beam_id];
            }
        }
    }
    else
    {
        for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
        {
            p->beam_weight[beam_id] = 0.125f;
            p->recur_beam_weight[beam_id] = p->cfg.alpha_weight_at * p->recur_beam_weight[beam_id] +
                                            alpha_weight_at_1 * p->beam_weight[beam_id];
        }
    }

    if (!p->complex_weighting)
    {
        int id_re_ref, id_im_ref;
        float amp, cos_val, sin_val;
        for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
        {
            id_re_ref = band_id + band_id;
            id_im_ref = id_re_ref + 1;
            sp_re = ref_sp[id_re_ref];
            sp_im = ref_sp[id_im_ref];
            amp = sqrtf(sp_re * sp_re + sp_im * sp_im + 1e-30f);
            cos_val = sp_re / amp;
            sin_val = sp_im / amp;
            for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            {
                weight = p->recur_beam_weight[beam_id];
                inst_psd_beam = p->inst_auto_psd + beam_id * psd_offset;
                if ( weight > 0 )
                {
                    amp = sqrtf(inst_psd_beam[band_id]) * weight;
                    output_sp[id_re_ref] += amp * cos_val;
                    output_sp[id_im_ref] += amp * sin_val;
                }
            }
        }
    }
    else
    {
        bf_sp_beam = input_sp;
        for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
        {
            weight = p->recur_beam_weight[beam_id];
            if ( weight > 0 )
            {
                for (int band_id = 0; band_id < AWI_FRAME_BAND_COUNT; band_id++)
                {
                    id_re = band_id + band_id;
                    id_im = id_re + 1;
                    output_sp[id_re] += weight * bf_sp_beam[id_re];
                    output_sp[id_im] += weight * bf_sp_beam[id_im];
                }
            }
            bf_sp_beam += sp_offset;
        }
    }

}

