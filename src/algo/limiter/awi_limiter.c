//
// Created by xzc on 19-7-1.
//

#include "awi_limiter.h"
#include "awi_filterbankparams.h"
#include <math.h>

void awi_limiter_init(awi_limiter_t *p)
{
    p->at_hold_counter      = 0;
    p->ar_hold_counter      = 0;
}


void awi_limiter_process(awi_limiter_t *p, float *audio, int is_speech)
{
    int i;
    float abs_audio, inst_alpha, inst_gain;
    float alpha_amp_inc_1 = 1 - p->cfg.alpha_amp_inc;
    float excess_gain;
    float max_amp = -100;
    float recur_amplifier_cur;

    for(i = 0; i < AWI_FB_FRAME_LENGTH; i++)
    {
        abs_audio = fabsf(audio[i]);
        max_amp = max_amp > abs_audio ? max_amp : abs_audio;
    }
    recur_amplifier_cur = p->cfg.limiter_am / (max_amp + AWI_EPS);

    p->cfg.recur_amplifier = p->cfg.recur_amplifier < recur_amplifier_cur ? p->cfg.recur_amplifier : recur_amplifier_cur;

    for(i = 0; i < AWI_FB_FRAME_LENGTH; i++)
    {
        /* expander */
        audio[i] *= p->cfg.recur_amplifier;
        abs_audio = fabsf(audio[i]);
        if (abs_audio < p->cfg.desired_am)
        {
            p->at_hold_counter++;
            if (p->at_hold_counter > p->cfg.at_hold_th && is_speech)
            {
                excess_gain = p->cfg.recur_amplifier + p->cfg.amp_step_inc;

                p->cfg.recur_amplifier = p->cfg.alpha_amp_inc * p->cfg.recur_amplifier +
                                          alpha_amp_inc_1 * excess_gain;
                p->ar_hold_counter = 0;
            }
        }
        else
        {
            p->ar_hold_counter++;

            if (p->ar_hold_counter > p->cfg.ar_hold_th)
            {
                p->cfg.recur_amplifier *= p->cfg.alpha_amp_dec;
                p->at_hold_counter = 0;
            }
        }

        if ( p->cfg.peak_am > p->cfg.limiter_am * 0.40f  && is_speech)
        {
            p->cfg.max_amplifier *= p->cfg.max_amp_dec;
            p->cfg.max_amplifier = p->cfg.max_amplifier > p->cfg.amp_low_th ?
                    p->cfg.max_amplifier : p->cfg.amp_low_th;
        }
        else
        {
            p->cfg.max_amplifier *= p->cfg.max_amp_inc;
            p->cfg.max_amplifier = p->cfg.max_amplifier > p->cfg.amp_high_th ?
                    p->cfg.amp_high_th : p->cfg.max_amplifier;
        }

        p->cfg.recur_amplifier = p->cfg.recur_amplifier > p->cfg.max_amplifier ?
                                  p->cfg.max_amplifier : p->cfg.recur_amplifier;

        /* limiter */
        if ( abs_audio > p->cfg.peak_am )
            inst_alpha = p->cfg.alpha_at;
        else
            inst_alpha = p->cfg.alpha_ar;

        p->cfg.peak_am = inst_alpha * p->cfg.peak_am + ( 1 - inst_alpha ) * abs_audio;

        inst_gain = p->cfg.limiter_am / p->cfg.peak_am;
        if( inst_gain >= 1 )
            inst_gain = 1;

        if( inst_gain > p->cfg.gain )
            inst_alpha = p->cfg.alpha_ar;
        else
            inst_alpha = p->cfg.alpha_at;

        p->cfg.gain = inst_alpha * p->cfg.gain + ( 1 - inst_alpha ) * inst_gain;

        if ( abs_audio > p->cfg.limiter_am )
            audio[i] = p->cfg.limiter_am;

        audio[i] *= p->cfg.gain;

    }
}

