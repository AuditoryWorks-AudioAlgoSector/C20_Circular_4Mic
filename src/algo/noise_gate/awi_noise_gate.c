//
// Created by xzc on 19-7-13.
//

#include "awi_noise_gate.h"
#include "awi_filterbankparams.h"
#include <math.h>

void awi_noise_gate_init(awi_noise_gate_t *p)
{
    p->recur_gain           = 1.f;
    p->recur_gain_prev      = 1.f;
    p->attack_hold_counter  = 0;
    p->release_hold_counter = 0;
}

void awi_noise_gate_process(awi_noise_gate_t *p, float *audio)
{
    int i;
    float audio_tmp, inst_gain;

    for(i = 0; i < AWI_FB_FRAME_LENGTH; i++)
    {

        audio_tmp = fabsf(audio[i]);
        if ( audio_tmp > p->cfg.threshold )
            inst_gain = 1.f;
        else
            inst_gain = p->cfg.scale_down;

      if ( p->recur_gain_prev > inst_gain )
      {
          p->attack_hold_counter++;
          if ( p->attack_hold_counter >= p->cfg.attack_hold_samples )
          {
              p->recur_gain = p->cfg.alpha_at * p->recur_gain_prev + ( 1 - p->cfg.alpha_at ) * inst_gain;
              p->release_hold_counter = 0;
          }
          else
              p->recur_gain = p->recur_gain_prev;
      }
      else
      {
          p->release_hold_counter++;
          if ( p->release_hold_counter >= p->cfg.release_hold_samples )
          {
              p->recur_gain = p->cfg.alpha_ar * p->recur_gain_prev + ( 1 - p->cfg.alpha_ar ) * inst_gain;
              p->attack_hold_counter = 0;
          }
          else
              p->recur_gain = p->recur_gain_prev;
      }


        p->recur_gain_prev = p->recur_gain;

        audio[i] = audio[i] * p->recur_gain;
    }
}

