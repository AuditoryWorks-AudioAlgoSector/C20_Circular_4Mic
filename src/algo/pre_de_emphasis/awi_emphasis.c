//
// Created by xzc on 19-7-26.
//

#include "awi_emphasis.h"
#include "awi_audioframe.h"

void awi_emphasis_init(awi_emphasis_t *p)
{
    memset(p->pre_emp_first_in_chs, 0, AWI_MIC_CHANNEL * sizeof(float));
    p->de_emp_first_out = 0.f;
    p->pre_emp_first_in_spk = 0.f;
    p->de_emp_mem = 0.f;
}


void awi_pre_emphasis_process(awi_emphasis_t *p, float *audio, float *first_in)
{

    float last_val;
    last_val = audio[AWI_FRAME_LENGTH-1];
    for (int i = AWI_FRAME_LENGTH-1; i > 0; i--)
    {
        audio[i] = audio[i] - p->cfg.alpha_pre * audio[i-1];

        if (audio[i] > 1)
            audio[i] = 0.9999f;

        if (audio[i] < -1)
            audio[i] = -0.9999f;
    }
    audio[0] = audio[0] - p->cfg.alpha_pre * ( *first_in );
    *first_in = last_val;

}

void awi_de_emphasis_process(awi_emphasis_t *p, float *audio)
{
    float tmp_result;
    for (int i = 0; i < AWI_FRAME_LENGTH; i++)
    {
        tmp_result = audio[i] + p->cfg.alpha_de * p->de_emp_mem;
        p->de_emp_mem = tmp_result;
        audio[i] = tmp_result;
    }
}

