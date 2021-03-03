#include "awi_fbf.h"
#include <string.h>


void awi_fbf_init(awi_fbf_t *fbf)
{
    return;
}

void awi_fbf_process_single(awi_fbf_t *fbf, float *frame_in, float *frame_out, int beam)
{
    const int sp_offset = 2 * AWI_FRAME_BAND_COUNT;
    float *params = fbf->cfg.mvdr_fbf_weight + beam * AWI_FBF_MIC_NUM * sp_offset;
    float *frame0, *frame1, *frame2, *frame3;
    float *out1 = frame_out;
    int id_params0, id_params1, id_params2, id_params3;
    int id_in, id_in_plus1;
    float a0, b0, c0, d0, sumR0, sumI0;
    float a1, b1, c1, d1, sumR1, sumI1;
    float a2, b2, c2, d2, sumR2, sumI2;
    float a3, b3, c3, d3, sumR3, sumI3;

    frame0 = frame_in;
    frame1 = frame0 + sp_offset;
    frame2 = frame1 + sp_offset;
    frame3 = frame2 + sp_offset;


    for(int band = 0; band < AWI_FRAME_BAND_COUNT; ++band)
    {
        id_params0 = band + band;
        id_in     = id_params0;
        id_in_plus1 = id_in + 1;

        a0 = params[id_params0];
        b0 = params[id_params0 + 1];
        c0 = frame0[id_in];
        d0 = frame0[id_in_plus1];


        sumR0 = (a0 * c0 - b0 * d0);
        sumI0 = (a0 * d0 + b0 * c0);

        id_params1 = id_params0 + sp_offset;

        a1 = params[id_params1];
        b1 = params[id_params1 + 1];
        c1 = frame1[id_in];
        d1 = frame1[id_in_plus1];

        sumR1 = (a1 * c1 - b1 * d1);
        sumI1 = (a1 * d1 + b1 * c1);


        id_params2 = id_params1 + sp_offset;

        a2 = params[id_params2];
        b2 = params[id_params2 + 1];
        c2 = frame2[id_in];
        d2 = frame2[id_in_plus1];

        sumR2 = (a2 * c2 - b2 * d2);
        sumI2 = (a2 * d2 + b2 * c2);


        id_params3 = id_params2 + sp_offset;

        a3 = params[id_params3];
        b3 = params[id_params3 + 1];
        c3 = frame3[id_in];
        d3 = frame3[id_in_plus1];

        sumR3 = (a3 * c3 - b3 * d3);
        sumI3 = (a3 * d3 + b3 * c3);

        *out1++ = (sumR0 + sumR1 + sumR2 + sumR3);
        *out1++ = (sumI0 + sumI1 + sumI2 + sumI3);

    }

}

