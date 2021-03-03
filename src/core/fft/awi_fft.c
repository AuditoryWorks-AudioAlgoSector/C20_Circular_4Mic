#include "awi_fft.h"

void awi_float2cpx(kiss_fft_cpx *dst, float *src, const int length)
{
    for(int i = 0; i < length; i++)
    {
        dst[i].r = src[2 * i];
        dst[i].i = src[2 * i + 1];
    }
}

void awi_cpx2float(float *dst, kiss_fft_cpx *src, const int length)
{
    for(int i = 0; i < length; i++)
    {
        dst[2 * i] = src[i].r;
        dst[2 * i + 1] = src[i].i;
    }
}

void awi_fft_1d(float *dst, float *src, int length)
{
    // memory
    kiss_fft_cpx *in = (kiss_fft_cpx *) malloc(sizeof(kiss_fft_cpx) * length);

    awi_float2cpx(in, src, length);

    // inplace fft
    kiss_fft_cfg cfg = kiss_fft_alloc(length, 0, NULL, NULL);
    if( cfg == 0)
    {
        kiss_fft_free(cfg);
        free(in);
        return;
    }
    kiss_fft(cfg, in, in);

    awi_cpx2float(dst, in, length);

    // clear workspace
    kiss_fft_free(cfg);
    free(in);
}

void awi_ifft_1d(float *dst, float *src, int length)
{
    // memory
    kiss_fft_cpx *in = (kiss_fft_cpx *) malloc(sizeof(kiss_fft_cpx) * length);

    awi_float2cpx(in, src, length);

    // inplace fft
    kiss_fft_cfg cfg = kiss_fft_alloc(length, 1, NULL, NULL);
    if( cfg == 0)
    {
        kiss_fft_free(cfg);
        free(in);
        return;
    }
    kiss_fft(cfg, in, in);

    awi_cpx2float(dst, in, length);

    // clear workspace
    kiss_fft_free(cfg);
    free(in);
}
