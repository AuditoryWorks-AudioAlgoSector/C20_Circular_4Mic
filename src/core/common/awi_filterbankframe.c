#include "awi_filterbankframe.h"

// filter bank frame is [sample][channel][band][2], which is a complex number

float awi_get_channel_energy(float *p,  int channel, int totalChannel, int frameBandCounts, int framesPerBlock)
{
    int offset = (framesPerBlock - 1) * totalChannel * frameBandCounts * 2;
    float sum = 0;

    for (int band = 1; band < frameBandCounts - 1; ++band)
    {
        float real = p[offset + channel * frameBandCounts * 2 + band * 2];
        float imag = p[offset + channel * frameBandCounts * 2 + band * 2 + 1];
        sum += real * real;
        sum += imag * imag;
    }

    return sum;
}

float awi_get_energy(float *p, int sample, int channel, int channelCounts, int frameBandCounts)
{
    int offset = sample * frameBandCounts * channelCounts * 2 + channel * frameBandCounts * 2;
    float re;
    float im;
    float sum = 0;

    for(int i = 1; i < frameBandCounts - 1; i++)
    {
        re = p[offset + i * 2];
        im = p[offset + i * 2 + 1];
        sum += re * re + im * im;
    }

    return sum;
}

float *awi_get_band_pointer(float *p, int sample, int channel, int channelCounts, int frameBandCounts)
{
    return p + sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2;
}

float *awi_get_pointer(float *p)
{
    return p;
}

float *awi_get_frame_pointer(float *p, int sample, int channelCounts, int frameBandCounts)
{
    return p + sample * channelCounts * frameBandCounts * 2;
}

float awi_get_real_value(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts)
{
    return p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2];
}

float awi_get_imag_value(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts)
{
    return p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2 + 1];
}

float awi_get_band_power(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts)
{
    float real = p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2];
    float imag = p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2 + 1];
    return real * real + imag * imag;
}

void awi_set_value(float *p, int sample, int channel, int band, float real, float imag, int channelCounts, int frameBandCounts)
{
    p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2] = real;
    p[sample * channelCounts * frameBandCounts * 2 + channel * frameBandCounts * 2 + band * 2 + 1] = imag;
}
