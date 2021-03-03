#include "awi_channel_convert.h"

void awi_split_channel(float *dst, float *src, int frameSize, int channels)
{
    for (int i = 0; i < frameSize; ++i)
        for (int j = 0; j < channels; ++j)
            dst[j * frameSize + i] = src[i * channels + j];
}

void awi_inplace_split_channel(float *arr, int frameSize, int channels)
{
    float *tmp = (float *) malloc(sizeof(float) * frameSize * channels);
    awi_split_channel(tmp, arr, frameSize, channels);
    memcpy(arr, tmp, sizeof(float) * frameSize * channels);
    free(tmp);
}

void awi_merge_channel(float *dst, float *src, int frameSize, int channels)
{
    for (int i = 0; i < frameSize; ++i)
        for (int j = 0; j < channels; ++j)
            dst[i * channels + j] = src[j * frameSize + i];
}

void awi_inplace_merge_channel(float *arr, int frameSize, int channels)
{
    float *tmp = (float *) malloc(sizeof(float) * frameSize * channels);
    awi_merge_channel(tmp, arr, frameSize, channels);
    memcpy(arr, tmp, sizeof(float) * frameSize * channels);
    free(tmp);
}

void awi_normal_split_channel(float *dst, short *src, int frameSize, int channels)
{
    float factor = 1.f / 32768;
    for (int i = 0; i < frameSize; ++i)
        for (int j = 0; j < channels; ++j)
            dst[j * frameSize + i] = factor * src[i * channels + j];
}


void awi_normal_24bit_split_channel(float *dst, int32_t *src, int frameSize, int channels)
{
    float factor = 1.f / 2147483648;
    for (int i = 0; i < frameSize; ++i) {
        for (int j = 0; j < channels; ++j) {
            dst[j * frameSize + i] = factor * ( src[i * channels + j] & -256 );
        }
    }
}

