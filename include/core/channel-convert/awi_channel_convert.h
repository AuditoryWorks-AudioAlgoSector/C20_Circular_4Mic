#ifndef AWI_CHANNEL_CONVERT_H
#define AWI_CHANNEL_CONVERT_H

#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void awi_split_channel(float *dst, float *src, int frameSize, int channels);

void awi_inplace_split_channel(float *arr, int frameSize, int channels);

void awi_merge_channel(float *dst, float *src, int frameSize, int channels);

void awi_inplace_merge_channel(float *arr, int frameSize, int channels);

void awi_normal_split_channel(float *dst, short *src, int frameSize, int channels);

void awi_normal_24bit_split_channel(float *dst, int32_t *src, int frameSize, int channels);

#ifdef __cplusplus
};
#endif

#endif // AWI_CHANNEL_CONVERT_H
