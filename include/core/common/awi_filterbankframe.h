#ifndef AUDIO_ALGORITHMS_FILTERBANKFRAME_H
#define AUDIO_ALGORITHMS_FILTERBANKFRAME_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "awi_constant.h"
#include "awi_filterbankparams.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef float filter_mic_bank_frame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_MIC_CHANNEL * 2];
typedef float filter_spk_bank_frame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_SPK_CHANNEL * 2];

float awi_get_channel_energy(float *p,  int channel, int totalChannel, int frameBandCounts, int framesPerBlock);

float awi_get_energy(float *p, int sample, int channel, int channelCounts, int frameBandCounts);

float *awi_get_band_pointer(float *p, int sample, int channel, int channelCounts, int frameBandCounts);

float *awi_get_pointer(float *p);

float *awi_get_frame_pointer(float *p, int sample, int channelCounts, int frameBandCounts);

float awi_get_real_value(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts);

float awi_get_imag_value(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts);

float awi_get_band_power(float *p, int sample, int channel, int band, int channelCounts, int frameBandCounts);

void awi_set_value(float *p, int sample, int channel, int band, float real, float imag, int channelCounts, int frameBandCounts);

#ifdef __cplusplus
};
#endif

#endif
