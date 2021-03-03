#ifndef AUDIO_ALGORITHMS_DFTSYNTHESISFILTERBANK_H
#define AUDIO_ALGORITHMS_DFTSYNTHESISFILTERBANK_H

#include "../common/awi_audioframe.h"
#include "../common/awi_filterbankframe.h"
#include "../common/awi_filterbankparams.h"
#include "../fft/awi_fft.h"
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct dft_synthesis_filter_bank_
{
    int channelCounts;
    int bufferLen;
    int overlaySampleSize;
    float fft[AWI_BAND_COUNT * 2];
    float newFrame[AWI_SYN_FILTER_LENGTH];
    
} dft_synthesis_filter_bank;

typedef float syn_mic_buffer[AWI_MIC_CHANNEL * AWI_SYN_FILTER_LENGTH];
typedef float syn_spk_buffer[AWI_SPK_CHANNEL * AWI_SYN_FILTER_LENGTH];

void awi_initialize_dft_synthesis_filter_bank(dft_synthesis_filter_bank *p, float *buffer,  int channelCounts, int bandCounts);
void awi_process_dft_synthesis_filter_bank(dft_synthesis_filter_bank *p, float *buffer,  float *fbFrameIn, float *audioFrameOut);
#ifdef __cplusplus
};
#endif

#endif
