#ifndef AUDIO_ALGORITHMS_DFTANALYSISFILTERBANK_H
#define AUDIO_ALGORITHMS_DFTANALYSISFILTERBANK_H

#include "awi_audioframe.h"
#include "awi_filterbankframe.h"
#include "awi_fft.h"

#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
typedef struct DftAnalysisFilterBank_
{
    int channelCounts;
    int bufferLen;
    int filterLen2bandRatio;
    float windowedValue[AWI_ANA_FILTER_LENGTH];
    float fftValue[AWI_BAND_COUNT * 2];
} dft_analysis_filter_bank;

typedef float ana_mic_buffer[AWI_MIC_CHANNEL * (AWI_FB_FRAME_LENGTH + AWI_ANA_FILTER_LENGTH)];
typedef float ana_spk_buffer[AWI_SPK_CHANNEL * (AWI_FB_FRAME_LENGTH + AWI_ANA_FILTER_LENGTH)];

void awi_initialize_dft_analysis_filter_bank(dft_analysis_filter_bank *p, float *buffer, int numChannels, int numBands);
void awi_free_dft_analysis_filter_bank(dft_analysis_filter_bank *p);
void awi_process_dft_analysis_filter_bank(dft_analysis_filter_bank *p, float *buffer, float *audioFrameIn, float *fbFrameOut);

#ifdef __cplusplus
};
#endif

#endif /* AUDIO_ALGORITHMS_DFTANALYSISFILTERBANK_H */
