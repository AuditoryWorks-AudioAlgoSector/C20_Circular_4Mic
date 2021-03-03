#include "awi_dft_analysis_filter_bank.h"

void awi_initialize_dft_analysis_filter_bank(dft_analysis_filter_bank *p, float *buffer, int numChannels, int numBands)
{
    p->channelCounts = numChannels;
    p->bufferLen = AWI_FB_FRAME_LENGTH + AWI_ANA_FILTER_LENGTH;
    p->filterLen2bandRatio = (int) (AWI_ANA_FILTER_LENGTH / AWI_BAND_COUNT);

    memset(buffer, 0, p->channelCounts * p->bufferLen * sizeof(float));
    memset(p->windowedValue, 0, AWI_ANA_FILTER_LENGTH * sizeof(float));
    memset(p->fftValue, 0, AWI_BAND_COUNT * 2 * sizeof(float));
}

void awi_process_dft_analysis_filter_bank(dft_analysis_filter_bank *p, float *buffer, float *audioFrameIn, float *fbFrameOut)
{
    float *arr;
    float *value;

    for (int channel = 0; channel < p->channelCounts; ++channel)
    {
        memmove(&buffer[channel * p->bufferLen], &buffer[channel * p->bufferLen + AWI_FB_FRAME_LENGTH], AWI_ANA_FILTER_LENGTH * sizeof(float));
        memcpy(&buffer[channel * p->bufferLen + AWI_ANA_FILTER_LENGTH - 1], audioFrameIn + channel * AWI_FB_FRAME_LENGTH, AWI_FB_FRAME_LENGTH * sizeof(float));
        for (int pos = 0; pos < AWI_FRAME_PER_BLOACK; ++pos)
        {
            value = &buffer[channel * p->bufferLen + pos * AWI_DECI_RATE];
            float *p_fftValue = p->fftValue;
            float *pSrcA = value;
            float *pSrcB = AWI_ANA_FILTER;
            float acc, bcc;

            for(int b = 0; b != AWI_BAND_COUNT; ++b)
            {
                acc = *(pSrcA + b);
                acc *= *(pSrcB + b);
                bcc = *(pSrcA + AWI_BAND_COUNT + b);
                bcc *= *(pSrcB + AWI_BAND_COUNT + b);
                *p_fftValue++ = acc + bcc;
                *p_fftValue++ = 0;
            }

#ifdef AWI_ARM
            awi_fft_256(p->fftValue);
#else
            awi_fft_1d(p->fftValue, p->fftValue, 256);
#endif
            arr = fbFrameOut + pos * p->channelCounts * AWI_FRAME_BAND_COUNT * 2 + channel * AWI_FRAME_BAND_COUNT * 2;
            memcpy(arr, p->fftValue, sizeof(float) * 2 * AWI_FRAME_BAND_COUNT);
            arr[1] = 0;
        }
    }
}
