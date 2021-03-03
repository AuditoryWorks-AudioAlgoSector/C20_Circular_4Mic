#include "awi_dft_synthesis_filter_bank.h"


void awi_initialize_dft_synthesis_filter_bank(dft_synthesis_filter_bank *p, float *buffer, int channelCounts, int bandCounts)
{
    p->channelCounts = channelCounts;
    p->bufferLen = AWI_SYN_FILTER_LENGTH;
    p->overlaySampleSize = AWI_SYN_FILTER_LENGTH - AWI_DECI_RATE;

    memset(buffer, 0, p->channelCounts * p->bufferLen * sizeof(float));
    memset(p->fft, 0, 2 * AWI_BAND_COUNT * sizeof(float));
    memset(p->newFrame, 0, AWI_SYN_FILTER_LENGTH * sizeof(float));
}


void awi_process_dft_synthesis_filter_bank(dft_synthesis_filter_bank *p, float *buffer, float *fbFrameIn, float *audioFrameOut)
{
    float *arr;

#ifdef AWI_ARM
    float factor = AWI_BAND_COUNT * AWI_BAND_COUNT / 2;
#else
    float factor = AWI_BAND_COUNT / 2;
#endif

    for (int channel = 0; channel < p->channelCounts; ++channel)
    {
        float *outPtr = audioFrameOut + channel * AWI_FRAME_LENGTH;

        for (int sample = 0; sample < AWI_FRAME_PER_BLOACK; ++sample)
        {
            arr = fbFrameIn + sample * p->channelCounts * AWI_FRAME_BAND_COUNT * 2 + channel * AWI_FRAME_BAND_COUNT * 2;
            memcpy(p->fft, arr, sizeof(float) * AWI_FRAME_BAND_COUNT * 2);

            for (int band = AWI_FRAME_BAND_COUNT; band < AWI_BAND_COUNT; ++band)
            {
                int index = AWI_BAND_COUNT - band;
                p->fft[2 * band] = p->fft[2 * index];
                p->fft[2 * band + 1] = - p->fft[2 * index + 1];
            }


#ifdef AWI_ARM
            awi_ifft_256(p->fft);
#else
            awi_ifft_1d(p->fft, p->fft, AWI_BAND_COUNT);
                        
#endif
            for (int m = 0; m < AWI_SYN_FILTER_LENGTH; ++m)
            {
                int fftIndex = m % AWI_BAND_COUNT;
                p->newFrame[m] = p->fft[2 * fftIndex] * AWI_SYN_FILTER[m];
            }

            float *oldPtr1 = &buffer[channel * p->bufferLen];
            float *oldPtr2 = oldPtr1 + AWI_DECI_RATE;
            float *newPtr = &p->newFrame[0];
            for (int k = 0; k < p->overlaySampleSize; ++k)
                *oldPtr1++ = (*newPtr++) + (*oldPtr2++);

            oldPtr2 = &buffer[channel * p->bufferLen];

            for (int k = 0; k < AWI_DECI_RATE; ++k)
            {
                *oldPtr1++ = *newPtr++;
                *outPtr++ = (factor * (*oldPtr2++));
            }

        }
    }
}
