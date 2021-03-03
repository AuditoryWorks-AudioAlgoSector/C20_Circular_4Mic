#ifndef _FFT_H_
#define _FFT_H_


#ifdef AWI_ARM

void awi_fft_512(float *arr);
void awi_ifft_512(float *arr);
void awi_fft_256(float *arr);
void awi_ifft_256(float *arr);

#else

#include "kiss_fft.h"
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif

void awi_float2cpx(kiss_fft_cpx *dst, float *src, const int length);
void awi_cpx2float(float *dst, kiss_fft_cpx *src, const int length);
void awi_fft_1d(float *dst, float *src, int length);
void awi_ifft_1d(float *dst, float *src, int length);

#endif

#ifdef __cplusplus
};
#endif
#endif // _FFT_H_
