#ifndef MFCC_H
#define MFCC_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fftw3.h"

#ifdef __cplusplus
extern "C"{
#endif


void mfcc(float *mfcc_result,float *powspectrum,float *data,int Sr,int n_filt, int n_cep, int n_fft, int cep_lifter, int append_Energy);

void get_filterbank_parameters(float **fbank, int n_filt, int sampleRate, int n_fft);

float hztomel(float hz);

float meltohz(float mel);



#ifdef __cplusplus
};
#endif


#endif//MFCC_H