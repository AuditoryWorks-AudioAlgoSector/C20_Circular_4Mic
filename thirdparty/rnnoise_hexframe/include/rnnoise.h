/* Copyright (c) 2018 Gregor Richards
 * Copyright (c) 2017 Mozilla */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef RNNOISE_H
#define RNNOISE_H 1

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef RNNOISE_BUILD
#define RNNOISE_BUILD
#endif // !RNNOISE_BUILD

# if defined(_WIN32)
#define DLL_EXPORT
#endif

#ifndef RNNOISE_EXPORT
# if defined(_WIN32)
#  if defined(RNNOISE_BUILD) && defined(DLL_EXPORT)
#   define RNNOISE_EXPORT __declspec(dllexport)
#  else
#   define RNNOISE_EXPORT
#  endif
# elif defined(__GNUC__) && defined(RNNOISE_BUILD)
#  define RNNOISE_EXPORT __attribute__ ((visibility ("default")))
# else
#  define RNNOISE_EXPORT
# endif
#endif

#define HEXMS_FRAME 1

#ifdef BAND_EXPAND
#define FRAME_SIZE_SHIFT 1
#define NB_BANDS 35
#define CEPS_MEM 8
#define NB_DELTA_CEPS 12
#elif HEXMS_FRAME
#define FRAME_SIZE_SHIFT 2
#define NB_BANDS 22
#define CEPS_MEM 8
#define NB_DELTA_CEPS 6
#else
#define FRAME_SIZE_SHIFT 2
#define NB_BANDS 18
#define CEPS_MEM 8
#define NB_DELTA_CEPS 6
#endif

#define BLOCK_SIZE 8000

#ifdef HEXMS_FRAME
#define FRAME_SIZE (256) // 256 means 16ms for 16k smaple rate
#define PITCH_MIN_PERIOD 32
#define PITCH_MAX_PERIOD 416
#define PITCH_FRAME_SIZE 512 //how to calculate pitch?
#define FREQ_DELTA (31.25)
#else
#define FRAME_SIZE (160) // 160 means 10ms for 16k sample rate
#define PITCH_MIN_PERIOD 20
#define PITCH_MAX_PERIOD 256
#define PITCH_FRAME_SIZE 320
#endif

#define WINDOW_SIZE (2*FRAME_SIZE)
#define FREQ_SIZE (FRAME_SIZE + 1)

#define PITCH_BUF_SIZE (PITCH_MAX_PERIOD+PITCH_FRAME_SIZE)
#define SQUARE(x) ((x)*(x))

#define NB_FEATURES (NB_BANDS+3*NB_DELTA_CEPS+2)

#ifndef TRAINING
#define TRAINING 1
#endif

typedef struct DenoiseState DenoiseState;
typedef struct RNNModel RNNModel;
typedef struct FeatureExtractor FeatureExtractor;

RNNOISE_EXPORT int rnnoise_get_size();

RNNOISE_EXPORT int rnnoise_init(DenoiseState *st, RNNModel *model);

RNNOISE_EXPORT DenoiseState *rnnoise_create(RNNModel *model);

RNNOISE_EXPORT void rnnoise_destroy(DenoiseState *st);

RNNOISE_EXPORT float rnnoise_process_frame(DenoiseState *st, float *out, const float *in);

RNNOISE_EXPORT float rnnoise_process_frame_vad(DenoiseState *st, const float *in);

RNNOISE_EXPORT RNNModel *rnnoise_model_from_file(FILE *f);

RNNOISE_EXPORT void rnnoise_model_free(RNNModel *model);

//feature extraction

RNNOISE_EXPORT float rnnoise_get_speech_gain(FeatureExtractor* st);

RNNOISE_EXPORT float rnnoise_get_noise_gain(FeatureExtractor* st);

RNNOISE_EXPORT float* rnnoise_get_Ln(FeatureExtractor* st);

RNNOISE_EXPORT float* rnnoise_get_features(FeatureExtractor* st);

RNNOISE_EXPORT float* rnnoise_get_gain(FeatureExtractor* st);

RNNOISE_EXPORT float* rnnoise_get_vad(FeatureExtractor* st);

RNNOISE_EXPORT int rnnoise_feature_extractor_init(FeatureExtractor* st);

RNNOISE_EXPORT FeatureExtractor* rnnoise_feature_extractor_create();

RNNOISE_EXPORT void rnnoise_feature_extractor_random_change(FeatureExtractor* st);

RNNOISE_EXPORT void rnnoise_feature_extract(FeatureExtractor* st, short* xi, short* ni);

RNNOISE_EXPORT void rnnoise_set_speech_energy(FeatureExtractor *st, float val);

RNNOISE_EXPORT void rnnoise_set_noise_energy(FeatureExtractor *st, float val);

RNNOISE_EXPORT void rnnoise_set_noise_factor(FeatureExtractor *st);

#ifdef __cplusplus
}
#endif

#endif
