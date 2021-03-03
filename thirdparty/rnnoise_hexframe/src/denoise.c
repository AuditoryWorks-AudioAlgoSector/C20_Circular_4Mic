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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "kiss_fft_opus.h"
#include "common.h"
#include <math.h>
#include "rnnoise.h"
#include "pitch.h"
#include "arch.h"
#include "rnn.h"
#include "rnn_data.h"
//#include "dr_wav.h"


/* The built-in model, used if no file is given as input */
extern const struct RNNModel rnnoise_model_orig;

#ifdef BAND_EXPAND
static const opus_int16 eband5ms[] = {
        /*0 100 200 300 400 500 600 700 800 900 1k 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2k 2.2 2.4 */
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24,
        /*  2.6 2.8 3.0 3.2 3.6 4k 4.4 4.8 5.2 5.6 6.2 6.8 7.4 8k 8.8 9.6 10.8 12k 13.8 15.6 17.8 20k*/
        26, 28, 30, 32, 36, 40, 44, 48, 52, 56, 62, 68, 74, 80, 88, 96, 108, 120, 138, 156, 178, 200
};
#elif HEXMS_FRAME
static const opus_int16 eband5ms[24] = {
        0, 100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270, 1480, 1720,
        2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000
};
#else
static const opus_int16 eband5ms[] = {
      /*0  200 400 600 800  1k 1.2 1.4 1.6  2k 2.4  2.8 3.2  4k 4.8 5.6 6.8  8k 9.6 12k 15.6 20k*/
        0, 1,  2,   3,  4,  5, 6,  7,   8,  10, 12, 14, 16, 20, 24, 28, 34, 40, 48, 60, 78, 100
};
#endif



typedef struct {
    int init;
    kiss_fft_state *kfft;
    float half_window[FRAME_SIZE];
    float dct_table[NB_BANDS * NB_BANDS];
} CommonState;

struct DenoiseState {
    float analysis_mem[FRAME_SIZE];
    float cepstral_mem[CEPS_MEM][NB_BANDS];
    int memid;
    float synthesis_mem[FRAME_SIZE];
    float pitch_buf[PITCH_BUF_SIZE];
    float pitch_enh_buf[PITCH_BUF_SIZE];
    float last_gain;
    int last_period;
    float mem_hp_x[2];
    float lastg[NB_BANDS];
    RNNState rnn;
};

CommonState common;

#if 1
int low_pass = FREQ_SIZE;
int band_lp = NB_BANDS;
#endif

short rnnoise_get_eband5ms(int i) {
    if (i >= 22) {
        return -1;
    }
    return eband5ms[i];
}

#ifdef HEXMS_FRAME
void compute_band_energy(float *bandE, const kiss_fft_cpx *X){
    float sum[NB_BANDS] = {0};
    for (int i = 0; i < NB_BANDS - 1; i++)
    {
        int freq_start = floor(eband5ms[i]/31.25) + 1;
        int freq_max = floor(eband5ms[i+1]/31.25);
        int band_size = freq_max - freq_start + 1;
        for(int j = freq_start; j <= freq_max; j++)
        {
            float tmp;
            float frac = (float)(j-freq_start) / band_size;
            tmp = SQUARE(X[j].r);
            tmp += SQUARE(X[j].i);
            sum[i] += (1 - frac) * tmp;
            sum[i + 1] += frac * tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS - 1] *= 2;
    for (int i = 0; i < NB_BANDS; i++) {
        bandE[i] = sum[i];
    }
}

void compute_band_corr(float *bandE, const kiss_fft_cpx *X, const kiss_fft_cpx *P){
    float sum[NB_BANDS] = {0};
    for (int i = 0; i < NB_BANDS - 1; i++)
    {
        int freq_start = floor(eband5ms[i]/31.25) + 1;
        int freq_max = floor(eband5ms[i+1]/31.25);
        int band_size = freq_max - freq_start + 1;
        for(int j = freq_start; j <= freq_max; j++)
        {
            float tmp;
            float frac = (float)(j-freq_start) / band_size;
            tmp = X[j].r * P[j].r;
            tmp += X[j].i * P[j].i;
            sum[i] += (1 - frac) * tmp;
            sum[i + 1] += frac * tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS - 1] *= 2;
    for (int i = 0; i < NB_BANDS; i++) {
        bandE[i] = sum[i];
    }
}

void interp_band_gain(float *g, const float *bandE){
    memset(g, 0, FREQ_SIZE);
    for (int i = 0; i < NB_BANDS - 1; i++) {
        int freq_start = floor(eband5ms[i]/31.25) + 1;
        int freq_max = floor(eband5ms[i+1]/31.25);
        int band_size = freq_max - freq_start + 1;
        for(int j = freq_start; j <= freq_max; j++) {
            float frac = (float) (j - freq_start) / band_size;
            g[j] = (1 - frac) * bandE[i] + frac * bandE[i + 1];
        }
    }
}
#else
void compute_band_energy(float *bandE, const kiss_fft_cpx *X) {
    int i;
    float sum[NB_BANDS] = {0};
    for (i = 0; i < NB_BANDS - 1; i++) {
        int j;
        int band_size;
        band_size = (eband5ms[i + 1] - eband5ms[i]) << FRAME_SIZE_SHIFT;
        for (j = 0; j < band_size; j++) {
            float tmp;
            float frac = (float) j / band_size;
            tmp = SQUARE(X[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].r);
            tmp += SQUARE(X[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].i);
            sum[i] += (1 - frac) * tmp;
            sum[i + 1] += frac * tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS - 1] *= 2;
    for (i = 0; i < NB_BANDS; i++) {
        bandE[i] = sum[i];
    }
}

void compute_band_corr(float *bandE, const kiss_fft_cpx *X, const kiss_fft_cpx *P) {
    int i;
    float sum[NB_BANDS] = {0};
    for (i = 0; i < NB_BANDS - 1; i++) {
        int j;
        int band_size;
        band_size = (eband5ms[i + 1] - eband5ms[i]) << FRAME_SIZE_SHIFT;
        for (j = 0; j < band_size; j++) {
            float tmp;
            float frac = (float) j / band_size;
            tmp = X[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].r * P[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].r;
            tmp += X[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].i * P[(eband5ms[i] << FRAME_SIZE_SHIFT) + j].i;
            sum[i] += (1 - frac) * tmp;
            sum[i + 1] += frac * tmp;
        }
    }
    sum[0] *= 2;
    sum[NB_BANDS - 1] *= 2;
    for (i = 0; i < NB_BANDS; i++) {
        bandE[i] = sum[i];
    }
}

void interp_band_gain(float *g, const float *bandE) {
    int i;
    memset(g, 0, FREQ_SIZE);
    for (i = 0; i < NB_BANDS - 1; i++) {
        int j;
        int band_size;
        band_size = (eband5ms[i + 1] - eband5ms[i]) << FRAME_SIZE_SHIFT;
        for (j = 0; j < band_size; j++) {
            float frac = (float) j / band_size;
            g[(eband5ms[i] << FRAME_SIZE_SHIFT) + j] = (1 - frac) * bandE[i] + frac * bandE[i + 1];
        }
    }
}
#endif

static void check_init() {
    int i;
    if (common.init) return;
    common.kfft = opus_fft_alloc_twiddles(2 * FRAME_SIZE, NULL, NULL, NULL, 0);
    for (i = 0; i < FRAME_SIZE; i++)
        common.half_window[i] = sin(
                .5 * M_PI * sin(.5 * M_PI * (i + .5) / FRAME_SIZE) * sin(.5 * M_PI * (i + .5) / FRAME_SIZE));
    for (i = 0; i < NB_BANDS; i++) {
        int j;
        for (j = 0; j < NB_BANDS; j++) {
            common.dct_table[i * NB_BANDS + j] = cos((i + .5) * j * M_PI / NB_BANDS);
            if (j == 0) common.dct_table[i * NB_BANDS + j] *= sqrt(.5);
        }
    }
    common.init = 1;
}

static void dct(float *out, const float *in) {
    int i;
    check_init();
    for (i = 0; i < NB_BANDS; i++) {
        int j;
        float sum = 0;
        for (j = 0; j < NB_BANDS; j++) {
            sum += in[j] * common.dct_table[j * NB_BANDS + i];
        }
        out[i] = sum * sqrt(2. / 22);
    }
}

#if 0
static void idct(float *out, const float *in) {
  int i;
  check_init();
  for (i=0;i<NB_BANDS;i++) {
    int j;
    float sum = 0;
    for (j=0;j<NB_BANDS;j++) {
      sum += in[j] * common.dct_table[i*NB_BANDS + j];
    }
    out[i] = sum*sqrt(2./22);
  }
}
#endif

static void forward_transform(kiss_fft_cpx *out, const float *in) {
    int i;
    kiss_fft_cpx x[WINDOW_SIZE];
    kiss_fft_cpx y[WINDOW_SIZE];
    check_init();
    for (i = 0; i < WINDOW_SIZE; i++) {
        x[i].r = in[i];
        x[i].i = 0;
    }
    opus_fft(common.kfft, x, y, 0);
    for (i = 0; i < FREQ_SIZE; i++) {
        out[i] = y[i];
    }
}

static void inverse_transform(float *out, const kiss_fft_cpx *in) {
    int i;
    kiss_fft_cpx x[WINDOW_SIZE];
    kiss_fft_cpx y[WINDOW_SIZE];
    check_init();
    for (i = 0; i < FREQ_SIZE; i++) {
        x[i] = in[i];
    }
    for (; i < WINDOW_SIZE; i++) {
        x[i].r = x[WINDOW_SIZE - i].r;
        x[i].i = -x[WINDOW_SIZE - i].i;
    }
    opus_fft(common.kfft, x, y, 0);
    /* output in reverse order for IFFT. */
    out[0] = WINDOW_SIZE * y[0].r;
    for (i = 1; i < WINDOW_SIZE; i++) {
        out[i] = WINDOW_SIZE * y[WINDOW_SIZE - i].r;
    }
}

static void apply_window(float *x) {
    int i;
    check_init();
    for (i = 0; i < FRAME_SIZE; i++) {
        x[i] *= common.half_window[i];
        x[WINDOW_SIZE - 1 - i] *= common.half_window[i];
    }
}

int rnnoise_get_size() {
    return sizeof(DenoiseState);
}

int rnnoise_init(DenoiseState *st, RNNModel *model) {
    memset(st, 0, sizeof(*st));
    if (model)
        st->rnn.model = model;
    else
        st->rnn.model = &rnnoise_model_orig;
    st->rnn.vad_gru_state = (float *) calloc(sizeof(float), st->rnn.model->vad_gru_size);
    st->rnn.noise_gru_state = (float *) calloc(sizeof(float), st->rnn.model->noise_gru_size);
    st->rnn.denoise_gru_state = (float *) calloc(sizeof(float), st->rnn.model->denoise_gru_size);
    return 0;
}

DenoiseState *rnnoise_create(RNNModel *model) {
    DenoiseState *st;
    st = (DenoiseState *) malloc(rnnoise_get_size());
    rnnoise_init(st, model);
    return st;
}

void rnnoise_destroy(DenoiseState *st) {
    free(st->rnn.vad_gru_state);
    free(st->rnn.noise_gru_state);
    free(st->rnn.denoise_gru_state);
    free(st);
}

// fft + band energy
static void frame_analysis(DenoiseState *st, kiss_fft_cpx *X, float *Ex, const float *in) {
    int i;
    float x[WINDOW_SIZE];
    RNN_COPY(x, st->analysis_mem, FRAME_SIZE);
    for (i = 0; i < FRAME_SIZE; i++) x[FRAME_SIZE + i] = in[i];
    RNN_COPY(st->analysis_mem, in, FRAME_SIZE);
    apply_window(x);
    forward_transform(X, x);
#if TRAINING
    for (i = low_pass; i < FREQ_SIZE; i++)
        X[i].r = X[i].i = 0;
#endif
    compute_band_energy(Ex, X);
}

static int compute_frame_features(DenoiseState *st, kiss_fft_cpx *X, kiss_fft_cpx *P,
                                  float *Ex, float *Ep, float *Exp, float *features, const float *in) {
    int i;
    float E = 0;
    float *ceps_0, *ceps_1, *ceps_2;
    float spec_variability = 0;
    float Ly[NB_BANDS];
    float p[WINDOW_SIZE];
    float pitch_buf[PITCH_BUF_SIZE >> 1];
    int pitch_index;
    float gain;
    float *(pre[1]);
    float tmp[NB_BANDS];
    float follow, logMax;

    // in -> X,Ex
    frame_analysis(st, X, Ex, in);

    // in -> pitch_index, p
    // pitch_buf save last PITCH_BUF_SIZE history data
    RNN_MOVE(st->pitch_buf, &st->pitch_buf[FRAME_SIZE], PITCH_BUF_SIZE - FRAME_SIZE);
    RNN_COPY(&st->pitch_buf[PITCH_BUF_SIZE - FRAME_SIZE], in, FRAME_SIZE);
    pre[0] = &st->pitch_buf[0];
    // down sample pitch _buf to half
    pitch_downsample(pre, pitch_buf, PITCH_BUF_SIZE, 1);
    // search the best pitch
    pitch_search(pitch_buf + (PITCH_MAX_PERIOD >> 1), pitch_buf, PITCH_FRAME_SIZE,
                 PITCH_MAX_PERIOD - 3 * PITCH_MIN_PERIOD, &pitch_index);
    pitch_index = PITCH_MAX_PERIOD - pitch_index;
    gain = remove_doubling(pitch_buf, PITCH_MAX_PERIOD, PITCH_MIN_PERIOD,
                           PITCH_FRAME_SIZE, &pitch_index, st->last_period, st->last_gain);
    st->last_period = pitch_index;
    st->last_gain = gain; // this gain is not used now
    for (i = 0; i < WINDOW_SIZE; i++)
        p[i] = st->pitch_buf[PITCH_BUF_SIZE - WINDOW_SIZE - pitch_index + i];

    // p -> P -> Ep
    apply_window(p);
    forward_transform(P, p);
    compute_band_energy(Ep, P);

    // X,Ex,P,Ep -> Exp
    compute_band_corr(Exp, X, P);
    for (i = 0; i < NB_BANDS; i++) Exp[i] = Exp[i] / sqrt(.001 + Ex[i] * Ep[i]);

    // dct(Exp) -> features[34:39]
    dct(tmp, Exp);
    for (i = 0; i < NB_DELTA_CEPS; i++) features[NB_BANDS + 2 * NB_DELTA_CEPS + i] = tmp[i];
    features[NB_BANDS + 2 * NB_DELTA_CEPS] -= 1.3;
    features[NB_BANDS + 2 * NB_DELTA_CEPS + 1] -= 0.9;

    // pitch_index -> features[40]
    features[NB_BANDS + 3 * NB_DELTA_CEPS] = .01 * (pitch_index - 300);

    // Ex -> Ly(log of band engergy), E(sum of band energy)
    // TODO: why limit the distnace-to-max and decrease-speed ?
    logMax = -2;
    follow = -2;
    for (i = 0; i < NB_BANDS; i++) {
        Ly[i] = log10(1e-2 + Ex[i]);
        Ly[i] = MAX16(logMax - 7, MAX16(follow - 1.5, Ly[i]));
        logMax = MAX16(logMax, Ly[i]);
        follow = MAX16(follow - 1.5, Ly[i]);
        E += Ex[i];
    }

    // when not training and energy low
    if ((!TRAINING) && E < 0.04) {
        /* If there's no audio, avoid messing up the state. */
        //fprintf(stderr, " If there's no audio, avoid messing up the state.\r");
        RNN_CLEAR(features, NB_FEATURES);
        return 1;
    }

    // Ly -> bfcc -> features[0:21]
    dct(features, Ly);
    features[0] -= 12;
    features[1] -= 4;

    // ceps pointer <- ceps ring buffer
    ceps_0 = st->cepstral_mem[st->memid];
    ceps_1 = (st->memid < 1) ? st->cepstral_mem[CEPS_MEM + st->memid - 1] : st->cepstral_mem[st->memid - 1];
    ceps_2 = (st->memid < 2) ? st->cepstral_mem[CEPS_MEM + st->memid - 2] : st->cepstral_mem[st->memid - 2];

    // features[0:21] -> ceps ring buffer
    for (i = 0; i < NB_BANDS; i++) ceps_0[i] = features[i];
    st->memid++;

    // ceps[0:5]*3 -> features[0:5]
    // d(ceps)     -> features[22:27]
    // d2(ceps)    -> features[28:33]
    for (i = 0; i < NB_DELTA_CEPS; i++) {
        features[i] = ceps_0[i] + ceps_1[i] + ceps_2[i];
        features[NB_BANDS + i] = ceps_0[i] - ceps_2[i];
        features[NB_BANDS + NB_DELTA_CEPS + i] = ceps_0[i] - 2 * ceps_1[i] + ceps_2[i];
    }

    /* Spectral variability features. */
    if (st->memid == CEPS_MEM) st->memid = 0; // this line should just after st->memid++;
    for (i = 0; i < CEPS_MEM; i++) {
        int j;
        float mindist = 1e15f; // min dist from i to others
        for (j = 0; j < CEPS_MEM; j++) {
            int k;
            float dist = 0; // sum of bands' dist^2
            for (k = 0; k < NB_BANDS; k++) {
                float tmp;
                tmp = st->cepstral_mem[i][k] - st->cepstral_mem[j][k];
                dist += tmp * tmp;
            }
            if (j != i)
                mindist = MIN32(mindist, dist);
        }
        spec_variability += mindist;
    }

    // spec_variability -> features[41]
    features[NB_BANDS + 3 * NB_DELTA_CEPS + 1] = spec_variability / CEPS_MEM - 2.1;

    // silence in training
    return TRAINING && E < 0.1;
}

static void frame_synthesis(DenoiseState *st, float *out, const kiss_fft_cpx *y) {
    float x[WINDOW_SIZE];
    int i;
    inverse_transform(x, y);
    apply_window(x);
    for (i = 0; i < FRAME_SIZE; i++) out[i] = x[i] + st->synthesis_mem[i];
    RNN_COPY(st->synthesis_mem, &x[FRAME_SIZE], FRAME_SIZE);
}

static void biquad(float *y, float mem[2], const float *x, const float *b, const float *a, int N) {
    int i;
    for (i = 0; i < N; i++) {
        float xi, yi;
        xi = x[i];
        yi = x[i] + mem[0];
        mem[0] = mem[1] + (b[0] * (double) xi - a[0] * (double) yi);
        mem[1] = (b[1] * (double) xi - a[1] * (double) yi);
        y[i] = yi;
    }
}

void rnnoise_biquad(float *y, float mem[2], const float *x, const float *b, const float *a, int N) {
    biquad(y, mem, x, b, a, N);
}

void pitch_filter(kiss_fft_cpx *X, const kiss_fft_cpx *P, const float *Ex, const float *Ep,
                  const float *Exp, const float *g) {
    int i;
    float r[NB_BANDS];
    float rf[FREQ_SIZE] = {0};
    for (i = 0; i < NB_BANDS; i++) {
#if 0
        if (Exp[i]>g[i]) r[i] = 1;
        else r[i] = Exp[i]*(1-g[i])/(.001 + g[i]*(1-Exp[i]));
        r[i] = MIN16(1, MAX16(0, r[i]));
#else
        if (Exp[i] > g[i]) r[i] = 1;
        else r[i] = SQUARE(Exp[i]) * (1 - SQUARE(g[i])) / (.001 + SQUARE(g[i]) * (1 - SQUARE(Exp[i])));
        r[i] = sqrt(MIN16(1, MAX16(0, r[i])));
#endif
        r[i] *= sqrt(Ex[i] / (1e-8 + Ep[i]));
    }
    interp_band_gain(rf, r);
    for (i = 0; i < FREQ_SIZE; i++) {
        X[i].r += rf[i] * P[i].r;
        X[i].i += rf[i] * P[i].i;
    }
    float newE[NB_BANDS];
    compute_band_energy(newE, X);
    float norm[NB_BANDS];
    float normf[FREQ_SIZE] = {0};
    for (i = 0; i < NB_BANDS; i++) {
        norm[i] = sqrt(Ex[i] / (1e-8 + newE[i]));
    }
    interp_band_gain(normf, norm);
    for (i = 0; i < FREQ_SIZE; i++) {
        X[i].r *= normf[i];
        X[i].i *= normf[i];
    }
}

float rnnoise_process_frame(DenoiseState *st, float *out, const float *in) {
    int i;
    kiss_fft_cpx X[FREQ_SIZE];
    kiss_fft_cpx P[WINDOW_SIZE];
    float x[FRAME_SIZE];
    float Ex[NB_BANDS], Ep[NB_BANDS];
    float Exp[NB_BANDS];
    float features[NB_FEATURES];
    float g[NB_BANDS];
    float gf[FREQ_SIZE] = {1};
    float vad_prob = 0;
    int silence;
    static const float a_hp[2] = {-1.99599, 0.99600};
    static const float b_hp[2] = {-2, 1};
    biquad(x, st->mem_hp_x, in, b_hp, a_hp, FRAME_SIZE);
    silence = compute_frame_features(st, X, P, Ex, Ep, Exp, features, x);

    if (!silence) {
        compute_rnn(&st->rnn, g, &vad_prob, features);
        pitch_filter(X, P, Ex, Ep, Exp, g);
        for (i = 0; i < NB_BANDS; i++) {
            float alpha = .6f;
            g[i] = MAX16(g[i], alpha * st->lastg[i]);
            st->lastg[i] = g[i];
        }
        interp_band_gain(gf, g);
#if 1
        for (i = 0; i < FREQ_SIZE; i++) {
            X[i].r *= gf[i];
            X[i].i *= gf[i];
        }
#endif
    }

    frame_synthesis(st, out, X);
    return vad_prob;
}

float rnnoise_process_frame_vad(DenoiseState *st, const float *in) {
    int i;
    kiss_fft_cpx X[FREQ_SIZE];
    kiss_fft_cpx P[WINDOW_SIZE];
    float x[FRAME_SIZE];
    float Ex[NB_BANDS], Ep[NB_BANDS];
    float Exp[NB_BANDS];
    float features[NB_FEATURES];
    float g[NB_BANDS];
    float gf[FREQ_SIZE] = {1};
    float vad_prob = 0;
    int silence;
    static const float a_hp[2] = {-1.99599, 0.99600};
    static const float b_hp[2] = {-2, 1};
    biquad(x, st->mem_hp_x, in, b_hp, a_hp, FRAME_SIZE);
    silence = compute_frame_features(st, X, P, Ex, Ep, Exp, features, x);

    if (!silence) {
        compute_rnn(&st->rnn, g, &vad_prob, features);
        pitch_filter(X, P, Ex, Ep, Exp, g);
        for (i = 0; i < NB_BANDS; i++) {
            float alpha = .6f;
            g[i] = MAX16(g[i], alpha * st->lastg[i]);
            st->lastg[i] = g[i];
        }
        interp_band_gain(gf, g);
    }
    return vad_prob;
}



static float uni_rand() {
    return rand() / (double) RAND_MAX - .5;
}

static void rand_resp(float *a, float *b) {
    a[0] = .75 * uni_rand();
    a[1] = .75 * uni_rand();
    b[0] = .75 * uni_rand();
    b[1] = .75 * uni_rand();
}

struct FeatureExtractor {
    float a_hp[2];
    float b_hp[2];
    float a_noise[2];
    float b_noise[2];
    float a_sig[2];
    float b_sig[2];
    float mem_hp_x[2];
    float mem_hp_n[2];
    float mem_resp_x[2];
    float mem_resp_n[2];
    float x[FRAME_SIZE];
    float n[FRAME_SIZE];
    float xn[FRAME_SIZE];
    int vad_cnt;
    //int gain_change_count;
    float speech_gain, noise_gain;
    DenoiseState *xst;
    DenoiseState *nst;
    DenoiseState *xnst;
    kiss_fft_cpx X[FREQ_SIZE], Y[FREQ_SIZE], N[FREQ_SIZE], P[WINDOW_SIZE];
    float Ex[NB_BANDS], Ey[NB_BANDS], En[NB_BANDS], Ep[NB_BANDS];
    float Exp[NB_BANDS];
    float Ln[NB_BANDS];
    float features[NB_FEATURES];
    float g[NB_BANDS];
    float vad;
    float E;
    float expected_snr;
    float speech_energy;
    float noise_energy;
    double noise_factor;
};

float rnnoise_get_speech_gain(FeatureExtractor *st) {
    return st->speech_gain;
}

float rnnoise_get_noise_gain(FeatureExtractor *st) {
    return st->noise_gain;
}

float *rnnoise_get_Ln(FeatureExtractor *st) {
    return &st->Ln[0];
}

float *rnnoise_get_features(FeatureExtractor *st) {
    return &st->features[0];
}

float *rnnoise_get_gain(FeatureExtractor *st) {
    return &st->g[0];
}

float *rnnoise_get_vad(FeatureExtractor *st) {
    return &st->vad;
}

void rnnoise_set_speech_energy(FeatureExtractor *st, float val) {
    st->speech_energy = val * 1.0f;
}

void rnnoise_set_noise_energy(FeatureExtractor *st, float val) {
    st->noise_energy = val * 1.0f;
}

void rnnoise_set_noise_factor(FeatureExtractor *st) {
    if (st->speech_gain != 0) {
        st->noise_factor = pow(10.,
                //(10 * log10(st->speech_energy * st->speech_gain * st->speech_gain / st->noise_energy) -
                               (10 * log10(st->speech_energy / st->noise_energy) - st->expected_snr) / 20);
    }

}

int rnnoise_feature_extractor_init(FeatureExtractor *st) {
    memset(st, 0, sizeof(*st));
    st->a_hp[0] = -1.99599;
    st->a_hp[1] = 0.99600;
    st->b_hp[0] = -2;
    st->b_hp[1] = 1;
    st->speech_gain = 1;
    st->noise_gain = 1;
    st->xst = rnnoise_create(NULL);
    st->nst = rnnoise_create(NULL);
    st->xnst = rnnoise_create(NULL);
    st->expected_snr = 0.0f;
    st->speech_energy = 0.0f;
    st->noise_energy = 0.0f;
    st->noise_factor = 1.0f;
    return 0;
}


FeatureExtractor *rnnoise_feature_extractor_create() {
    FeatureExtractor *st;
    st = (FeatureExtractor *) malloc(sizeof(FeatureExtractor));
    rnnoise_feature_extractor_init(st);
    return st;
}

void rnnoise_feature_extractor_random_change(FeatureExtractor *st) {
    int i;
    st->speech_gain = pow(10., (-20 + (rand() % 30)) / 20.);// speech gain [0.1, 3]
    st->noise_gain = pow(10., (-30 + (rand() % 30)) / 20.);// noise gain [0.1 3]
    if (rand() % 10 == 0) st->noise_gain = 0;
    st->noise_gain *= st->speech_gain;
    if (rand() % 10 == 0) st->speech_gain = 0;

    st->expected_snr = (rand() % 41) - 20;
    //fprintf(stderr, "expected_snr changed to value  =  %f\n", st->expected_snr);
    rand_resp(st->a_noise, st->b_noise);
    rand_resp(st->a_sig, st->b_sig);
#ifdef USE_LOW_PASS
    low_pass = FREQ_SIZE * 3000. / 8000. * pow(50., rand() / (double) RAND_MAX);
    for (i = 0; i < NB_BANDS; i++) {
        if (eband5ms[i] << FRAME_SIZE_SHIFT > low_pass) {
            band_lp = i;
            break;
        }
    }
#endif
}

//#define USE_SNR 1

void rnnoise_feature_extract(FeatureExtractor *st, short *xi, short *ni) {
    int i;
    st->E = 0;
    if (st->speech_gain != 0) {
#if USE_SNR
        for (i = 0; i < FRAME_SIZE; i++) st->x[i] = xi[i];
        for (i = 0; i < FRAME_SIZE; i++) st->E += (xi[i]) * (float) (xi[i]); // no gain involved
#else
        for (i = 0; i < FRAME_SIZE; i++) st->x[i] = st->speech_gain * xi[i];
        for (i = 0; i < FRAME_SIZE; i++) st->E += (xi[i]) * (float) (xi[i]); // no gain involved
#endif
    } else {
        for (i = 0; i < FRAME_SIZE; i++) st->x[i] = 0;
        st->E = 0;
    }

    if (st->noise_gain != 0) {
#if USE_SNR
        //st->noise_factor = pow(10., (10 * log10(st->speech_energy / st->noise_energy) - st->expected_snr) / 20);
        for (i = 0; i < FRAME_SIZE; i++) st->n[i] = st->noise_factor * ni[i];
#else
        for (i = 0; i < FRAME_SIZE; i++) st->n[i] = st->noise_gain * ni[i];
#endif
    } else {
        for (i = 0; i < FRAME_SIZE; i++) st->n[i] = 0;
    }

    // filters
    biquad(st->x, st->mem_hp_x, st->x, st->b_hp, st->a_hp, FRAME_SIZE);
    biquad(st->x, st->mem_resp_x, st->x, st->b_sig, st->a_sig, FRAME_SIZE);
    biquad(st->n, st->mem_hp_n, st->n, st->b_hp, st->a_hp, FRAME_SIZE);
    biquad(st->n, st->mem_resp_n, st->n, st->b_noise, st->a_noise, FRAME_SIZE);

    // x + n -> xn
    for (i = 0; i < FRAME_SIZE; i++) st->xn[i] = st->x[i] + st->n[i];

    // vad_cnt
    //   > e9f   : = 0
    // e8f ~ e9f : -= 5
    // e7f ~ d8f : += 1
    //   < e7f   : += 2
    if (st->E > 1e9f / 3) {
        st->vad_cnt = 0;
    } else if (st->E > 1e8f / 3) {
        st->vad_cnt -= 5;
    } else if (st->E > 1e7f / 3) {
        st->vad_cnt++;
    } else {
        st->vad_cnt += 2;
    }
    // limit to 0~15
    if (st->vad_cnt < 0) st->vad_cnt = 0;
    if (st->vad_cnt > 15) st->vad_cnt = 15;

    // vad_cnt -> vad
    // 0   : 1
    // 1-9 : 0.5
    // 10  : 0
    if (st->vad_cnt >= 10) st->vad = 0;
    else if (st->vad_cnt > 0) st->vad = 0.5f;
    else st->vad = 1.f;

    // x -> Y, Ey(band energy)
    frame_analysis(st->xst, st->Y, st->Ey, st->x);
    // n -> N, En(band energy)
    frame_analysis(st->nst, st->N, st->En, st->n);
    // En -> Ln (log - band Engery)
    for (i = 0; i < NB_BANDS; i++) st->Ln[i] = log10(1e-2 + st->En[i]);

    // as func name, xn -> X(in), P(pitch), E.., vad
    int silence = compute_frame_features(st->xnst, st->X, st->P, st->Ex, st->Ep, st->Exp, st->features, st->xn);
    // pitch filter to X
    pitch_filter(st->X, st->P, st->Ex, st->Ep, st->Exp, st->g);

    // Ex,Ey -> g
    //printf("%f %d\n", noisy->last_gain, noisy->last_period);
    for (i = 0; i < NB_BANDS; i++) {
        st->g[i] = sqrt((st->Ey[i] + 1e-3) / (st->Ex[i] + 1e-3));
        //st->g[i] = sqrt((st->Ey[i] + 1e-3) / (st->En[i] + st->Ey[i] + 1e-3));
        if (st->g[i] > 1) st->g[i] = 1;
        if (silence || i > band_lp) st->g[i] = -1;
        if (st->Ey[i] < 5e-2 && st->Ex[i] < 5e-2) st->g[i] = -1;
        if (st->vad == 0 && st->noise_gain == 0) st->g[i] = -1;
#define VAD_SILENCE 1
#if VAD_SILENCE
        //if (st->vad == 0) st->g[0] = 0.0316;
#endif
    }
}
