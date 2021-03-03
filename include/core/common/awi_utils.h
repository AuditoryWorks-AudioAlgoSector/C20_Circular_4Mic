#ifndef AUDIO_ALGORITHMS_UTILS_H
#define AUDIO_ALGORITHMS_UTILS_H

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef PI
#define PI 3.1415926535897932f
#endif

#ifndef AWI_MAX
#define AWI_MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif // AWI_MAX

#ifndef AWI_MIN
#define AWI_MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif // AWI_MIN
#ifdef __cplusplus
extern "C" {
#endif

float awi_calculate_vss(float corr);

float awi_sigmoid(float corr, int midPoint, int slope);

float awi_db2pow(float ydB);

float awi_get_radius(int sampling_rate);

void awi_filter_dc_notch16(const float *in, float radius, float *out, int len, float *mem);

void awi_cross_fadef(float *output, float *newChannel, float *oldChannel, int n, int id, int totalFrame);

int awi_is_active(float *arr, int length, float threshold);

void awi_normalize(float *dst, short *src, int length);
void awi_unnormalize(short *dst, float *src, int length);

int awi_get_max_index_float(float *arr, int length);

int awi_get_min_index_float(float *arr, int length);

float awi_get_max_float(float *arr, int length);

float awi_get_min_float(float *arr, int length);

float awi_get_mean_float(float *arr, int length);

int awi_get_max_index_int(int *arr, int length);

float awi_get_median_float(float *a, int n, int i);

void awi_insertion_sort_float(float *a, int n);

int awi_median_partition_float(float *a, int n, float x);

void awi_correct_delay(float *buffer, float *data, int dataSize, int channel, int delaySize);


float NormC(float *z);

void CXC(float *z, float *x, float *y);

void CPlusC(float *z, float *x, float *y);

void CSubtractC(float *z, float *x, float *y);

void awi_init_step_size(float *arr, int bins);

float fminf_local(float a, float b);

float fmaxf_local(float a, float b);

void awi_smoothing_1D_vector(float *dst, const float *src, const float *window, int data_length, int window_length);

#ifdef __cplusplus
};
#endif

#endif //AUDIO_ALGORITHMS_UTILS_H
