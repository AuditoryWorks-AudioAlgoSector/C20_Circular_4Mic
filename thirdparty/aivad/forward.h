#ifndef FORWARD_H
#define FORWARD_H

#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

void DensewithBN(float *input, int array_len, float *param_w, int num_unit, float *n, float *output);

void layer_splice(float *inputx, int inputx_dim, float *inputy, int inputy_dim, float *layer_splice_out);

void prelu(float *input, int num_unit, float *alpha, float *output);

void Densewithsigmoid(float *input, int array_len, float *param_w, int num_unit, float *output);

void forward(float *input, float *prediction);

#ifdef __cplusplus
};
#endif

#endif//FORWARD_H
