#include <math.h>
#include "forward.h"
#include "string.h"


void DensewithBN(float *input, int array_len, float *param_w, int num_unit, float *n, float *output) {

    int stride = num_unit;
    for (int i = 0; i < num_unit; i++) {
        float sum = n[i];
        for (int j = 0; j < array_len; j++) {
            sum += input[j] * param_w[j * stride + i];
        }
        output[i] = sum;
    }
}

void Densewithsigmoid(float *input, int array_len,  float *param_w,int num_unit, float *output) {

    int stride = num_unit;
    for (int i = 0; i < num_unit; i++) {
        float sum = 0;
        for (int j = 0; j < array_len; j++) {
            sum += input[j] * param_w[j * stride + i];
        }
        output[i] = 1.0f / (1 + expf(-sum));
    }
}


void prelu(float *input, int num_unit, float *alpha, float *output) {
    for (int i = 0; i < num_unit; i++) {
        if (input[i] > 0)
            output[i] = input[i];
        else
            output[i] = input[i] * alpha[i];
    }
}

void layer_splice(float *inputx, int inputx_dim, float *inputy, int inputy_dim, float *layer_splice_out) {

    memcpy(layer_splice_out, inputx, sizeof(float) * inputx_dim);
    memcpy(layer_splice_out + inputx_dim, inputy, sizeof(float) * inputy_dim);

}

void forward(float *input, float *prediction) {

    float h1_out[N_h1] = {0};
    float h1relu_out[N_h1] = {0};
    float layer_splice_1[NUM_FRAMES * 26 + N_h1] = {0};
    float h2_out[N_h2] = {0};
    float h2relu_out[N_h2] = {0};
    float layer_splice_2[N_h1 + N_h2] = {0};
    float h3_out[N_h3] = {0};
    float h3relu_out[N_h3] = {0};
    float layer_splice_3[N_h2 + N_h3] = {0};


    DensewithBN(input, NUM_FRAMES * 26, h1_combine, N_h1, h1bn_n, h1_out);
    prelu(h1_out, N_h1, h1_alpha, h1relu_out);

    layer_splice(input, NUM_FRAMES * 26, h1relu_out, N_h1, layer_splice_1);

    DensewithBN(layer_splice_1, NUM_FRAMES * 26 + N_h1, h2_combine, N_h2, h2bn_n, h2_out);
    prelu(h2_out, N_h2, h2_alpha, h2relu_out);

    layer_splice(h1relu_out, N_h1, h2relu_out, N_h2, layer_splice_2);

    DensewithBN(layer_splice_2, N_h1 + N_h2, h3_combine, N_h3, h3bn_n, h3_out);
    prelu(h3_out, N_h3, h3_alpha, h3relu_out);

    layer_splice(h2relu_out, N_h2, h3relu_out, N_h3, layer_splice_3);

    Densewithsigmoid(layer_splice_3, N_h2 + N_h3, out_w, NUM_FRAMES, prediction);


}
