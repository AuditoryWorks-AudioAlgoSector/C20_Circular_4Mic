#ifndef PARAMS_H
#define PARAMS_H

#define NUM_FRAMES 5
#define N_h1 64
#define N_h2 64
#define N_h3 64 

extern float h1_alpha[64];
extern float h1_combine[8320];
extern float h1bn_n[64];
extern float h2_alpha[64];
extern float h2_combine[12416];
extern float h2bn_n[64];
extern float h3_alpha[64];
extern float h3_combine[8192];
extern float h3bn_n[64];
extern float out_w[640];

#endif//PARAMS_H