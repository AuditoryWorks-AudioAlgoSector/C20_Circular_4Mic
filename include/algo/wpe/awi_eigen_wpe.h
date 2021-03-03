#ifndef _AWI_EIGEN_WPE_H_
#define _AWI_EIGEN_WPE_H_

// #include <Eigen/Dense>
// #include <Eigen/Eigen>
#include "awi_wpe_cfg.h"
#include <stdbool.h>

#ifdef __cplusplus
 extern "C"{
 #endif

typedef struct eigen_wpe awi_eigen_wpe_t;

awi_eigen_wpe_t* awi_eigen_wpe_create();

void awi_eigen_wpe_init(awi_eigen_wpe_t *wpe, awi_wpe_cfg_t *wpe_cfg);

void awi_eigen_wpe_process(awi_eigen_wpe_t *wpe, float *input_sp, float *out_sp);

void awi_eigen_wpe_process_freq_dep(awi_eigen_wpe_t *wpe, float *input_sp, float *out_sp);


#ifdef __cplusplus
 };
 #endif

#endif

