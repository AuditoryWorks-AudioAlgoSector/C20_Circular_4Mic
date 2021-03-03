#ifndef _awi_abf_H_
#define _awi_abf_H_

#include "awi_abf_cfg.h"
#include "awi_filterbankframe.h"
#include "awi_constant.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct awi_abf_part1
{
    awi_abf_cfg_t cfg;
    float wfb[AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * AWI_DOA_BIN_NUM]; // 18.14kb
    float sd[AWI_FRAME_BAND_COUNT * AWI_DOA_BIN_NUM]; // 3.02kb
    float se[AWI_FRAME_BAND_COUNT * AWI_DOA_BIN_NUM]; // 3.02kb
    float sx[AWI_FRAME_BAND_COUNT * AWI_ABF_REFERENCE * AWI_DOA_BIN_NUM];         // 9.07kb
    float recur_snr[AWI_FRAME_BAND_COUNT * AWI_ABF_REFERENCE * AWI_DOA_BIN_NUM];  // 9.07kb
    float recur_global_snr[AWI_ABF_REFERENCE * AWI_DOA_BIN_NUM];                  // 0.07kb

} awi_abf_part1_t;


typedef struct awi_abf_part2
{
    float xfm[AWI_FRAME_BAND_COUNT * AWI_ABF_TAP * AWI_ABF_REFERENCE * 2 * AWI_DOA_BIN_NUM]; // 18.14kb
    float mEk[AWI_FRAME_BAND_COUNT * 2]; // 1.01kb
    float yfm[AWI_FRAME_BAND_COUNT * 2]; // 1.01kb

} awi_abf_part2_t;



void awi_abf_init(awi_abf_part1_t *p1, awi_abf_part2_t *p2);

void awi_abf_process(awi_abf_part1_t *p1, awi_abf_part2_t *p2, float *frameIn, float *ref1_sp, float *ref2_sp,
        float *ref3_sp, float *frameOut, int beamId);

void awi_abf_reset(awi_abf_part1_t *p1, awi_abf_part2_t *p2);


#ifdef __cplusplus
};
#endif

#endif
