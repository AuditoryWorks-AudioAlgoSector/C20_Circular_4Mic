#ifndef AUDIO_ALGORITHMS_FILTERBANK_H
#define AUDIO_ALGORITHMS_FILTERBANK_H
#include "awi_constant.h"

#ifdef __cplusplus
extern "C" {
#endif

extern float FILTER_256[ANA_FILTER_LENGTH_256];

extern float AWI_F[FRAME_BAND_COUNTS_256];

extern float AWI_CRITICAL_BAND[AWI_FREQ_BANDS * 3];

extern int AWI_FREQ_LOW[AWI_FREQ_BANDS];

extern int AWI_FREQ_HIGH[AWI_FREQ_BANDS];

extern int AWI_AES_TABLE[FRAME_BAND_COUNTS_256];

extern float AWI_NSE_ALPHA_UPDATE[FRAME_BAND_COUNTS_256];

extern float DOA_BINS[AWI_DOA_BIN_NUM * 3];

extern int WPE_FREQ_TAPS[FRAME_BAND_COUNTS_256];

extern int FREQ_TAPS[FRAME_BAND_COUNTS_256];

extern int TAP_FREQS[AWI_AEC_TAP];

extern float EQ_COMPEN[FRAME_BAND_COUNTS_256];

extern float LOCAL_SMOOTHING_WINDOW[NS_LOCAL_WIN_LEN];

extern float GLOBAL_SMOOTHING_WINDOW[NS_GLOBAL_WIN_LEN];

#ifdef __cplusplus
};
#endif

#endif
