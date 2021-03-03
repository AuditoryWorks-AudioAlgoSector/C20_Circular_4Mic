#ifndef AUDIO_ALGORITHMS_AUDIOFRAME_H
#define AUDIO_ALGORITHMS_AUDIOFRAME_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "awi_constant.h"


typedef float audio_mic_frame[AWI_FRAME_LENGTH * AWI_MIC_CHANNEL];

typedef float audio_spk_frame[AWI_FRAME_LENGTH * AWI_SPK_CHANNEL];

#endif /* AUDIO_ALGORITHMS_AUDIOFRAME_H */
