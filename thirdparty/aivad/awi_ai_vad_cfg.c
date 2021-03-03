#include "awi_ai_vad_cfg.h"


void awi_ai_vad_cfg_init(awi_ai_vad_cfg_t *cfg) {

    cfg->th = 0.5f;
    cfg->num_mfcc_frames = NUM_FRAMES;
    cfg->tempFrames = NUM_FRAMES *256;
    cfg->vad_flag = 0;
    cfg->num_total_pred = NUM_FRAMES*NUM_FRAMES;
    cfg->num_mfcc_features = 13;
    cfg->features_size = 26;
    cfg->nfilt = 26;
    cfg->numcep = 13;
    cfg->ceplifter = 22;
    cfg->nfft = 512;
    cfg->sr = 16000;
    cfg->frame_size = 512;
    cfg->frame_step = 256;
    cfg->appendEnergy = 0;
    cfg->preem_tempdata = 0.0f;
    cfg->Preemphasis = 0.97f;
}
