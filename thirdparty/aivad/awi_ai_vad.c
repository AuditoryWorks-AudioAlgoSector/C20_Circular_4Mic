#include <stdlib.h>
#include "awi_ai_vad.h"


void awi_ai_vad_init(awi_ai_vad_t *ai_vad, awi_ai_vad_cfg_t *cfg) {

    ai_vad->cfg = cfg;
    rb_s16_init(&ai_vad->temp_inputdata_frames, cfg->tempFrames);
    // init mfcc features
    rb_init(&ai_vad->num_frames_mfcc, cfg->num_mfcc_frames * cfg->num_mfcc_features);
    float add_mfcc_featres[(NUM_FRAMES - 1) / 2 * 13] = {-36.043653389f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                         -36.043653389f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    rb_write(&ai_vad->num_frames_mfcc, add_mfcc_featres, (NUM_FRAMES - 1) / 2 * 13);
    //
    rb_init(&ai_vad->pred_array, cfg->num_total_pred);
    float add_pred[NUM_FRAMES * (NUM_FRAMES - 1) / 2] = {0};
    rb_write(&ai_vad->pred_array, add_pred, NUM_FRAMES * (NUM_FRAMES - 1) / 2);

    rb_init(&ai_vad->one_frames_input, cfg->frame_size);
    float add_zeros[256] = {0};
    rb_write(&ai_vad->one_frames_input, add_zeros, cfg->frame_step);

    ai_vad->temp_feature_input = (float *) malloc(sizeof(float) * cfg->num_mfcc_features);
    ai_vad->delta = (float *) malloc(sizeof(float) * cfg->num_mfcc_frames * cfg->num_mfcc_features);
    ai_vad->features = (float *) malloc(sizeof(float) * cfg->num_mfcc_frames * cfg->features_size);
    ai_vad->input_frame = (float *) malloc(sizeof(float) * cfg->frame_size);
    ai_vad->total_mfcc_out = (float *) malloc(sizeof(float) * cfg->num_mfcc_frames * cfg->num_mfcc_features);
    ai_vad->mfcc_features = (float *) malloc(sizeof(float) * cfg->num_mfcc_features);
    ai_vad->powspectrum = (float *) malloc((cfg->nfft / 2 + 1) * sizeof(float));
    ai_vad->local_data = (float *) malloc(sizeof(float) * cfg->frame_size);

    memset(ai_vad->temp_feature_input, 0, sizeof(float) * cfg->num_mfcc_features);
}

short awi_ai_vad_process(awi_ai_vad_t *ai_vad, short *input, short *output) {

    float step_pred[NUM_FRAMES] = {0};
    float total_pred[NUM_FRAMES * NUM_FRAMES] = {0};
    float prediction = 0;
    float input_256[256] = {0};
    for (int i = 0; i < ai_vad->cfg->frame_step; i++) {
        input_256[i] = input[i] * 3.0517578125e-5f;
    }

    if (ai_vad->pred_array.valid_size < (NUM_FRAMES - 1) * NUM_FRAMES) {
        if (ai_vad->num_frames_mfcc.valid_size < 13 * (NUM_FRAMES - 1)) {
            // save the input frame data
            rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);
            // calculate and save mfcc features
            input_to_mfcc(ai_vad, input_256);
            return -1;
        } else {
            // save the input frame data
            rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);
            input_to_mfcc(ai_vad, input_256);
            // read five frames mfcc features and calculate prediction
            mfcc_to_pred(ai_vad, step_pred);
            return -1;
        }
    } else {
        // save the input frame data
        rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);
        input_to_mfcc(ai_vad, input_256);
        mfcc_to_pred(ai_vad, step_pred);
        // read all five framse pred value for calculate the result
        rb_read(&ai_vad->pred_array, total_pred, NUM_FRAMES * NUM_FRAMES);
        // calculate prediction
        for (int i = 1; i <= NUM_FRAMES; i++) {
            prediction = prediction + total_pred[(NUM_FRAMES - 1) * i];
        }
        prediction = prediction / NUM_FRAMES;
        prediction = (prediction > ai_vad->cfg->th) ? 1 : 0;
        ai_vad->cfg->vad_flag = (short) prediction;
        // output the predicted frame data after completing the first complete nn calculation
        rb_s16_read(&ai_vad->temp_inputdata_frames, output, ai_vad->cfg->frame_step);
    }
    return ai_vad->cfg->vad_flag;
}

void input_to_mfcc(awi_ai_vad_t *ptr, float *input) {

    rb_write(&ptr->one_frames_input, input, ptr->cfg->frame_step);
    // read the input data in rb and calculate mfcc
    rb_read(&ptr->one_frames_input, ptr->input_frame, ptr->cfg->frame_size);
    // calculate mfcc features
    get_mfcc(ptr, ptr->input_frame, ptr->mfcc_features);
    // save the obtained mfcc features into rb
    rb_write(&ptr->num_frames_mfcc, ptr->mfcc_features, ptr->cfg->num_mfcc_features);
}

//
void mfcc_to_pred(awi_ai_vad_t *ptr, float *step_pred) {
    // load five frames mfcc feature
    rb_read(&ptr->num_frames_mfcc, ptr->total_mfcc_out, NUM_FRAMES * ptr->cfg->num_mfcc_features);
    // calculate delta features
    get_delta(ptr,ptr->total_mfcc_out, ptr->delta, ptr->temp_feature_input);
    // concatenate mfcc and delta features
    concatenate(ptr, ptr->total_mfcc_out, ptr->delta, ptr->features);
    // nn inference
    forward(ptr->features, step_pred);
    // write pred result to pred_array
    rb_write(&ptr->pred_array, step_pred, NUM_FRAMES);
}


void get_mfcc(awi_ai_vad_t *ptr, float *input, float *mfcc_out) {
    // pre-emphasis
    for (int i = ptr->cfg->frame_size - 1; i > 0; i--) {
        ptr->local_data[i] = input[i] * 1.0f - input[i - 1] * ptr->cfg->Preemphasis;
    }
    // save temp value
    ptr->local_data[0] = 1.0f * input[0] - ptr->cfg->preem_tempdata * ptr->cfg->Preemphasis;
    ptr->cfg->preem_tempdata = input[ptr->cfg->frame_step - 1];
    // cal mfcc features
    mfcc(mfcc_out, ptr->powspectrum, ptr->local_data, ptr->cfg->sr, ptr->cfg->nfilt, ptr->cfg->numcep, ptr->cfg->nfft,
         ptr->cfg->ceplifter, ptr->cfg->appendEnergy);
}

void get_delta(awi_ai_vad_t *ptr, float *input, float *delta, float *temp_feature_input) {
    // calculate delta features
    // delta = f(x)-f(x-1)
    // i=0
    for (int j = 0; j < ptr->cfg->num_mfcc_features; j++) {
        delta[j] = input[j] - temp_feature_input[j];
        // save temp value for the next calculation
        temp_feature_input[j] = input[j];
    }
    for (int i = 1; i < NUM_FRAMES; i++) {
        for (int j = 0; j < ptr->cfg->num_mfcc_features; j++) {
            delta[i * ptr->cfg->num_mfcc_features + j] =
                    input[i * ptr->cfg->num_mfcc_features + j] - input[(i - 1) * ptr->cfg->num_mfcc_features + j];
        }
    }
}

void concatenate(awi_ai_vad_t *ptr, float *mfcc, float *delta, float *features) {
    // concatenate mfcc and delta
    for (int i = 0; i < ptr->cfg->num_mfcc_frames; i++) {
        for (int j = 0; j < ptr->cfg->num_mfcc_features; j++) {
            features[i * ptr->cfg->features_size + j] = mfcc[i * ptr->cfg->num_mfcc_features + j];
            features[i * ptr->cfg->features_size + j + 13] = delta[i * ptr->cfg->num_mfcc_features + j];
        }
    }
}

short awi_ai_vad_access(awi_ai_vad_t *ai_vad, short *input, short *output){
    rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);
    if (ai_vad->pred_array.valid_size < (NUM_FRAMES - 1)*NUM_FRAMES) {
        return -1;
    }else{
        rb_s16_read(&ai_vad->temp_inputdata_frames, output, ai_vad->cfg->frame_step);
    }
    return ai_vad->cfg->vad_flag;
}


short awi_ai_vad_determine(awi_ai_vad_t *ai_vad, float *input){
    float step_pred[NUM_FRAMES] = {0};
    float total_pred[NUM_FRAMES*NUM_FRAMES] = {0};
    float prediction = 0;
    //float input_256[256] = {0};
    // for (int i = 0; i < ai_vad->cfg->frame_step; i++) {
    //     input_256[i] = input[i];
    // }

    if (ai_vad->pred_array.valid_size < (NUM_FRAMES - 1)*NUM_FRAMES) {
        if (ai_vad->num_frames_mfcc.valid_size < 13 *(NUM_FRAMES - 1) ) {
            //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
            input_to_mfcc(ai_vad, input);
        } else {
            //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
            input_to_mfcc(ai_vad, input);
            mfcc_to_pred(ai_vad, step_pred);
        }
        return -1;
    } else {
        //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
        input_to_mfcc(ai_vad, input);
        mfcc_to_pred(ai_vad, step_pred);
        rb_read(&ai_vad->pred_array, total_pred,NUM_FRAMES*NUM_FRAMES);
        for (int i = 1; i <= NUM_FRAMES ; i++) {
            prediction = prediction + total_pred[(NUM_FRAMES - 1)*i];
        }
        prediction = prediction / NUM_FRAMES;
        prediction = (prediction > ai_vad->cfg->th) ? 1 : 0;
        ai_vad->cfg->vad_flag = (short) prediction;
        //rb_s16_read(&ai_vad->temp_inputdata_frames, output, ai_vad->cfg->frame_step);//完成第一次完整NN计算后将预测的帧输出
    }
    return ai_vad->cfg->vad_flag;
}


float awi_ai_vad_determine_f(awi_ai_vad_t *ai_vad, float *input){
    float step_pred[NUM_FRAMES] = {0};
    float total_pred[NUM_FRAMES*NUM_FRAMES] = {0};
    float prediction = 0;
    //float input_256[256] = {0};
    // for (int i = 0; i < ai_vad->cfg->frame_step; i++) {
    //     input_256[i] = input[i];
    // }

    if (ai_vad->pred_array.valid_size <(NUM_FRAMES - 1)*NUM_FRAMES) {
        if (ai_vad->num_frames_mfcc.valid_size < 13 * (NUM_FRAMES - 1)) {
            //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
            input_to_mfcc(ai_vad, input);
        } else {
            //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
            input_to_mfcc(ai_vad, input);
            mfcc_to_pred(ai_vad, step_pred);
        }
        return -1;
    } else {
        //rb_s16_write(&ai_vad->temp_inputdata_frames, input, ai_vad->cfg->frame_step);//将输入的帧保存起来
        input_to_mfcc(ai_vad, input);
        mfcc_to_pred(ai_vad, step_pred);
        rb_read(&ai_vad->pred_array, total_pred,NUM_FRAMES*NUM_FRAMES);
        for (int i = 1; i <= NUM_FRAMES; i++) {
            prediction = prediction + total_pred[(NUM_FRAMES - 1)*i];
        }
        prediction = prediction / NUM_FRAMES;
        ai_vad->cfg->vad_flag = (short)(prediction > ai_vad->cfg->th) ? 1 : 0;
        //rb_s16_read(&ai_vad->temp_inputdata_frames, output, ai_vad->cfg->frame_step);//完成第一次完整NN计算后将预测的帧输出
    }
    return prediction;
}

void awi_ai_vad_deinit(awi_ai_vad_t *ai_vad) {

    rb_deinit(&ai_vad->num_frames_mfcc);
    rb_deinit(&ai_vad->pred_array);
    rb_deinit(&ai_vad->one_frames_input);
    rb_s16_deinit(&ai_vad->temp_inputdata_frames);
    free(ai_vad->total_mfcc_out);
    free(ai_vad->input_frame);
    free(ai_vad->powspectrum);
    free(ai_vad->local_data);
    free(ai_vad->features);
    free(ai_vad->delta);
    free(ai_vad->temp_feature_input);
}

