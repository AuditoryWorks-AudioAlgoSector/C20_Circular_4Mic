//
// Created by chado on 2020/8/17.
//
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include "rnnoise.h"
#include "dr_wav.h"
#include <sys/time.h>
#include <iostream>
#include "modules/audio_processing/agc/legacy/gain_control.h"
#include "../utility/media_buffer.h"
#include <unistd.h>

#define FRAME_SIZE 160
//#define FRAME_SIZE 480
#define SAMPLE_RATE 16000
#define SECONDS_PER_HOUR (3600)
#define MINUTES_PER_HOUR 60
#define SECONDS_PER_MINUTE 60
int16_t rnnoise_samples[SAMPLE_RATE * SECONDS_PER_HOUR / 60] = {0};
int16_t agc_samples[SAMPLE_RATE * SECONDS_PER_HOUR / 60] = {0};
using namespace webrtc;

size_t on_file_seek_and_write(void *pUserData, const void *pData, size_t bytesToWrite) {
    return bytesToWrite;
}

int wavWrite_int16(char *filename, int16_t *buffer, int sampleRate, uint32_t totalSampleCount) {
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_PCM;
    format.channels = 1;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 16;
    drwav pWav;
    drwav_init_file_write(&pWav, filename, &format, NULL);

    drwav_uint64 samplesWritten = drwav_write_raw(&pWav, totalSampleCount * sizeof(int16_t), buffer);
    drwav_uninit(&pWav);
    if (samplesWritten != totalSampleCount * 2) {
        fprintf(stderr, "ERROR\n");
        return -1;
    }
    return 0;
}

awi::MediaBuffer buffer(6400);
bool audio_process_runnable = true;

//typedef void(int i) call_list;
int audio_emit(int n) {
    printf("%d samples audio buffer emitted\n");
    return 0;
}

int audio_process(void *callback) {
    int i;
    int first = 1;
    float x[FRAME_SIZE]; //input
    FILE *f1, *fout;
    FILE *fvad = NULL;
    DenoiseState *st;
    char wav_header[44];
    short vad_tmp[FRAME_SIZE];

    struct timeval t1, t2, t3, t4;
    uint32_t totalSampleCount = 0;

    printf("exe start11111111111111111111111111111");
    //那么函数f运行所花的时间为
    //deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec 微秒
    st = rnnoise_create(NULL);

    int16_t split_band_data[3][FRAME_SIZE];
    int16_t *split_bands[3] = {split_band_data[0], split_band_data[1], split_band_data[2]};
    int32_t gains[11] = {0};
    std::cout << "Hello World!\n";
    void *agcInst = WebRtcAgc_Create();
    WebRtcAgc_Init(agcInst, 0, 255, 2, 16000);

    WebRtcAgcConfig agcConfig;
    agcConfig.compressionGaindB = 20;//在Fixed模式下，越大声音越大
    agcConfig.limiterEnable = 1;
    agcConfig.targetLevelDbfs = 3;   //dbfs表示相对于full scale的下降值，0表示full scale，越小声音越大
    WebRtcAgc_set_config(agcInst, agcConfig);

    int32_t capture_levels = 0;
    int32_t new_capture_level = 0;
    uint8_t saturation_warning = 0;
    bool stream_has_echo = false;
    int32_t  count = 1;
    while (audio_process_runnable) {

        //int len = fread(split_bands[0], sizeof(short), FRAME_SIZE, input);
        if (buffer.readableSize() >= 320){
            buffer.read((int8_t *) split_bands[0], 320);
        } else{
            usleep(5000);
            continue;
        }


        for (i = 0; i < FRAME_SIZE; i++) x[i] = split_bands[0][i];
        gettimeofday(&t1, NULL);
        float vad_prob = rnnoise_process_frame(st, x, x);
        gettimeofday(&t2, NULL);

        for (int j = 0; j < FRAME_SIZE; ++j) {
            vad_tmp[j] = vad_prob * 1000;
        }
        if (!first && fvad != NULL) fwrite(vad_tmp, sizeof(short), FRAME_SIZE, fvad);

        for (i = 0; i < FRAME_SIZE; i++) split_bands[0][i] = x[i];
        //if (!first) fwrite(tmp, sizeof(short), FRAME_SIZE, fout);
        memcpy(rnnoise_samples + totalSampleCount,
               (char *) split_bands[0],
               sizeof(int16_t) * FRAME_SIZE);

        gettimeofday(&t3, NULL);
        WebRtcAgc_Analyze(
                agcInst, split_bands, 1,
                FRAME_SIZE, capture_levels, &new_capture_level,
                stream_has_echo, &saturation_warning, gains);

        capture_levels = new_capture_level;

        WebRtcAgc_Process(agcInst, gains,
                          (const int16_t *const *) split_bands,
                          1,
                          (int16_t *const *) split_bands);
        gettimeofday(&t4, NULL);
        long deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        long deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
        printf("\rframe num %d，RNNoise process one frame(10ms) cost %ld us , Agc cost %ld us\r",
                count, deltaT1, deltaT2);
        //fwrite(split_bands[0], sizeof(short), FRAME_SIZE, output);
        memcpy(agc_samples + totalSampleCount,
               (char *) split_bands[0],
               sizeof(int16_t) * FRAME_SIZE);
        first = 0;
        count++;
        totalSampleCount += 160;
        usleep(5000 - (deltaT1 + deltaT2));
    }
    wavWrite_int16("rnnoise_out.wav", rnnoise_samples, SAMPLE_RATE, totalSampleCount);
    wavWrite_int16("agc_out.wav", agc_samples, SAMPLE_RATE, totalSampleCount);

    rnnoise_destroy(st);
    if (f1 != NULL) fclose(f1);
    if (fvad != NULL)fclose(fvad);
    std::cout << "Hello World!\n";
    return 1;

}

int main(int argc, char **argv) {
#define FRAME_MS 16
#define FRAME_16_SIZE 16*FRAME_MS
    short tmp[FRAME_16_SIZE];

    std::thread audio_process_thread(audio_process, (void *) &audio_emit);
    //audio_process_thread.detach();
    FILE *input = fopen(argv[1], "rb");
    while (1) {
        int len = fread(tmp, sizeof(short), FRAME_16_SIZE, input);
        if (len < FRAME_16_SIZE) {
            printf("frame stop");
            audio_process_runnable = false;
            break;
        }
        buffer.write((const int8_t *) tmp, len * 2);
        usleep(16000);
    }
    audio_process_runnable = false;
    if (audio_process_thread.joinable()) {
        audio_process_thread.join();
    }
    return 0;
}


int process_frame(int argc, char **argv) {
    int i;
    int first = 1;
    float x[FRAME_SIZE]; //input
    FILE *f1, *fout;
    FILE *fvad = NULL;
    DenoiseState *st;
    char wav_header[44];
    short vad_tmp[FRAME_SIZE];

    struct timeval t1, t2, t3, t4;
    uint32_t totalSampleCount = 0;

    printf("exe start11111111111111111111111111111");
    //那么函数f运行所花的时间为
    //deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec 微秒
    st = rnnoise_create(NULL);
    if (argc <= 1) {
        fprintf(stderr, "usage: %s <noisy speech> <output denoised>\n", argv[0]);
        return 1;
    }

    if (strstr(argv[1], ".pk")) {
        fprintf(stderr, "\rthis is a pk file, %s\n", argv[1]);
        return -1;
    }

    if (argc == 4) {

        if (strcmp(argv[3], "vad") == 0) {
            char filename[256];
            sprintf(filename, "%s_vad.raw", argv[2]);
            fvad = fopen(filename, "wb");
        }
    }

    FILE *input = fopen(argv[1], "rb");
    FILE *output = fopen(argv[2], "wb");
    int16_t split_band_data[3][FRAME_SIZE];
    int16_t *split_bands[3] = {split_band_data[0], split_band_data[1], split_band_data[2]};
    int32_t gains[11] = {0};
    std::cout << "Hello World!\n";
    void *agcInst = WebRtcAgc_Create();
    WebRtcAgc_Init(agcInst, 0, 255, 2, 16000);

    WebRtcAgcConfig agcConfig;
    agcConfig.compressionGaindB = 20;//在Fixed模式下，越大声音越大
    agcConfig.limiterEnable = 1;
    agcConfig.targetLevelDbfs = 3;   //dbfs表示相对于full scale的下降值，0表示full scale，越小声音越大
    WebRtcAgc_set_config(agcInst, agcConfig);

    int32_t capture_levels = 0;
    int32_t new_capture_level = 0;
    uint8_t saturation_warning = 0;
    bool stream_has_echo = false;

    fread(wav_header, sizeof(char), 44, input);
    while (!feof(input)) {

        int len = fread(split_bands[0], sizeof(short), FRAME_SIZE, input);
        if (len < FRAME_SIZE) {
            break;
        }
        for (i = 0; i < FRAME_SIZE; i++) x[i] = split_bands[0][i];
        gettimeofday(&t1, NULL);
        float vad_prob = rnnoise_process_frame(st, x, x);
        gettimeofday(&t2, NULL);


        for (int j = 0; j < FRAME_SIZE; ++j) {
            vad_tmp[j] = vad_prob * 1000;
        }
        if (!first && fvad != NULL) fwrite(vad_tmp, sizeof(short), FRAME_SIZE, fvad);

        for (i = 0; i < FRAME_SIZE; i++) split_bands[0][i] = x[i];
        //if (!first) fwrite(tmp, sizeof(short), FRAME_SIZE, fout);
        memcpy(rnnoise_samples + totalSampleCount,
               (char *) split_bands[0],
               sizeof(int16_t) * FRAME_SIZE);

        gettimeofday(&t3, NULL);
        WebRtcAgc_Analyze(
                agcInst, split_bands, 1,
                FRAME_SIZE, capture_levels, &new_capture_level,
                stream_has_echo, &saturation_warning, gains);

        capture_levels = new_capture_level;

        WebRtcAgc_Process(agcInst, gains,
                          (const int16_t *const *) split_bands,
                          1,
                          (int16_t *const *) split_bands);
        gettimeofday(&t4, NULL);
        long deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        long deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
        fprintf(stderr, "\rRNNoise process one frame(10ms) cost %ld us , Agc cost %ld us\r", deltaT1, deltaT2);
        //fwrite(split_bands[0], sizeof(short), FRAME_SIZE, output);
        memcpy(agc_samples + totalSampleCount,
               (char *) split_bands[0],
               sizeof(int16_t) * FRAME_SIZE);
        first = 0;
        totalSampleCount += len;
    }
    char rnnoise_out[256] = {0};
    char agc_out[256] = {0};
    sprintf(rnnoise_out, "rnnoise_out_%s", argv[1]);
    wavWrite_int16(rnnoise_out, rnnoise_samples, SAMPLE_RATE, totalSampleCount);
    sprintf(agc_out, "agc_out_%s", argv[1]);
    wavWrite_int16(agc_out, agc_samples, SAMPLE_RATE, totalSampleCount);

    rnnoise_destroy(st);
    if (f1 != NULL) fclose(f1);
    if (fvad != NULL)fclose(fvad);
    std::cout << "Hello World!\n";
    return 1;

}
