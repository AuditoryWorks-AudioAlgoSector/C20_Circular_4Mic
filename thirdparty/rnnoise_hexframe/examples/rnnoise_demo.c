/* Copyright (c) 2018 Gregor Richards
 * Copyright (c) 2017 Mozilla */
/*
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include "rnnoise.h"
#include "dr_wav.h"
#ifdef linux
#include <sys/time.h>
#elif defined(WIN32)
#include <stdint.h>
#include <time.h>
#else
#include <stdint.h>
#include <time.h>
#endif

//#define FRAME_SIZE 480
#define SAMPLE_RATE 16000
#define SECONDS_PER_HOUR (3600)
int16_t samples[SAMPLE_RATE * SECONDS_PER_HOUR] = {0};

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

#if 1//!TRAINING
int main(int argc, char **argv) {
    int i;
    int first = 1;
    float x[FRAME_SIZE];
    FILE *f1, *fout;
    FILE *fvad = NULL;
    DenoiseState *st;
    char wav_header[44];
    short vad_tmp[FRAME_SIZE];

#ifdef linux
    struct timeval t1, t2;
#endif
    uint32_t totalSampleCount = 0;
    //那么函数f运行所花的时间为
    //deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec 微秒
    st = rnnoise_create(NULL);
    if (argc <= 2) {
        fprintf(stderr, "usage: %s <noisy speech> <output denoised>\n", argv[0]);
        return 1;
    }

    if (strstr(argv[1],".pk")){
        fprintf(stderr, "\rthis is a pk file, %s\n", argv[1]);
        return -1;
    }

    f1 = fopen(argv[1], "rb");
    //fout = fopen(argv[2], "wb");

    if (argc == 4) {

        if (strcmp(argv[3], "vad") == 0) {
            char filename[256];
            sprintf(filename, "%s_vad.raw", argv[2]);
            fvad = fopen(filename, "wb");
        }
    }

    fread(wav_header, sizeof(char), 44, f1);
    while (1) {
        short tmp[FRAME_SIZE];
        int len = fread(tmp, sizeof(short), FRAME_SIZE, f1);
        if (feof(f1)) break;
        for (i = 0; i < FRAME_SIZE; i++) x[i] = tmp[i];
#ifdef linux
        gettimeofday(&t1, NULL);
#endif
        float vad_prob = rnnoise_process_frame(st, x, x);
#ifdef linux
        gettimeofday(&t2, NULL);
        long deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        fprintf(stderr, "\rRNNoise process one frame(10ms) cost %ld u\r", deltaT);
#endif
        for (int j = 0; j < FRAME_SIZE; ++j) {
            vad_tmp[j] = vad_prob * 1000;
        }

        if (!first && fvad != NULL) fwrite(vad_tmp, sizeof(short), FRAME_SIZE, fvad);

        for (i = 0; i < FRAME_SIZE; i++) tmp[i] = x[i];
        //if (!first) fwrite(tmp, sizeof(short), FRAME_SIZE, fout);
        memcpy(samples + totalSampleCount, (char *) tmp, sizeof(int16_t) * FRAME_SIZE);

        first = 0;
        totalSampleCount += len;
    }

    wavWrite_int16(argv[2], samples, SAMPLE_RATE, totalSampleCount);

    rnnoise_destroy(st);
    fclose(f1);
    //fclose(fout);
    if (fvad != NULL)fclose(fvad);


    return 0;
}

#endif
