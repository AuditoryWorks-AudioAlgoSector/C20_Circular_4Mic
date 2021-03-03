/* Copyright (c) 2017 Mozilla */
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "rnnoise.h"
#include <dirent.h>
#include <time.h>
#include <sys/stat.h>        // 提供属性操作函数
#include <sys/types.h>         // 提供mode_t 类型
#include "dr_wav.h"


#define MAX_FILE 999999
#define MAX_FILE_NAME 1024
#define SAMPLE_RATE 16000
#define DR_WAV_IMPLEMENTATION

static short speech_file_flag[MAX_FILE] = {0};
static short noise_file_flag[MAX_FILE] = {0};

//#if TRAINING
#if 1

float calculate_wav_energy(int16_t *wav_buffer, uint64_t sampleCount) {
    int16_t *input = wav_buffer;
    size_t blocks = sampleCount / BLOCK_SIZE;
    size_t frames = sampleCount / FRAME_SIZE;
    int flag = 0;
    float E_max = 0.00001;
    float E_speech = 0.0;
    if (sampleCount > BLOCK_SIZE) {
        while (flag < blocks) {
            //printf("blocks ==== %d, samplecount = %d, flag == %d\n", blocks, sampleCount, flag);
            for (int j = 0; j < BLOCK_SIZE; j++) {
                E_speech += input[j] * input[j];

            }

            if (E_speech > E_max) E_max = E_speech;
            input += BLOCK_SIZE;
            E_speech = 0.0;
            flag++;
        }
    } else if (sampleCount < BLOCK_SIZE && sampleCount > FRAME_SIZE) {
        while (flag < frames) {
            for (int j = 0; j < FRAME_SIZE; j++) {
                E_speech += input[j] * input[j];
            }

            if (E_speech > E_max) E_max = E_speech;
            input += FRAME_SIZE;
            E_speech = 0.0;
            flag++;
        }
    } else {
        for (int i = 0; i < sampleCount; i++) {
            E_speech = input[i] * input[i];
            if (E_speech > E_max) E_max = E_speech;
        }
    }
    return E_max * 1.0f;
}

void get_file_list(DIR *dir, char *path, char **file_list, int *count) {
    struct dirent *entry;
    char single_file[1024];
    char *pFile;
    char file_path[1024];
    sprintf(file_path, "%s", path);
    //chdir (path);
    while ((entry = readdir(dir)) != NULL) {
        //struct stat stat_buf;
        //lstat(entry->d_name, &stat_buf);
        //if(S_IFDIR & stat_buf.st_mode)    // 判断下一级成员是否是目录
        if (entry->d_type == 4) {
            if (strcmp(".", entry->d_name) == 0 || strcmp("..", entry->d_name) == 0) {
                continue;
            }
            sprintf(file_path, "%s/%s", path, entry->d_name);
            //printf("%s/\n", "", entry->d_name);  // 输出目录名称
            DIR *tmp_dir = opendir(file_path);
            if (tmp_dir != NULL) {
                get_file_list(tmp_dir, file_path, file_list, count);              // 递归调用自身，扫描下一级目录的内容
            }
        } else {
            //printf("%*s%s\n", depth, "", entry->d_name);  // 输出属性不是目录的成
            pFile = strchr(entry->d_name, '.');
            if (pFile != NULL) {
                //printf("file things ==== %s\n", pFile);
                if (strcmp(pFile, ".wav") == 0) {
                    //printf("check %d ===== %s\n", *count, entry->d_name);
                    sprintf(single_file, "%s/%s", path, entry->d_name); //构成文件全路径v
                    //printf("count %d, wav single file full path name === %s \n", *count, single_file);
                    strcpy(file_list[*count], single_file);
                    //printf("file list [%d] = %s\n",*count, file_list[*count]);
                    (*count)++;
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }
    }
    sprintf(file_path, "%s", path);
    //fprintf(stderr, "num of file in directory ====== %d\n", *count);
    fprintf(stderr, "num of file in directory %s ====== %d\n", file_path, *count);
}

drwav_int16 *wavRead_s16(const char *filename, uint32_t *sampleRate, uint64_t *sampleCount, uint32_t *channels) {
    drwav_uint64 totalSampleCount = 0;
    //fprintf(stderr, "read file name ==== %s\n", filename);

    drwav_int16 *input = drwav_open_file_and_read_pcm_frames_s16(filename, channels, sampleRate, &totalSampleCount,
                                                                 NULL);

    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    if (*channels != 1) {
        fprintf(stderr, "read file [%s] error.channels = [%d]\n", filename, *channels);
        exit(1);
    }
    *sampleCount = totalSampleCount * (*channels);
//    for (int32_t i = 0; i < *sampleCount; ++i) {
//        input[i] = input[i] * 32768.0f;
//    }
    return input;
}

float *wavRead_f32(const char *filename, uint32_t *sampleRate, uint64_t *sampleCount, uint32_t *channels) {
    drwav_uint64 totalSampleCount = 0;
    fprintf(stderr, "read file name ==== %s\n", filename);

    float *input = drwav_open_file_and_read_pcm_frames_f32(filename, channels, sampleRate, &totalSampleCount, NULL);
    //fprintf(stderr, "input ==== %d\n", input);
    if (input == NULL) {

    }
    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    *sampleCount = totalSampleCount * (*channels);
    for (int32_t i = 0; i < *sampleCount; ++i) {
        input[i] = input[i] * 32768.0f;
    }
    return input;
}

int speech_rand = 0;
int noise_rand = 0;

int16_t *open_file_s16(char **dir_list, int file_list_len, uint32_t *sampleRate,
                       uint64_t *sampleCount, uint32_t *channels, int flag) {

    //srand((unsigned)time(NULL));
    //int i = rand() % (file_list_len);
    if (flag == 0) {
        speech_rand = rand() % (file_list_len);
    } else {
        noise_rand = rand() % (file_list_len);
    }

//    if (flag ==1) {
//        //printf(" This should be speech file, %d th file === File_name ===  %s \n",i, dir_list[i]);
//    }else{
//        //printf(" This should be noise file , %d th file === File_name ===  %s \n",i, dir_list[i]);
//    }
    //fprintf(stderr, "open file name =================== %s\n", dir_list[i]);
    if (flag == 0) //speech
    {
        speech_file_flag[speech_rand]++;
        //fprintf(stderr, "speech file used: %s\n",dir_list[i]);
    } else {
        noise_file_flag[noise_rand]++;
        //fprintf(stderr, "noise file use: %s\n", dir_list[i]);
    }

    //char *filename = "/media/yongyug/DATA/new_train/noise_dir/noise-free-sound-0711.wav";
    int16_t *buffer = NULL;
    if (flag == 0){
        buffer = wavRead_s16(dir_list[speech_rand], sampleRate, sampleCount, channels);
    } else{
        buffer = wavRead_s16(dir_list[noise_rand], sampleRate, sampleCount, channels);
    }


    //printf("Energy ============= %f\n", *Energy);
    return buffer;
}

float *open_file_f32(char **dir_list, int file_list_len, uint32_t *sampleRate,
                     uint64_t *sampleCount, uint32_t *channels, int flag) {

    //srand((unsigned)time(NULL));
    int i = rand() % (file_list_len);
//    if (flag ==1) {
//        //printf(" This should be speech file, %d th file === File_name ===  %s \n",i, dir_list[i]);
//    }else{
//        //printf(" This should be noise file , %d th file === File_name ===  %s \n",i, dir_list[i]);
//    }
    //fprintf(stderr, "open file name =================== %s\n", dir_list[i]);

    //char *filename = "/media/yongyug/DATA/new_train/noise_dir/noise-free-sound-0711.wav";
    float *buffer = wavRead_f32(dir_list[i], sampleRate, sampleCount, channels);

    //printf("Energy ============= %f\n", *Energy);
    return buffer;
}

#endif

//define TRAINING 1
#if 1

int main(int argc, char **argv) {
    int i;
    int file_rand_i;
    int count = 0;
    int first = 1;
    FILE *f1 = NULL, *f2 = NULL, *fout = NULL;
    int maxCount = 0;
    int gain_change_count = 0;
    short xi[FRAME_SIZE], ni[FRAME_SIZE];
    short tmp[FRAME_SIZE];
    FeatureExtractor *st = rnnoise_feature_extractor_create();

    char **speech_array = (char **) malloc(sizeof(char *) * MAX_FILE);
    char **noise_array = (char **) malloc(sizeof(char *) * MAX_FILE);

    int k;
    for (k = 0; k < MAX_FILE; k++) {
        speech_array[k] = (char *) malloc(sizeof(char) * MAX_FILE_NAME);
        noise_array[k] = (char *) malloc(sizeof(char) * MAX_FILE_NAME);
    }

    int speech_count = 0;
    int noise_count = 0;
    float speech_energy = 0.0f;
    float noise_energy = 0.0f;

    DIR *dirSpeech, *dirNoise = NULL;

    f1 = fopen("clean_speech.raw", "wb");
    f2 = fopen("speech_vad.raw", "wb");

    fprintf(stderr, "This is the random speech version, with static speech gain\n");

    if (argc != 5) {
        fprintf(stderr, "usage: %s <speech_dir> <noise_dir> <count> <output file>\n", argv[0]);
        return 1;
    }
    if ((dirSpeech = opendir(argv[1])) == NULL || (dirNoise = opendir(argv[2])) == NULL) {
        return -1;
    } else {
        get_file_list(dirSpeech, argv[1], speech_array, &speech_count);
        get_file_list(dirNoise, argv[2], noise_array, &noise_count);
    }
    closedir(dirSpeech);
    closedir(dirNoise);
    srand((unsigned) time(NULL));

    fout = fopen(argv[4], "wb");
    maxCount = atoi(argv[3]);

    int16_t *speech_buffer = NULL;
    int16_t *noise_buffer = NULL;
    uint32_t speech_sampleRate = 0;
    uint64_t speech_sampleCount = 0;
    uint32_t speech_channels = 0;
    size_t speech_frames = 0;

    uint32_t noise_sampleRate = 0;
    uint64_t noise_sampleCount = 0;
    uint32_t noise_channels = 0;
    size_t noise_frames = 0;

    short *speech_input = NULL;
    short *noise_input = NULL;

    while (1) {

        if (count == maxCount) break;
        if ((count % 1000) == 0) fprintf(stderr, "%d\r", count);

        // random gain/filter/freq_range for speech and noise
        if (++gain_change_count > 2821) { // 7 x 13 x 31
            rnnoise_feature_extractor_random_change(st);
            gain_change_count = 0;
        }

        if (speech_frames == 0) {
            drwav_free(speech_buffer, NULL);
            speech_buffer = open_file_s16(speech_array, speech_count,
                                          &speech_sampleRate, &speech_sampleCount,
                                          &speech_channels, 0);
            speech_energy = calculate_wav_energy(speech_buffer, speech_sampleCount);

            rnnoise_set_speech_energy(st, speech_energy);
            rnnoise_set_noise_factor(st);
            speech_frames = speech_sampleCount / FRAME_SIZE;
            speech_input = speech_buffer;
        }

        if (noise_frames == 0) {
            drwav_free(noise_buffer, NULL);
            noise_buffer = open_file_s16(noise_array, noise_count,
                                         &noise_sampleRate, &noise_sampleCount,
                                         &noise_channels, 1);
            noise_energy = calculate_wav_energy(noise_buffer, noise_sampleCount);
            rnnoise_set_noise_energy(st, noise_energy);
            rnnoise_set_noise_factor(st);
            noise_frames = noise_sampleCount / FRAME_SIZE;
            noise_input = noise_buffer;
        }

        //memcpy(xi, speech_input, FRAME_SIZE * sizeof(int16_t));
        //memcpy(ni, noise_input, FRAME_SIZE * sizeof(int16_t));

        //rnnoise_feature_extract(st, &xi[0], &ni[0]);
        rnnoise_feature_extract(st, speech_input, noise_input);

        if (rnnoise_get_speech_gain(st) != 0) {
            speech_frames--;
            speech_input += FRAME_SIZE;
        }
        if (rnnoise_get_noise_gain(st) != 0) {
            noise_frames--;
            noise_input += FRAME_SIZE;
        }
        count++;

#define VAD_TEST 0
#if VAD_TEST
        {
            short vad_tmp[FRAME_SIZE];
            for (int j = 0; j < FRAME_SIZE; ++j) {
                float *vad_prob = rnnoise_get_vad(st);
                vad_tmp[j] = 1000 * (*vad_prob);
            }
            fwrite(vad_tmp, sizeof(short), FRAME_SIZE, f2);
            fwrite(speech_input, sizeof(short), FRAME_SIZE, f1);
        };
#endif

#if 1
        fwrite(rnnoise_get_features(st), sizeof(float), NB_FEATURES, fout);//0:73
        fwrite(rnnoise_get_gain(st), sizeof(float), NB_BANDS, fout);//73:108
        //fwrite(rnnoise_get_Ln(st), sizeof(float), NB_BANDS, fout);//108:143
        fwrite(rnnoise_get_vad(st), sizeof(float), 1, fout);//143:144
#endif
    }
    int multi_read = 0;
    int no_read = 0;
    for (i = 0; i < speech_count; ++i) {
        if (speech_file_flag[i] >= 2) {
            multi_read++;
            //fprintf(stderr, "This file was readed %d times, name [%s]\n", speech_file_flag[i], speech_array[i]);
        }
        if (speech_file_flag[i] == 0) {
            no_read++;
        }
    }
    fprintf(stderr, "speech files statics: multi readed files [%d]/[%f], no readed files [%d]/[%f%]\n",
            multi_read, multi_read * 100 / (float) speech_count, no_read, no_read * 100 / (float) speech_count);

    multi_read = 0;
    no_read = 0;
    for (i = 0; i < noise_count; ++i) {
        if (noise_file_flag[i] >= 2) {
            multi_read++;
            //fprintf(stderr, "This file was readed %d times, name [%s]\n", noise_file_flag[i], speech_array[i]);
        }
        if (noise_file_flag[i] == 0) {
            no_read++;
        }
    }
    fprintf(stderr, "noise files statics: multi readed files [%d]/[%f], no readed files [%d]/[%f%]\n",
            multi_read, multi_read * 100 / (float) noise_count, no_read, no_read * 100 / (float) noise_count);

    for (k = 0; k < MAX_FILE; k++) {
        free(speech_array[k]);
        free(noise_array[k]);
    }

    free(speech_array);
    free(noise_array);

    fprintf(stderr, "matrix size: %d x %d\n", count, NB_FEATURES + NB_BANDS + 1);
    drwav_free(noise_buffer, NULL);
    drwav_free(speech_buffer, NULL);
    fclose(f1);
    fclose(f2);
    fclose(fout);
    return 0;
}

#endif
