// rnnoise_train.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <stdio.h>
#include "rnnoise.h"

#if TRAINING

int main(int argc, char **argv) {
    int i;
    int count = 0;
    FILE *f1, *f2, *fout;
    int maxCount = 0;
    int gain_change_count = 0;
    short xi[FRAME_SIZE], ni[FRAME_SIZE];
    FeatureExtractor *st = rnnoise_feature_extractor_create();
    // command line para
    if (argc != 5) {
        fprintf(stderr, "usage: %s <speech> <noise> <count> <output file>\n", argv[0]);
        return 1;
    }
    fprintf(stderr, "%s,%s\n",argv[1],argv[2]);
    f1 = fopen(argv[1], "rb");
    f2 = fopen(argv[2], "rb");
    fout = fopen(argv[4], "wb");
    maxCount = atoi(argv[3]);
    short tmp_wav[FRAME_SIZE];
    fread(tmp_wav,sizeof(char),44,f1);
    fread(tmp_wav,sizeof(char),44,f2);

    // pre skip in noise file
    for (i = 0; i < 150; i++) {
        short tmp[FRAME_SIZE];
        fread(tmp, sizeof(short), FRAME_SIZE, f2);
    }

    // main loop
    while (1) {
        if (count == maxCount) break;
        if ((count % 1000) == 0) fprintf(stderr, "%d\r", count);

        // random gain/filter/freq_range for speech and noise
        if (++gain_change_count > 2821) { // 7 x 13 x 31
            rnnoise_feature_extractor_random_change(st);
            gain_change_count = 0;
        }

        // f1 (speech file)*speech_gain -> x
        // x*x -> E(energy)
        if (rnnoise_get_speech_gain(st) != 0) {

            int len = fread(xi, sizeof(short), FRAME_SIZE, f1);
            if (feof(f1)) {
                rewind(f1);
                int len = fread(xi, sizeof(short), FRAME_SIZE, f1);
            }
        }

        // f2 (noise file)*noise_gain -> n
        if (rnnoise_get_noise_gain(st) != 0) {
            int len = fread(ni, sizeof(short), FRAME_SIZE, f2);
            if (feof(f2)) {
                rewind(f2);
                int len = fread(ni, sizeof(short), FRAME_SIZE, f2);
            }
        }

        rnnoise_feature_extract(st, &xi[0], &ni[0]);
        // feature write out
        count++;
#if 1
        fwrite(rnnoise_get_features(st), sizeof(float), NB_FEATURES, fout);
        fwrite(rnnoise_get_gain(st), sizeof(float), NB_BANDS, fout);
        fwrite(rnnoise_get_Ln(st), sizeof(float), NB_BANDS, fout);
        fwrite(rnnoise_get_vad(st), sizeof(float), 1, fout);
#endif
    }
    fprintf(stderr, "matrix size: %d x %d\n", count, NB_FEATURES + 2 * NB_BANDS + 1);
    fclose(f1);
    fclose(f2);
    fclose(fout);
    return 0;
}

#endif