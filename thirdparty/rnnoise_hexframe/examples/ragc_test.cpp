
#include <fstream>
#include "rnnoise_agc.h"
#include <sys/time.h>


int main(int argc, char **argv) {
    if (argc < 2) return 0;
    FILE *finput = fopen(argv[1], "rb");
    FILE *foutput = fopen("ragc_out.pcm", "wb");
    RnnoiseAgcState *st = awi_rnnoise_agc_create();
    awi_rnnoise_agc_init(st);

    struct timeval t1, t2;

    short intmp[256];
    float infloat[256];
    short outtmp[256];
    fread(intmp,sizeof(char),44,finput);
    while (!feof(finput)) {
        //int len = fread(infloat, sizeof(float), 256, finput);
        int len = fread(intmp, sizeof(short ), 256, finput);
        if (len < 256) {
            break;
        }
        for (int i = 0; i < 256; ++i) {
            infloat[i] = intmp[i] * 1.0f;
        }
        gettimeofday(&t1, nullptr);
        awi_rnnoise_agc_process(st, infloat, outtmp);
        gettimeofday(&t2, nullptr);
        long deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        //long deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
        if (deltaT1 >= 500) {
            printf("\rframe num ，RNNoise process one frame(10ms) cost %ld us\r", deltaT1);
        }
        fwrite(outtmp, sizeof(short), 256, foutput);
    }
    awi_rnnoise_agc_free(st);
    printf("hello world stop exec!");
    return 0;
}