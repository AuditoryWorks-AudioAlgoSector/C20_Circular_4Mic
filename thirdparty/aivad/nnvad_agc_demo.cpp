#include <fstream>
#include <sys/time.h>
#include "nn_vad_agc.h"
int vad_origin(int argc, char **argv);
int vad_split(int argc, char **argv);
int main(int argc, char **argv){
    //vad_origin(argc,argv);
    vad_split(argc,argv);
    return 0;
}
int vad_origin(int argc, char **argv) {
    //if (argc < 2) return 0;
    FILE *finput = fopen("./ns_out.wav", "rb");//需要拿去做agc的数据
    FILE *forigin = fopen("./z-02.wav", "rb");//需要拿去做vad的数据
    

    FILE *fvadout1 = fopen("./determine_out.wav", "wb");
    //FILE *fvadout2 = fopen("./origin_out.wav", "rb");

    FILE *foutput = fopen("./ragc_out.wav", "wb");
    NnvadAgcState *st = awi_nnvad_agc_create();
    awi_nnvad_agc_init(st);

    struct timeval t1, t2;
    int count = 0;
    short intmp[256];
    short tmpvad[256];
    float infloat[256];
    float vadfloat[256];
    short outtmp[256];
    fread(intmp,sizeof(char),44,finput);
    fread(tmpvad, sizeof(char), 44, forigin);
    while (!feof(finput)) {
        //int len = fread(infloat, sizeof(float), 256, finput);
        int len = fread(intmp, sizeof(short ), 256, finput);
        int vadlen = fread(tmpvad, sizeof(short), 256, forigin);
        count++;
        if (len < 256) {
            break;
        }
        for (int i = 0; i < 256; ++i) {
            infloat[i] = intmp[i] / 32767.0f;
            vadfloat[i] = tmpvad[i] / 32767.0f;
        }

//         int vad_result = awi_nnvad_determine(st, vadfloat);
// {
//     short vad_tmp[256] = {0};
//     for (size_t i = 0; i < 256; i++)
//     {
//         vad_tmp[i] = vad_result * 20000;
//     }
//     fwrite(vad_tmp, sizeof(short), 256, fvadout1);
// }

        gettimeofday(&t1, nullptr);

        awi_nnvad_agc_process(st, infloat, outtmp);

        gettimeofday(&t2, nullptr);
    long deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    //long deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
    if (deltaT1 >= 5000) {
    //if(1){
        
        printf("\n\rframe num ，RNNoise process one frame(16ms) cost %ld us, count : %d\r\n", deltaT1, count);
    }
        

        fwrite(outtmp, sizeof(short), 256, foutput);
    }

    awi_nnvad_agc_free(st);
    printf("//hello world stop exec!");
    return 0;
}
int vad_split(int argc, char **argv) {
    //if (argc < 2) return 0;
    FILE *finput = fopen("./ns_out.wav", "rb");//需要拿去做agc的数据
    FILE *forigin = fopen("./c-02.wav", "rb");//需要拿去做vad的数据
    

    FILE *fvadout1 = fopen("./determine_out.wav", "wb");
    //FILE *fvadout2 = fopen("./origin_out.wav", "rb");

    FILE *foutput = fopen("./ragc_out.wav", "wb");
    NnvadAgcState *st = awi_nnvad_agc_create();
    awi_nnvad_agc_init(st);

    struct timeval t1, t2;
    int count = 0;
    short intmp[256];
    short tmpvad[256];
    float infloat[256];
    float vadfloat[256];
    short outtmp[256];
    fread(intmp,sizeof(char),44,finput);
    fread(tmpvad, sizeof(char), 44, forigin);
    while (!feof(finput)) {
        //int len = fread(infloat, sizeof(float), 256, finput);
        int len = fread(intmp, sizeof(short ), 256, finput);
        int vadlen = fread(tmpvad, sizeof(short), 256, forigin);
        count++;
        if (len < 256) {
            break;
        }
        for (int i = 0; i < 256; ++i) {
            infloat[i] = intmp[i] / 32767.0f;
            vadfloat[i] = tmpvad[i] / 32767.0f;
        }

        int vad_result = awi_nnvad_determine(st, infloat);
{
    short vad_tmp[256] = {0};
    for (size_t i = 0; i < 256; i++)
    {
        vad_tmp[i] = vad_result * 20000;
    }
    fwrite(vad_tmp, sizeof(short), 256, fvadout1);
}

        awi_nnvad_access(st, infloat, outtmp);    

        gettimeofday(&t1, nullptr);

        awi_nnvad_agc_process(st, infloat, outtmp);

        gettimeofday(&t2, nullptr);
    long deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    //long deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
    if (deltaT1 >= 5000) {
    //if(1){
        
        printf("\n\rframe num ，RNNoise process one frame(16ms) cost %ld us, count : %d\r\n", deltaT1, count);
    }
        

        fwrite(outtmp, sizeof(short), 256, foutput);
    }

    awi_nnvad_agc_free(st);
    printf("//hello world stop exec!");
    return 0;
}