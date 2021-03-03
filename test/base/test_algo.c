//
// Created by xzc on 19-11-24.
//
#include <sndfile.h>
#include <stdlib.h>
#include "awi_algo.h"
#include "awi_algo_module_accu.h"
#include "awi_audioframe.h"
#include "awi_filterbankframe.h"

void test_algo(char *test_audio_file);


int main(int argc, char **argv)
{
    test_algo(argv[1]);

    return 0;
}


void test_algo(char *test_audio_file)
{
    const sf_count_t frameSize = AWI_FRAME_LENGTH;
    const int total_chs_size = frameSize*(AWI_MIC_CHANNEL+AWI_SPK_CHANNEL+1);

    float float_audio_input[total_chs_size];
    short short_audio_input[total_chs_size];
    short short_audio_output[frameSize];

    SNDFILE *test_wav_file;
    SF_INFO TestFileInfo;

    if(!(test_wav_file = sf_open(test_audio_file, SFM_READ, &TestFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", test_audio_file);
        exit(1);
    }

    SNDFILE *recordOutFile = sf_open("record_out.wav", SFM_WRITE, &TestFileInfo);

    awi_algo_modules_accu_t *test_algo;
    test_algo = (awi_algo_modules_accu_t*)awi_algo_init();

    int block_num = 0;
    sf_count_t read_size = 0;
    while(1)
    {
        read_size = sf_readf_float(test_wav_file, float_audio_input, frameSize);
        if (read_size != frameSize)
            break;

        for(int i = 0; i < total_chs_size; i++)
        {
            if ( float_audio_input[i] < 0 )
                short_audio_input[i] = (short)(float_audio_input[i] * 32768 - 0.5f);
            else
                short_audio_input[i] = (short)(float_audio_input[i] * 32767 + 0.5f);

        }

        awi_algo_process(test_algo, short_audio_input, short_audio_output);
        block_num++;

        sf_writef_short(recordOutFile, short_audio_input, frameSize);

    }

    awi_algo_deinit(test_algo);

    sf_close(test_wav_file);

    sf_close(recordOutFile);

}
