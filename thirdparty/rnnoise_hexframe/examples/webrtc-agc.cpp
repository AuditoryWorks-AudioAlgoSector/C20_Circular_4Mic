// webrtc.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <string>
#include "modules/audio_processing/agc/legacy/gain_control.h"
#include "common_audio/resampler/include/resampler.h"
using namespace webrtc;
using namespace std;

#define IN_FILE "./ns_out.wav"
#define OUT_FILE "./ns_out_webrtc.wav"
#define FRAME_SIZE 160

void audioResample(const char* path, const char* out_path, int in_freq, int out_freq) {
	webrtc::Resampler resampler;
	resampler.Reset(in_freq, out_freq, 1);
	FILE* fin = fopen(path, "rb");
	FILE* fout = fopen(out_path, "wb");
	int16_t inbuf[160];
	int16_t outbuf[480];
	size_t outlen = 0;
    fread(inbuf, sizeof(char), 44, fin);
	while (!feof(fin))
	{
		int len = fread(inbuf, sizeof(int16_t), 160, fin);
		//printf("read len = %d\n", len);
		resampler.Push(inbuf, len, outbuf, 480, outlen);
		//printf("resample out len = %d\n", outlen);
		fwrite(outbuf, sizeof(int16_t), outlen, fout);
	}

    printf("resample hello world\n");

	fclose(fin);
	fclose(fout);
}

int main(int argc, char** argv){
    audioResample(IN_FILE, "./resample_out.pcm", 160, 480);
    return 0;
}

int main1(int argc, char** argv)
{
    FILE* input = fopen(IN_FILE, "rb");
    FILE* output = fopen(OUT_FILE, "wb");
    short tmp[1][FRAME_SIZE], out[1][FRAME_SIZE];
    int32_t gains[11] = { 0 };
    std::cout << "Hello World!\n";
    void* agcInst = WebRtcAgc_Create();
    WebRtcAgc_Init(agcInst, 0, 255, 2, 16000);


    WebRtcAgcConfig agcConfig;
    agcConfig.compressionGaindB = 40;//在Fixed模式下，越大声音越大
    agcConfig.limiterEnable = 1;
    agcConfig.targetLevelDbfs = 6;   //dbfs表示相对于full scale的下降值，0表示full scale，越小声音越大
    WebRtcAgc_set_config(agcInst, agcConfig);

    int32_t capture_levels = 0;
    int32_t new_capture_level = 0;
    uint8_t saturation_warning = 0;
    bool stream_has_echo = false;
    while (!feof(input))
    {
        bool stream_is_saturated_ = false;
        bool error_reported = false;

        int16_t split_band_data[3][FRAME_SIZE];
        int16_t* split_bands[3] = {split_band_data[0], split_band_data[1], split_band_data[2] };

        int len = fread(split_bands[0], sizeof(short), FRAME_SIZE, input);
        if (len < FRAME_SIZE)
        {
            break;
        }
        for (size_t i = 0; i < len; i++)
        {
            /* code */
            split_bands[0][i]*=16;
        }
        
        // The call to stream_has_echo() is ok from a deadlock perspective
        // as the capture lock is allready held.

        int err_analyze = WebRtcAgc_Analyze(
                agcInst, split_bands, 1,
                FRAME_SIZE, capture_levels, &new_capture_level,
                stream_has_echo, &saturation_warning, gains);
        capture_levels = new_capture_level;

        WebRtcAgc_Process(agcInst, gains, (const int16_t* const*)split_bands, 1, (int16_t* const*)split_bands);
        fwrite(split_bands[0], sizeof(short), FRAME_SIZE, output);
    }
    std::cout << "Hello World!\n";
    return 1;

}
