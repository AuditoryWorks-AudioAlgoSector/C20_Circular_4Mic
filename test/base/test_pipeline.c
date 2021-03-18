#include <sndfile.h>
#include <stdlib.h>
#include "awi_audioframe.h"
#include "awi_filterbankframe.h"

#include "awi_dft_analysis_filter_bank.h"
#include "awi_dft_synthesis_filter_bank.h"
#include "awi_channel_convert.h"
#include "awi_utils.h"
#include "awi_subband_aec_cfg.h"
#include "awi_subband_aec.h"
#include "awi_aes_cfg.h"
#include "awi_aes.h"
#include "awi_fbf_cfg.h"
#include "awi_fbf.h"
#include "awi_abf_cfg.h"
#include "awi_abf.h"
#include "awi_nse_cfg.h"
#include "awi_nse.h"
#include "awi_ns_cfg.h"
#include "awi_ns.h"
#include "awi_agc_cfg.h"
#include "awi_agc.h"
#include "awi_simple_vad_cfg.h"
#include "awi_simple_vad.h"
#include "awi_cng_cfg.h"
#include "awi_cng.h"
#include "awi_limiter.h"
#include "awi_noise_gate.h"
#include "awi_emphasis.h"
#include "awi_beam_fusion.h"

#define NNVAD_USE 0
#define RNNOISE_USE 0
#define ORIG_USE 1
#define ORIG_USE_MOD 0

#if RNNOISE_USE
#include "rnnoise.h"
DenoiseState *rnnoise_st;
static FILE *vad_file = NULL;
#endif

#if  NNVAD_USE
#include "nn_vad_agc.h"
NnvadAgcState* nnvad_st = NULL;
static FILE *nnvad_file = NULL;
#endif

float vad_fullband_snr[AWI_FRAME_PER_BLOACK];
float vad_critical_band_snr[AWI_FRAME_PER_BLOACK][AWI_FREQ_BANDS];

void pipeline_test(char *microphone, char *speaker);
void pipeline_test_mod(char *microphone, char *speaker);
void pipeline_test_rnnoise(char *microphone, char *speaker);
void pipeline_test_nnvad(char *microphone, char *speaker);

int main(int argc, char **argv)
{

#if ORIG_USE
    pipeline_test(argv[1], argv[2]);
#endif

#if RNNOISE_USE
    pipeline_test_rnnoise(argv[1], argv[2]);
#endif

#if NNVAD_USE
    pipeline_test_nnvad(argv[1], argv[2]);
#endif

#if ORIG_USE_MOD
    pipeline_test_mod(argv[1], argv[2]);
#endif
    return 0;
}



#if NNVAD_USE
void pipeline_test_nnvad(char *microphone, char *speaker)
{
    printf("pipeline NNVAD VERSION USE MODEL\n");
    const int frameSize = 256;

    char *mic_file_name = microphone;
    char *spk_file_name = speaker;
    SNDFILE *MicSignalFile, *SpkSignalFile;
    SF_INFO MicFileInfo, SpkFileInfo;

    if(!(MicSignalFile = sf_open(mic_file_name, SFM_READ, &MicFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", mic_file_name);
        exit(1);
    }

    if(!(SpkSignalFile = sf_open(spk_file_name, SFM_READ, &SpkFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", spk_file_name);
        exit(1);
    }

    SNDFILE *fbfOutFile = sf_open("fbf_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aecOutFile = sf_open("aec_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *gscOutFile = sf_open("gsc_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aesOutFile = sf_open("aes_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *bfOutFile  = sf_open("bf_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *cngOutFile = sf_open("cng_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *nsOutFile  = sf_open("ns_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *agcOutFile = sf_open("agc_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *MonoOutFile = sf_open("mono_out_nnvad.wav", SFM_WRITE, &SpkFileInfo); 

    float notch_filter_radius = 0.982f;
    float *mem = (float *) malloc(2 * AWI_MIC_CHANNEL * sizeof(float));
    memset(mem, 0, 2 * AWI_MIC_CHANNEL * sizeof(float));

    awi_fbf_t fbf;
    awi_fbf_cfg_init(&fbf.cfg);
    awi_fbf_init(&fbf);

    awi_subband_aec_part1_t aec_part1;
    awi_subband_aec_part2_t aec_part2;
    awi_subband_aec_cfg_init(&aec_part1.cfg);
    awi_subband_aec_init(&aec_part1, &aec_part2);

    awi_abf_cfg_t abf_cfg;
    awi_abf_cfg_init(&abf_cfg);
    awi_abf_part1_t abf_part1;
    awi_abf_part2_t abf_part2;
    awi_abf_init( &abf_part1, &abf_part2);
    
    awi_aes_t aes;
    awi_aes_cfg_init(&aes.cfg);
    awi_aes_init(&aes);

    awi_nse_t nse_gsc;
    awi_nse_cfg_init(&nse_gsc.cfg, kGlobalWindow120);
    awi_nse_init(&nse_gsc);

    awi_simple_vad_t vad_gsc;
    awi_simple_vad_cfg_init(&vad_gsc.cfg);
    awi_simple_vad_init(&vad_gsc);

    awi_bf_t bf;
    awi_bf_cfg_init(&bf.cfg);
    awi_bf_init(&bf, 0);

    awi_cng_t cng;
    awi_cng_cfg_init(&cng.cfg);
    awi_cng_init(&cng);

    awi_nse_t nse_cng;
    awi_nse_cfg_init(&nse_cng.cfg, kGlobalWindow120);
    awi_nse_init(&nse_cng);

    awi_ns_t ns;
    awi_ns_cfg_init(&ns.cfg);
    awi_ns_init(&ns);

    awi_agc_t agc;
    awi_agc_cfg_init(&agc.cfg);
    awi_agc_init(&agc);

    awi_limiter_t limiter;
    awi_limiter_cfg_init(&limiter.cfg);
    awi_limiter_init(&limiter);

    awi_noise_gate_t noise_gate;
    awi_noise_gate_cfg_init(&noise_gate.cfg);
    awi_noise_gate_init(&noise_gate);

    awi_emphasis_t emphasis;
    awi_emphasis_cfg_init(&emphasis.cfg);
    awi_emphasis_init(&emphasis);

    awi_limiter_t auxi_limiter;
    awi_limiter_cfg_init(&auxi_limiter.cfg);
    auxi_limiter.cfg.max_amplifier   = 30;
    auxi_limiter.cfg.recur_amplifier = 30;
    auxi_limiter.cfg.amp_low_th      = 2;
    auxi_limiter.cfg.amp_high_th     = 30;
    awi_limiter_init(&auxi_limiter);

#if NNVAD_USE
    nnvad_st =  awi_nnvad_agc_create();
    awi_nnvad_agc_init(nnvad_st);
#endif

    audio_mic_frame MicAudioFrame;
    audio_spk_frame SpkAudioFrame;
    
    filter_mic_bank_frame MicFilterBankFrame;
    filter_mic_bank_frame MicFilterBankFrame_Prev;
    filter_mic_bank_frame FbfOutFilterBankFrame;
    filter_mic_bank_frame AecOutFilterBankFrame;
    filter_mic_bank_frame GscOutFilterBankFrame;
    filter_mic_bank_frame AesOutFilterBankFrame;
    filter_spk_bank_frame BfOutFilterBankFrame;
    filter_spk_bank_frame CngOutFilterBankFrame;
    filter_spk_bank_frame NsOutFilterBankFrame;
    filter_spk_bank_frame AgcOutFilterBankFrame;

    filter_spk_bank_frame SpkFilterBankFrame;
    filter_spk_bank_frame SpkFilterBankFrame_Prev;

    dft_analysis_filter_bank MicSignalFilterBank;
    ana_mic_buffer anaMicBuffer;

    dft_analysis_filter_bank SpkSignalFilterBank;
    ana_spk_buffer anaSpkBuffer;

    dft_synthesis_filter_bank FbfOutSynthesisFilterBank;
    dft_synthesis_filter_bank AecOutSynthesisFilterBank;
    dft_synthesis_filter_bank GscOutSynthesisFilterBank;
    dft_synthesis_filter_bank AesOutSynthesisFilterBank;
    dft_synthesis_filter_bank BfOutSynthesisFilterBank;
    dft_synthesis_filter_bank CngOutSynthesisFilterBank;
    dft_synthesis_filter_bank NsOutSynthesisFilterBank;
    dft_synthesis_filter_bank AgcOutSynthesisFilterBank;

    syn_mic_buffer FbfOutSynBuffer;
    syn_mic_buffer AecOutSynBuffer;
    syn_mic_buffer GscOutSynBuffer;
    syn_mic_buffer AesOutSynBuffer;
    syn_spk_buffer BfOutSynBuffer;
    syn_spk_buffer CngOutSynBuffer;
    syn_spk_buffer NsOutSynBuffer;
    syn_spk_buffer AgcOutSynBuffer;

    audio_mic_frame FbfOutAudioFrame;
    audio_mic_frame AecOutAudioFrame;
    audio_mic_frame GscOutAudioFrame;
    audio_mic_frame AesOutAudioFrame;
    audio_spk_frame BfOutAudioFrame;
    audio_spk_frame CngOutAudioFrame;
    audio_spk_frame NsOutAudioFrame;
    audio_spk_frame AgcOutAudioFrame;

    float cng_fd_delay_buffer[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * 2 * NNVAD_DELAY_FRAMES];
    memset(cng_fd_delay_buffer, 0, AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * 2 * NNVAD_DELAY_FRAMES * sizeof(float));

    float bkg_est_psd[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * NNVAD_DELAY_FRAMES];
    memset(bkg_est_psd, 0, AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * NNVAD_DELAY_FRAMES * sizeof(float));

    int aec_on_flag[AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES];
    memset(aec_on_flag, 0, AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES * sizeof(float));

    float bkg_energy[AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES];
    memset(bkg_energy, 0, AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES * sizeof(float));

    awi_initialize_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

    awi_initialize_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank,  BfOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank,  NsOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);


    // printf(" aec_part1:%.2f;\n aec_part2:%.2f;\n fbf:%.2f;\n abf_part1:%.2f;\n abf_part2:%.2f;\n aes:%.2f;\n, cng:%.2f;\n agc:%.2f;\n",
    //         sizeof(awi_subband_aec_part1_t)/1024.f, sizeof(awi_subband_aec_part2_t)/1024.f, sizeof(fbf)/1024.f, sizeof(abf_part1)/1024.f, sizeof(abf_part2)/1024.f, sizeof(aes)/1024.f,
    //        sizeof(cng)/1024.f, sizeof(agc)/1024.f);

    // printf(" nse:%.2f;\n ns:%.2f;\n vad:%.2f;\n bf:%.2f;\n",
    //        sizeof(nse_gsc)/1024.f*2, sizeof(ns)/1024.f, sizeof(vad_gsc)/1024.f, sizeof(bf)/1024.f);


    int block_num = 0;
    int offset = AWI_FRAME_BAND_COUNT * 2;
    memset(SpkFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
    memset(MicFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));
    while(1)
    {
        if(frameSize != sf_readf_float(MicSignalFile, MicAudioFrame, frameSize))
            break;

        if(frameSize != sf_readf_float(SpkSignalFile, SpkAudioFrame, frameSize))
            break;

        awi_inplace_split_channel(SpkAudioFrame, frameSize, AWI_SPK_CHANNEL);

        awi_inplace_split_channel(MicAudioFrame, frameSize, AWI_MIC_CHANNEL);

        for(int i = 0; i < AWI_MIC_CHANNEL; i++)
        {
            awi_filter_dc_notch16(MicAudioFrame + i * frameSize, notch_filter_radius, MicAudioFrame + i * frameSize, frameSize, mem + 2 * i);
        }

        if (emphasis.cfg.enable_flag)
        {
            for(int i = 0; i < AWI_MIC_CHANNEL; i++)
                awi_pre_emphasis_process(&emphasis, MicAudioFrame + i * frameSize, &(emphasis.pre_emp_first_in_chs[i]));

            awi_pre_emphasis_process(&emphasis, SpkAudioFrame, &(emphasis.pre_emp_first_in_spk));
        }

        awi_process_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, MicAudioFrame, MicFilterBankFrame);
        
        awi_process_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, SpkAudioFrame, SpkFilterBankFrame);

        float *mic_sp, *fbf_sp, *aec_sp, *gsc_sp, *aes_sp, *bf_sp, *cng_sp, *ns_sp, *agc_sp;
        float spk_sp[2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2];
        
        if (aec_part1.init_aec_flag)
        {
            memset(spk_sp, 0, 2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
        }

        for(int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;

            mic_sp = MicFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            memcpy(spk_sp, SpkFilterBankFrame + pos * offset, offset * sizeof(float));

            // if (pos == 0)
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame_Prev + ( pos + 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // } else
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame + ( pos - 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // }

            fbf_sp = FbfOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aec_sp = AecOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            gsc_sp = GscOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aes_sp = AesOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            bf_sp  = BfOutFilterBankFrame  + pos * 1 * offset;
            cng_sp = CngOutFilterBankFrame + pos * 1 * offset;

            for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            {
                awi_fbf_process_single(&fbf, mic_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
            }

            awi_subband_aec_process(&aec_part1, &aec_part2, fbf_sp, spk_sp, aec_sp);

            // int refId[AWI_ABF_REFERENCE];
            // for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            // {
            //     for (int ref_id = 0; ref_id < AWI_ABF_REFERENCE; ref_id++ )
            //         refId[ref_id] = fbf.cfg.block_angle[beam_id * AWI_ABF_REFERENCE + ref_id];

            //     float *ref1_sp = aec_sp + refId[0] * offset;
            //     float *ref2_sp = aec_sp + refId[1] * offset;
            //     float *ref3_sp = aec_sp + refId[2] * offset;

            //     awi_abf_process(&abf_part1, &abf_part2, aec_sp + beam_id * offset, ref1_sp, ref2_sp, ref3_sp,
            //                          gsc_sp + beam_id * offset, beam_id);
            // }

            memcpy(gsc_sp, aec_sp, AWI_BEAM_CHANNEL * offset * sizeof(float));

            awi_nse_process_multichannel(&nse_gsc, gsc_sp, 0, 0);

            awi_simple_vad_process(&vad_gsc, nse_gsc.bkg_est_psd, nse_gsc.recur_block_psd, nse_gsc.avg_inst_noisy_psd);

            awi_aes_process(&aes, &aec_part1, fbf_sp, aec_sp, gsc_sp, aes_sp, vad_gsc.bkg_est_psd_floor);

            awi_bf_process(&bf, mic_sp, aes_sp, bf_sp);

            awi_cng_process(&cng, aes.min_aes_gain, bf_sp, cng_sp, 
                            vad_gsc.bkg_est_psd_floor, vad_gsc.recur_fullband_bkg_energy);

            aec_on_flag[AWI_FRAME_PER_BLOACK*(NNVAD_DELAY_FRAMES-1)+pos] = aec_part1.is_aec_on;

            bkg_energy[AWI_FRAME_PER_BLOACK*(NNVAD_DELAY_FRAMES-1)+pos] = vad_gsc.recur_fullband_bkg_energy;

            memcpy(bkg_est_psd + AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * (NNVAD_DELAY_FRAMES-1) + pos * AWI_FRAME_BAND_COUNT, 
                    vad_gsc.bkg_est_psd_floor, AWI_FRAME_BAND_COUNT * sizeof(float));

        }

        memcpy(SpkFilterBankFrame_Prev, SpkFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
        memcpy(MicFilterBankFrame_Prev, MicFilterBankFrame, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));

        awi_process_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AecOutFilterBankFrame, AecOutAudioFrame);
        awi_inplace_merge_channel(AecOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aecOutFile, AecOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, FbfOutFilterBankFrame, FbfOutAudioFrame);
        awi_inplace_merge_channel(FbfOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(fbfOutFile, FbfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, GscOutFilterBankFrame, GscOutAudioFrame);
        awi_inplace_merge_channel(GscOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(gscOutFile, GscOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AesOutFilterBankFrame, AesOutAudioFrame);
        awi_inplace_merge_channel(AesOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aesOutFile, AesOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank, BfOutSynBuffer, BfOutFilterBankFrame, BfOutAudioFrame);
        sf_writef_float(bfOutFile, BfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer, CngOutFilterBankFrame, CngOutAudioFrame);
        awi_limiter_process(&auxi_limiter, CngOutAudioFrame, 1);
        sf_writef_float(cngOutFile, CngOutAudioFrame, frameSize);

        memcpy(cng_fd_delay_buffer + (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * offset, CngOutFilterBankFrame, 
                AWI_FRAME_PER_BLOACK * offset * sizeof(float));
        
        float vad_prob = awi_nnvad_determine_f(nnvad_st, CngOutAudioFrame);

        if ( vad_prob < 0 )
        {
            printf("vad_flag:%.4f\n", vad_prob);
        }
        
        if (nnvad_file == NULL)
        {
            nnvad_file = fopen("nnvad_out.wav", "wb");
            short tmp[256];
            for (size_t i = 0; i < frameSize; i++)
            {
                tmp[i] = vad_prob * 32768.0f;
            }
            fwrite(tmp, sizeof(short), frameSize, nnvad_file);
        }
        else
        {
            short tmp[256];
            for (size_t i = 0; i < frameSize; i++)
            {
                tmp[i] = vad_prob * 32768.0f;
            }
            fwrite(tmp, sizeof(short), frameSize, nnvad_file);
        }
        
        float cng_sp_again[AWI_FRAME_BAND_COUNT * 2];
        int inst_nn_vad, com_vad;
        memset(cng_sp_again, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
        for (int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;
            cng_sp = cng_fd_delay_buffer   + pos * 1 * offset;
            ns_sp  = NsOutFilterBankFrame  + pos * 1 * offset;
            agc_sp = AgcOutFilterBankFrame + pos * 1 * offset;
            inst_nn_vad = vad_prob > 0.5f ? 1 : 0;

            if ( inst_nn_vad )
            {
                vad_gsc.nn_speech_counter++;
                if ( vad_gsc.nn_speech_counter >= vad_gsc.cfg.nn_speech_counter_th )
                {
                    vad_gsc.nn_post_vad_flag = 1;
                    vad_gsc.nn_non_speech_counter = 0;
                }
            }
            else
            {
                vad_gsc.nn_non_speech_counter++;
                if ( vad_gsc.nn_non_speech_counter >= vad_gsc.cfg.nn_non_speech_counter_th )
                {
                    vad_gsc.nn_post_vad_flag = 0;
                    vad_gsc.nn_speech_counter = 0;
                }
            }
            com_vad = vad_gsc.nn_post_vad_flag;
                    
            awi_cng_process_mod(&cng, cng_sp, cng_sp_again, com_vad,
                        bkg_est_psd + pos * AWI_FRAME_BAND_COUNT, bkg_energy[pos]);

            awi_nse_process_mono(&nse_cng, cng_sp_again, com_vad, 1);

            awi_ns_process(&ns, cng_sp_again, ns_sp, nse_cng.avg_inst_noisy_psd, nse_cng.recur_block_psd, nse_cng.bkg_est_psd, cng.pure_noise_flag && com_vad);

            awi_agc_process(&agc, &aes, ns_sp, agc_sp, vad_gsc.fullband_snr, vad_gsc.critical_band_snr,
                            com_vad, cng.pure_noise_flag && com_vad, ns.inst_post_snr_cng);

        }

        memmove(cng_fd_delay_buffer, cng_fd_delay_buffer + AWI_FRAME_PER_BLOACK * offset, 
                (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * offset * sizeof(float));

        memmove(bkg_est_psd, bkg_est_psd + AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT, 
                (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * sizeof(float));

        memmove(aec_on_flag, aec_on_flag + AWI_FRAME_PER_BLOACK, 
                (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * sizeof(float));

        memmove(bkg_energy, bkg_energy + AWI_FRAME_PER_BLOACK, 
                (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * sizeof(float));

        awi_process_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank, NsOutSynBuffer, NsOutFilterBankFrame, NsOutAudioFrame);
        sf_writef_float(nsOutFile, NsOutAudioFrame, frameSize);
        
        awi_process_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer, AgcOutFilterBankFrame, AgcOutAudioFrame);
        sf_writef_float(agcOutFile, AgcOutAudioFrame, frameSize);

        if (emphasis.cfg.enable_flag)
        {
            awi_de_emphasis_process(&emphasis, AgcOutAudioFrame);
        }

        awi_limiter_process(&limiter, AgcOutAudioFrame, cng.pure_noise_flag && vad_gsc.is_speech_triggered);

        awi_noise_gate_process(&noise_gate, AgcOutAudioFrame);

        sf_writef_float(MonoOutFile, AgcOutAudioFrame, frameSize);

    }

    awi_nnvad_agc_free(nnvad_st);

    free(mem);
    sf_close(MicSignalFile);
    sf_close(SpkSignalFile);
    sf_close(aecOutFile);
    sf_close(fbfOutFile);
    sf_close(gscOutFile);
    sf_close(aesOutFile);
    sf_close(bfOutFile);
    sf_close(cngOutFile);
    sf_close(nsOutFile);
    sf_close(agcOutFile);
    sf_close(MonoOutFile);
}
#endif

#if ORIG_USE_MOD
void pipeline_test_mod(char *microphone, char *speaker)
{
    printf("Use modification of original pipeline!\n");
    const int frameSize = 256;

    char *mic_file_name = microphone;
    char *spk_file_name = speaker;
    SNDFILE *MicSignalFile, *SpkSignalFile;
    SF_INFO MicFileInfo, SpkFileInfo;

    if(!(MicSignalFile = sf_open(mic_file_name, SFM_READ, &MicFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", mic_file_name);
        exit(1);
    }

    if(!(SpkSignalFile = sf_open(spk_file_name, SFM_READ, &SpkFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", spk_file_name);
        exit(1);
    }

    SNDFILE *fbfOutFile = sf_open("fbf_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aecOutFile = sf_open("aec_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *gscOutFile = sf_open("gsc_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aesOutFile = sf_open("aes_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *bfOutFile  = sf_open("bf_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *cngOutFile = sf_open("cng_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *nsOutFile  = sf_open("ns_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *agcOutFile = sf_open("agc_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *MonoOutFile = sf_open("mono_out_origin_mod.wav", SFM_WRITE, &SpkFileInfo); 

    float notch_filter_radius = 0.982f;
    float *mem = (float *) malloc(2 * AWI_MIC_CHANNEL * sizeof(float));
    memset(mem, 0, 2 * AWI_MIC_CHANNEL * sizeof(float));

    awi_fbf_t fbf;
    awi_fbf_cfg_init(&fbf.cfg);
    awi_fbf_init(&fbf);

    awi_subband_aec_part1_t aec_part1;
    awi_subband_aec_part2_t aec_part2;
    awi_subband_aec_cfg_init(&aec_part1.cfg);
    awi_subband_aec_init(&aec_part1, &aec_part2);

    awi_abf_cfg_t abf_cfg;
    awi_abf_cfg_init(&abf_cfg);
    awi_abf_part1_t abf_part1;
    awi_abf_part2_t abf_part2;
    awi_abf_init( &abf_part1, &abf_part2);
    
    awi_aes_t aes;
    awi_aes_cfg_init(&aes.cfg);
    awi_aes_init(&aes);

    awi_nse_t nse_gsc;
    awi_nse_cfg_init(&nse_gsc.cfg, kGlobalWindow120);
    awi_nse_init(&nse_gsc);

    awi_simple_vad_t vad_gsc;
    awi_simple_vad_cfg_init(&vad_gsc.cfg);
    awi_simple_vad_init(&vad_gsc);

    awi_bf_t bf;
    awi_bf_cfg_init(&bf.cfg);
    awi_bf_init(&bf, 0);

    awi_cng_t cng;
    awi_cng_cfg_init(&cng.cfg);
    awi_cng_init(&cng);

    awi_nse_t nse_cng;
    awi_nse_cfg_init(&nse_cng.cfg, kGlobalWindow120);
    awi_nse_init(&nse_cng);

    awi_ns_t ns;
    awi_ns_cfg_init(&ns.cfg);
    awi_ns_init(&ns);

    awi_agc_t agc;
    awi_agc_cfg_init(&agc.cfg);
    awi_agc_init(&agc);

    awi_limiter_t limiter;
    awi_limiter_cfg_init(&limiter.cfg);
    awi_limiter_init(&limiter);

    awi_noise_gate_t noise_gate;
    awi_noise_gate_cfg_init(&noise_gate.cfg);
    awi_noise_gate_init(&noise_gate);

    awi_emphasis_t emphasis;
    awi_emphasis_cfg_init(&emphasis.cfg);
    awi_emphasis_init(&emphasis);

    audio_mic_frame MicAudioFrame;
    audio_spk_frame SpkAudioFrame;
    
    filter_mic_bank_frame MicFilterBankFrame;
    filter_mic_bank_frame MicFilterBankFrame_Prev;
    filter_mic_bank_frame FbfOutFilterBankFrame;
    filter_mic_bank_frame AecOutFilterBankFrame;
    filter_mic_bank_frame GscOutFilterBankFrame;
    filter_mic_bank_frame AesOutFilterBankFrame;
    filter_spk_bank_frame BfOutFilterBankFrame;
    filter_spk_bank_frame CngOutFilterBankFrame;
    filter_spk_bank_frame NsOutFilterBankFrame;
    filter_spk_bank_frame AgcOutFilterBankFrame;

    filter_spk_bank_frame SpkFilterBankFrame;
    filter_spk_bank_frame SpkFilterBankFrame_Prev;

    dft_analysis_filter_bank MicSignalFilterBank;
    ana_mic_buffer anaMicBuffer;

    dft_analysis_filter_bank SpkSignalFilterBank;
    ana_spk_buffer anaSpkBuffer;

    dft_synthesis_filter_bank FbfOutSynthesisFilterBank;
    dft_synthesis_filter_bank AecOutSynthesisFilterBank;
    dft_synthesis_filter_bank GscOutSynthesisFilterBank;
    dft_synthesis_filter_bank AesOutSynthesisFilterBank;
    dft_synthesis_filter_bank BfOutSynthesisFilterBank;
    dft_synthesis_filter_bank CngOutSynthesisFilterBank;
    dft_synthesis_filter_bank NsOutSynthesisFilterBank;
    dft_synthesis_filter_bank AgcOutSynthesisFilterBank;

    syn_mic_buffer FbfOutSynBuffer;
    syn_mic_buffer AecOutSynBuffer;
    syn_mic_buffer GscOutSynBuffer;
    syn_mic_buffer AesOutSynBuffer;
    syn_spk_buffer BfOutSynBuffer;
    syn_spk_buffer CngOutSynBuffer;
    syn_spk_buffer NsOutSynBuffer;
    syn_spk_buffer AgcOutSynBuffer;

    audio_mic_frame FbfOutAudioFrame;
    audio_mic_frame AecOutAudioFrame;
    audio_mic_frame GscOutAudioFrame;
    audio_mic_frame AesOutAudioFrame;
    audio_spk_frame BfOutAudioFrame;
    audio_spk_frame CngOutAudioFrame;
    audio_spk_frame NsOutAudioFrame;
    audio_spk_frame AgcOutAudioFrame;

    awi_initialize_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

    awi_initialize_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank,  BfOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank,  NsOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);


    // printf(" aec_part1:%.2f;\n aec_part2:%.2f;\n fbf:%.2f;\n abf_part1:%.2f;\n abf_part2:%.2f;\n aes:%.2f;\n, cng:%.2f;\n agc:%.2f;\n",
    //         sizeof(awi_subband_aec_part1_t)/1024.f, sizeof(awi_subband_aec_part2_t)/1024.f, sizeof(fbf)/1024.f, sizeof(abf_part1)/1024.f, sizeof(abf_part2)/1024.f, sizeof(aes)/1024.f,
    //        sizeof(cng)/1024.f, sizeof(agc)/1024.f);

    // printf(" nse:%.2f;\n ns:%.2f;\n vad:%.2f;\n bf:%.2f;\n",
    //        sizeof(nse_gsc)/1024.f*2, sizeof(ns)/1024.f, sizeof(vad_gsc)/1024.f, sizeof(bf)/1024.f);


    int block_num = 0;
    int offset = AWI_FRAME_BAND_COUNT * 2;
    memset(SpkFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
    memset(MicFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));
    while(1)
    {
        if(frameSize != sf_readf_float(MicSignalFile, MicAudioFrame, frameSize))
            break;

        if(frameSize != sf_readf_float(SpkSignalFile, SpkAudioFrame, frameSize))
            break;

        awi_inplace_split_channel(SpkAudioFrame, frameSize, AWI_SPK_CHANNEL);

        awi_inplace_split_channel(MicAudioFrame, frameSize, AWI_MIC_CHANNEL);

        for(int i = 0; i < AWI_MIC_CHANNEL; i++)
        {
            awi_filter_dc_notch16(MicAudioFrame + i * frameSize, notch_filter_radius, MicAudioFrame + i * frameSize, frameSize, mem + 2 * i);
        }

        if (emphasis.cfg.enable_flag)
        {
            for(int i = 0; i < AWI_MIC_CHANNEL; i++)
                awi_pre_emphasis_process(&emphasis, MicAudioFrame + i * frameSize, &(emphasis.pre_emp_first_in_chs[i]));

            awi_pre_emphasis_process(&emphasis, SpkAudioFrame, &(emphasis.pre_emp_first_in_spk));
        }

        awi_process_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, MicAudioFrame, MicFilterBankFrame);
        
        awi_process_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, SpkAudioFrame, SpkFilterBankFrame);

        float *mic_sp, *fbf_sp, *aec_sp, *gsc_sp, *aes_sp, *bf_sp, *cng_sp, *ns_sp, *agc_sp;
        float spk_sp[2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2];
        float cng_sp_again[AWI_FRAME_BAND_COUNT * 2];
        memset(cng_sp_again, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));

        if (aec_part1.init_aec_flag)
        {
            memset(spk_sp, 0, 2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
        }

        for(int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;

            mic_sp = MicFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            memcpy(spk_sp, SpkFilterBankFrame + pos * offset, offset * sizeof(float));

            // if (pos == 0)
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame_Prev + ( pos + 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // } else
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame + ( pos - 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // }

            fbf_sp = FbfOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aec_sp = AecOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            gsc_sp = GscOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aes_sp = AesOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            bf_sp  = BfOutFilterBankFrame  + pos * 1 * offset;
            cng_sp = CngOutFilterBankFrame + pos * 1 * offset;
            ns_sp  = NsOutFilterBankFrame  + pos * 1 * offset;
            agc_sp = AgcOutFilterBankFrame + pos * 1 * offset;

            for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            {
                awi_fbf_process_single(&fbf, mic_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
            }

            awi_subband_aec_process(&aec_part1, &aec_part2, fbf_sp, spk_sp, aec_sp);

            // int refId[AWI_ABF_REFERENCE];
            // for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            // {
            //     for (int ref_id = 0; ref_id < AWI_ABF_REFERENCE; ref_id++ )
            //         refId[ref_id] = fbf.cfg.block_angle[beam_id * AWI_ABF_REFERENCE + ref_id];

            //     float *ref1_sp = aec_sp + refId[0] * offset;
            //     float *ref2_sp = aec_sp + refId[1] * offset;
            //     float *ref3_sp = aec_sp + refId[2] * offset;

            //     awi_abf_process(&abf_part1, &abf_part2, aec_sp + beam_id * offset, ref1_sp, ref2_sp, ref3_sp,
            //                          gsc_sp + beam_id * offset, beam_id);
            // }

            memcpy(gsc_sp, aec_sp, AWI_BEAM_CHANNEL * offset * sizeof(float));

            awi_nse_process_multichannel(&nse_gsc, gsc_sp, 0, 0);

            awi_simple_vad_process(&vad_gsc, nse_gsc.bkg_est_psd, nse_gsc.recur_block_psd, nse_gsc.avg_inst_noisy_psd);

            awi_aes_process(&aes, &aec_part1, fbf_sp, aec_sp, gsc_sp, aes_sp, vad_gsc.bkg_est_psd_floor);

            awi_bf_process(&bf, mic_sp, aes_sp, bf_sp);

            awi_cng_process(&cng, aes.min_aes_gain, bf_sp, cng_sp, 
                            vad_gsc.bkg_est_psd_floor, vad_gsc.recur_fullband_bkg_energy);

            awi_cng_process_mod(&cng, cng_sp, cng_sp_again, vad_gsc.is_speech_triggered,
                            vad_gsc.bkg_est_psd_floor, vad_gsc.recur_fullband_bkg_energy);

            awi_nse_process_mono(&nse_cng, cng_sp_again, vad_gsc.is_speech_triggered, 1);

            awi_ns_process(&ns, cng_sp_again, ns_sp, nse_cng.avg_inst_noisy_psd, nse_cng.recur_block_psd, nse_cng.bkg_est_psd, cng.pure_noise_flag);

            awi_agc_process(&agc, &aes, ns_sp, agc_sp, vad_gsc.fullband_snr, vad_gsc.critical_band_snr,
                            vad_gsc.is_speech_triggered, cng.pure_noise_flag, ns.inst_post_snr_cng);

        }

        memcpy(SpkFilterBankFrame_Prev, SpkFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
        memcpy(MicFilterBankFrame_Prev, MicFilterBankFrame, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));

        awi_process_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AecOutFilterBankFrame, AecOutAudioFrame);
        awi_inplace_merge_channel(AecOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aecOutFile, AecOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, FbfOutFilterBankFrame, FbfOutAudioFrame);
        awi_inplace_merge_channel(FbfOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(fbfOutFile, FbfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, GscOutFilterBankFrame, GscOutAudioFrame);
        awi_inplace_merge_channel(GscOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(gscOutFile, GscOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AesOutFilterBankFrame, AesOutAudioFrame);
        awi_inplace_merge_channel(AesOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aesOutFile, AesOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank, BfOutSynBuffer, BfOutFilterBankFrame, BfOutAudioFrame);
        sf_writef_float(bfOutFile, BfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer, CngOutFilterBankFrame, CngOutAudioFrame);
        sf_writef_float(cngOutFile, CngOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank, NsOutSynBuffer, NsOutFilterBankFrame, NsOutAudioFrame);
        sf_writef_float(nsOutFile, NsOutAudioFrame, frameSize);
        
        awi_process_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer, AgcOutFilterBankFrame, AgcOutAudioFrame);
        sf_writef_float(agcOutFile, AgcOutAudioFrame, frameSize);

        if (emphasis.cfg.enable_flag)
        {
            awi_de_emphasis_process(&emphasis, AgcOutAudioFrame);
        }

        awi_limiter_process(&limiter, AgcOutAudioFrame, cng.pure_noise_flag && vad_gsc.is_speech_triggered);

        awi_noise_gate_process(&noise_gate, AgcOutAudioFrame);

        sf_writef_float(MonoOutFile, AgcOutAudioFrame, frameSize);

    }

    free(mem);
    sf_close(MicSignalFile);
    sf_close(SpkSignalFile);
    sf_close(aecOutFile);
    sf_close(fbfOutFile);
    sf_close(gscOutFile);
    sf_close(aesOutFile);
    sf_close(bfOutFile);
    sf_close(cngOutFile);
    sf_close(nsOutFile);
    sf_close(agcOutFile);
    sf_close(MonoOutFile);
}
#endif


#if ORIG_USE
void pipeline_test(char *microphone, char *speaker)
{
    const int frameSize = 256;

    char *mic_file_name = microphone;
    char *spk_file_name = speaker;
    SNDFILE *MicSignalFile, *SpkSignalFile;
    SF_INFO MicFileInfo, SpkFileInfo;

    if(!(MicSignalFile = sf_open(mic_file_name, SFM_READ, &MicFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", mic_file_name);
        exit(1);
    }

    if(!(SpkSignalFile = sf_open(spk_file_name, SFM_READ, &SpkFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", spk_file_name);
        exit(1);
    }

    SNDFILE *fbfOutFile = sf_open("fbf_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aecOutFile = sf_open("aec_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *gscOutFile = sf_open("gsc_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aesOutFile = sf_open("aes_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *bfOutFile  = sf_open("bf_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *cngOutFile = sf_open("cng_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *nsOutFile  = sf_open("ns_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *agcOutFile = sf_open("agc_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *MonoOutFile = sf_open("mono_out_orig.wav", SFM_WRITE, &SpkFileInfo); 

    float notch_filter_radius = 0.982f;
    float *mem = (float *) malloc(2 * AWI_MIC_CHANNEL * sizeof(float));
    memset(mem, 0, 2 * AWI_MIC_CHANNEL * sizeof(float));

    awi_fbf_t fbf;
    awi_fbf_cfg_init(&fbf.cfg);
    awi_fbf_init(&fbf);

    awi_subband_aec_part1_t aec_part1;
    awi_subband_aec_part2_t aec_part2;
    awi_subband_aec_cfg_init(&aec_part1.cfg);
    awi_subband_aec_init(&aec_part1, &aec_part2);

    awi_abf_cfg_t abf_cfg;
    awi_abf_cfg_init(&abf_cfg);
    awi_abf_part1_t abf_part1;
    awi_abf_part2_t abf_part2;
    awi_abf_init( &abf_part1, &abf_part2);
    
    awi_aes_t aes;
    awi_aes_cfg_init(&aes.cfg);
    awi_aes_init(&aes);

    awi_nse_t nse_gsc;
    awi_nse_cfg_init(&nse_gsc.cfg, kGlobalWindow120);
    awi_nse_init(&nse_gsc);

    awi_simple_vad_t vad_gsc;
    awi_simple_vad_cfg_init(&vad_gsc.cfg);
    awi_simple_vad_init(&vad_gsc);

    awi_bf_t bf;
    awi_bf_cfg_init(&bf.cfg);
    awi_bf_init(&bf, 0);

    awi_cng_t cng;
    awi_cng_cfg_init(&cng.cfg);
    awi_cng_init(&cng);

    awi_nse_t nse_cng;
    awi_nse_cfg_init(&nse_cng.cfg, kGlobalWindow120);
    awi_nse_init(&nse_cng);

    awi_ns_t ns;
    awi_ns_cfg_init(&ns.cfg);
    awi_ns_init(&ns);

    awi_agc_t agc;
    awi_agc_cfg_init(&agc.cfg);
    awi_agc_init(&agc);

    awi_limiter_t limiter;
    awi_limiter_cfg_init(&limiter.cfg);
    awi_limiter_init(&limiter);

    awi_noise_gate_t noise_gate;
    awi_noise_gate_cfg_init(&noise_gate.cfg);
    awi_noise_gate_init(&noise_gate);

    awi_emphasis_t emphasis;
    awi_emphasis_cfg_init(&emphasis.cfg);
    awi_emphasis_init(&emphasis);

    audio_mic_frame MicAudioFrame;
    audio_spk_frame SpkAudioFrame;
    
    filter_mic_bank_frame MicFilterBankFrame;
    filter_mic_bank_frame MicFilterBankFrame_Prev;
    filter_mic_bank_frame FbfOutFilterBankFrame;
    filter_mic_bank_frame AecOutFilterBankFrame;
    filter_mic_bank_frame GscOutFilterBankFrame;
    filter_mic_bank_frame AesOutFilterBankFrame;
    filter_spk_bank_frame BfOutFilterBankFrame;
    filter_spk_bank_frame CngOutFilterBankFrame;
    filter_spk_bank_frame NsOutFilterBankFrame;
    filter_spk_bank_frame AgcOutFilterBankFrame;

    filter_spk_bank_frame SpkFilterBankFrame;
    filter_spk_bank_frame SpkFilterBankFrame_Prev;

    dft_analysis_filter_bank MicSignalFilterBank;
    ana_mic_buffer anaMicBuffer;

    dft_analysis_filter_bank SpkSignalFilterBank;
    ana_spk_buffer anaSpkBuffer;

    dft_synthesis_filter_bank FbfOutSynthesisFilterBank;
    dft_synthesis_filter_bank AecOutSynthesisFilterBank;
    dft_synthesis_filter_bank GscOutSynthesisFilterBank;
    dft_synthesis_filter_bank AesOutSynthesisFilterBank;
    dft_synthesis_filter_bank BfOutSynthesisFilterBank;
    dft_synthesis_filter_bank CngOutSynthesisFilterBank;
    dft_synthesis_filter_bank NsOutSynthesisFilterBank;
    dft_synthesis_filter_bank AgcOutSynthesisFilterBank;

    syn_mic_buffer FbfOutSynBuffer;
    syn_mic_buffer AecOutSynBuffer;
    syn_mic_buffer GscOutSynBuffer;
    syn_mic_buffer AesOutSynBuffer;
    syn_spk_buffer BfOutSynBuffer;
    syn_spk_buffer CngOutSynBuffer;
    syn_spk_buffer NsOutSynBuffer;
    syn_spk_buffer AgcOutSynBuffer;

    audio_mic_frame FbfOutAudioFrame;
    audio_mic_frame AecOutAudioFrame;
    audio_mic_frame GscOutAudioFrame;
    audio_mic_frame AesOutAudioFrame;
    audio_spk_frame BfOutAudioFrame;
    audio_spk_frame CngOutAudioFrame;
    audio_spk_frame NsOutAudioFrame;
    audio_spk_frame AgcOutAudioFrame;

    awi_initialize_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

    awi_initialize_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank,  BfOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank,  NsOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);


    // printf(" aec_part1:%.2f;\n aec_part2:%.2f;\n fbf:%.2f;\n abf_part1:%.2f;\n abf_part2:%.2f;\n aes:%.2f;\n, cng:%.2f;\n agc:%.2f;\n",
    //         sizeof(awi_subband_aec_part1_t)/1024.f, sizeof(awi_subband_aec_part2_t)/1024.f, sizeof(fbf)/1024.f, sizeof(abf_part1)/1024.f, sizeof(abf_part2)/1024.f, sizeof(aes)/1024.f,
    //        sizeof(cng)/1024.f, sizeof(agc)/1024.f);

    // printf(" nse:%.2f;\n ns:%.2f;\n vad:%.2f;\n bf:%.2f;\n",
    //        sizeof(nse_gsc)/1024.f*2, sizeof(ns)/1024.f, sizeof(vad_gsc)/1024.f, sizeof(bf)/1024.f);


    int block_num = 0;
    int offset = AWI_FRAME_BAND_COUNT * 2;
    memset(SpkFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
    memset(MicFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));
    while(1)
    {
        if(frameSize != sf_readf_float(MicSignalFile, MicAudioFrame, frameSize))
            break;

        if(frameSize != sf_readf_float(SpkSignalFile, SpkAudioFrame, frameSize))
            break;

        awi_inplace_split_channel(SpkAudioFrame, frameSize, AWI_SPK_CHANNEL);

        awi_inplace_split_channel(MicAudioFrame, frameSize, AWI_MIC_CHANNEL);

        for(int i = 0; i < AWI_MIC_CHANNEL; i++)
        {
            awi_filter_dc_notch16(MicAudioFrame + i * frameSize, notch_filter_radius, MicAudioFrame + i * frameSize, frameSize, mem + 2 * i);
        }

        if (emphasis.cfg.enable_flag)
        {
            for(int i = 0; i < AWI_MIC_CHANNEL; i++)
                awi_pre_emphasis_process(&emphasis, MicAudioFrame + i * frameSize, &(emphasis.pre_emp_first_in_chs[i]));

            awi_pre_emphasis_process(&emphasis, SpkAudioFrame, &(emphasis.pre_emp_first_in_spk));
        }

        awi_process_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, MicAudioFrame, MicFilterBankFrame);
        
        awi_process_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, SpkAudioFrame, SpkFilterBankFrame);

        float *mic_sp, *fbf_sp, *aec_sp, *gsc_sp, *aes_sp, *bf_sp, *cng_sp, *ns_sp, *agc_sp;
        float spk_sp[2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2];

        if (aec_part1.init_aec_flag)
        {
            memset(spk_sp, 0, 2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
        }

        for(int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;
            mic_sp = MicFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            memcpy(spk_sp, SpkFilterBankFrame + pos * offset, offset * sizeof(float));

            // if (pos == 0)
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame_Prev + ( pos + 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // } else
            // {
            //     memcpy(spk_sp, SpkFilterBankFrame + ( pos - 1 ) * offset, offset * sizeof(float));
            //     memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            // }

            fbf_sp = FbfOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aec_sp = AecOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            gsc_sp = GscOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aes_sp = AesOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            bf_sp  = BfOutFilterBankFrame  + pos * 1 * offset;
            cng_sp = CngOutFilterBankFrame + pos * 1 * offset;
            ns_sp  = NsOutFilterBankFrame  + pos * 1 * offset;
            agc_sp = AgcOutFilterBankFrame + pos * 1 * offset;

            for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            {
                awi_fbf_process_single(&fbf, mic_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
            }

            awi_subband_aec_process(&aec_part1, &aec_part2, fbf_sp, spk_sp, aec_sp);

            // int refId[AWI_ABF_REFERENCE];
            // for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            // {
            //     for (int ref_id = 0; ref_id < AWI_ABF_REFERENCE; ref_id++ )
            //         refId[ref_id] = fbf.cfg.block_angle[beam_id * AWI_ABF_REFERENCE + ref_id];

            //     float *ref1_sp = aec_sp + refId[0] * offset;
            //     float *ref2_sp = aec_sp + refId[1] * offset;
            //     float *ref3_sp = aec_sp + refId[2] * offset;

            //     awi_abf_process(&abf_part1, &abf_part2, aec_sp + beam_id * offset, ref1_sp, ref2_sp, ref3_sp,
            //                          gsc_sp + beam_id * offset, beam_id);
            // }

            memcpy(gsc_sp, aec_sp, AWI_BEAM_CHANNEL * offset * sizeof(float));

            awi_nse_process_multichannel(&nse_gsc, gsc_sp, 0, 0);

            awi_simple_vad_process(&vad_gsc, nse_gsc.bkg_est_psd, nse_gsc.recur_block_psd, nse_gsc.avg_inst_noisy_psd);

            awi_aes_process(&aes, &aec_part1, fbf_sp, aec_sp, gsc_sp, aes_sp, vad_gsc.bkg_est_psd_floor);

            awi_bf_process(&bf, mic_sp, aes_sp, bf_sp);

            awi_cng_process(&cng, aes.min_aes_gain, bf_sp, cng_sp, 
                            vad_gsc.bkg_est_psd_floor, vad_gsc.recur_fullband_bkg_energy);

            awi_nse_process_mono(&nse_cng, cng_sp, vad_gsc.is_speech_triggered, 1);

            awi_ns_process(&ns, cng_sp, ns_sp, nse_cng.avg_inst_noisy_psd, nse_cng.recur_block_psd, nse_cng.bkg_est_psd, cng.pure_noise_flag);

            awi_agc_process(&agc, &aes, ns_sp, agc_sp, vad_gsc.fullband_snr, vad_gsc.critical_band_snr,
                            vad_gsc.is_speech_triggered, cng.pure_noise_flag, ns.ns_gain_seq);

        }

        memcpy(SpkFilterBankFrame_Prev, SpkFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
        memcpy(MicFilterBankFrame_Prev, MicFilterBankFrame, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));

        awi_process_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AecOutFilterBankFrame, AecOutAudioFrame);
        awi_inplace_merge_channel(AecOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aecOutFile, AecOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, FbfOutFilterBankFrame, FbfOutAudioFrame);
        awi_inplace_merge_channel(FbfOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(fbfOutFile, FbfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, GscOutFilterBankFrame, GscOutAudioFrame);
        awi_inplace_merge_channel(GscOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(gscOutFile, GscOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AesOutFilterBankFrame, AesOutAudioFrame);
        awi_inplace_merge_channel(AesOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aesOutFile, AesOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank, BfOutSynBuffer, BfOutFilterBankFrame, BfOutAudioFrame);
        sf_writef_float(bfOutFile, BfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer, CngOutFilterBankFrame, CngOutAudioFrame);
        sf_writef_float(cngOutFile, CngOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank, NsOutSynBuffer, NsOutFilterBankFrame, NsOutAudioFrame);
        sf_writef_float(nsOutFile, NsOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer, AgcOutFilterBankFrame, AgcOutAudioFrame);
        sf_writef_float(agcOutFile, AgcOutAudioFrame, frameSize);

        if (emphasis.cfg.enable_flag)
        {
            awi_de_emphasis_process(&emphasis, AgcOutAudioFrame);
        }

        awi_limiter_process(&limiter, AgcOutAudioFrame, cng.pure_noise_flag && vad_gsc.is_speech_triggered);

        awi_noise_gate_process(&noise_gate, AgcOutAudioFrame);

        sf_writef_float(MonoOutFile, AgcOutAudioFrame, frameSize);

    }

    free(mem);
    sf_close(MicSignalFile);
    sf_close(SpkSignalFile);
    sf_close(aecOutFile);
    sf_close(fbfOutFile);
    sf_close(gscOutFile);
    sf_close(aesOutFile);
    sf_close(bfOutFile);
    sf_close(cngOutFile);
    sf_close(nsOutFile);
    sf_close(agcOutFile);
    sf_close(MonoOutFile);
}
#endif

#if RNNOISE_USE
void pipeline_test_rnnoise(char *microphone, char *speaker)
{
    printf("RNNOISE VERSION USE MODEL 2020-10-19 58");
    const int frameSize = 256;

    char *mic_file_name = microphone;
    char *spk_file_name = speaker;
    SNDFILE *MicSignalFile, *SpkSignalFile;
    SF_INFO MicFileInfo, SpkFileInfo;

    if(!(MicSignalFile = sf_open(mic_file_name, SFM_READ, &MicFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", mic_file_name);
        exit(1);
    }

    if(!(SpkSignalFile = sf_open(spk_file_name, SFM_READ, &SpkFileInfo)))
    {
        fprintf(stderr, "Invalid file %s\n", spk_file_name);
        exit(1);
    }

    SNDFILE *fbfOutFile = sf_open("fbf_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aecOutFile = sf_open("aec_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *gscOutFile = sf_open("gsc_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *aesOutFile = sf_open("aes_out.wav", SFM_WRITE, &MicFileInfo);
    SNDFILE *bfOutFile  = sf_open("bf_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *cngOutFile = sf_open("cng_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *nsOutFile  = sf_open("ns_out.wav",  SFM_WRITE, &SpkFileInfo);
    SNDFILE *agcOutFile = sf_open("agc_out.wav", SFM_WRITE, &SpkFileInfo);
    SNDFILE *MonoOutFile = sf_open("mono_out_rnnoise.wav", SFM_WRITE, &SpkFileInfo); 

    float notch_filter_radius = 0.982f;
    float *mem = (float *) malloc(2 * AWI_MIC_CHANNEL * sizeof(float));
    memset(mem, 0, 2 * AWI_MIC_CHANNEL * sizeof(float));

    awi_fbf_t fbf;
    awi_fbf_cfg_init(&fbf.cfg);
    awi_fbf_init(&fbf);

    awi_subband_aec_part1_t aec_part1;
    awi_subband_aec_part2_t aec_part2;
    awi_subband_aec_cfg_init(&aec_part1.cfg);
    awi_subband_aec_init(&aec_part1, &aec_part2);

    awi_abf_cfg_t abf_cfg;
    awi_abf_cfg_init(&abf_cfg);
    awi_abf_part1_t abf_part1;
    awi_abf_part2_t abf_part2;
    awi_abf_init( &abf_part1, &abf_part2);
    
    awi_aes_t aes;
    awi_aes_cfg_init(&aes.cfg);
    awi_aes_init(&aes);

    awi_nse_t nse_gsc;
    awi_nse_cfg_init(&nse_gsc.cfg, kGlobalWindow120);
    awi_nse_init(&nse_gsc);

    awi_simple_vad_t vad_gsc;
    awi_simple_vad_cfg_init(&vad_gsc.cfg);
    vad_gsc.cfg.subband_snr_ratio      = 3.0f;
    vad_gsc.cfg.fullband_snr_ratio     = 1.6f;
    vad_gsc.cfg.snr_counter_th         = 10;
    vad_gsc.cfg.speech_counter_th      = 8;
    vad_gsc.cfg.non_speech_counter_th  = 25;
    awi_simple_vad_init(&vad_gsc);

    awi_bf_t bf;
    awi_bf_cfg_init(&bf.cfg);
    awi_bf_init(&bf, 0);

    awi_cng_t cng;
    awi_cng_cfg_init(&cng.cfg);
    awi_cng_init(&cng);

    awi_nse_t nse_cng;
    awi_nse_cfg_init(&nse_cng.cfg, kGlobalWindow120);
    awi_nse_init(&nse_cng);

    awi_ns_t ns;
    awi_ns_cfg_init(&ns.cfg);
    awi_ns_init(&ns);

    awi_agc_t agc;
    awi_agc_cfg_init(&agc.cfg);
    awi_agc_init(&agc);

    awi_limiter_t limiter;
    awi_limiter_cfg_init(&limiter.cfg);
    awi_limiter_init(&limiter);

    awi_noise_gate_t noise_gate;
    awi_noise_gate_cfg_init(&noise_gate.cfg);
    awi_noise_gate_init(&noise_gate);

    awi_emphasis_t emphasis;
    awi_emphasis_cfg_init(&emphasis.cfg);
    awi_emphasis_init(&emphasis);

#if RNNOISE_USE
    rnnoise_st =  rnnoise_create(NULL);
#endif

    audio_mic_frame MicAudioFrame;
    audio_spk_frame SpkAudioFrame;
    
    filter_mic_bank_frame MicFilterBankFrame;
    filter_mic_bank_frame MicFilterBankFrame_Prev;
    filter_mic_bank_frame FbfOutFilterBankFrame;
    filter_mic_bank_frame AecOutFilterBankFrame;
    filter_mic_bank_frame GscOutFilterBankFrame;
    filter_mic_bank_frame AesOutFilterBankFrame;
    filter_spk_bank_frame BfOutFilterBankFrame;
    filter_spk_bank_frame CngOutFilterBankFrame;
    filter_spk_bank_frame NsOutFilterBankFrame;
    filter_spk_bank_frame AgcOutFilterBankFrame;

    filter_spk_bank_frame SpkFilterBankFrame;
    filter_spk_bank_frame SpkFilterBankFrame_Prev;

    dft_analysis_filter_bank MicSignalFilterBank;
    ana_mic_buffer anaMicBuffer;

    dft_analysis_filter_bank SpkSignalFilterBank;
    ana_spk_buffer anaSpkBuffer;

    dft_synthesis_filter_bank FbfOutSynthesisFilterBank;
    dft_synthesis_filter_bank AecOutSynthesisFilterBank;
    dft_synthesis_filter_bank GscOutSynthesisFilterBank;
    dft_synthesis_filter_bank AesOutSynthesisFilterBank;
    dft_synthesis_filter_bank BfOutSynthesisFilterBank;
    dft_synthesis_filter_bank CngOutSynthesisFilterBank;
    dft_synthesis_filter_bank NsOutSynthesisFilterBank;
    dft_synthesis_filter_bank AgcOutSynthesisFilterBank;

    syn_mic_buffer FbfOutSynBuffer;
    syn_mic_buffer AecOutSynBuffer;
    syn_mic_buffer GscOutSynBuffer;
    syn_mic_buffer AesOutSynBuffer;
    syn_spk_buffer BfOutSynBuffer;
    syn_spk_buffer CngOutSynBuffer;
    syn_spk_buffer NsOutSynBuffer;
    syn_spk_buffer AgcOutSynBuffer;

    audio_mic_frame FbfOutAudioFrame;
    audio_mic_frame AecOutAudioFrame;
    audio_mic_frame GscOutAudioFrame;
    audio_mic_frame AesOutAudioFrame;
    audio_spk_frame BfOutAudioFrame;
    audio_spk_frame CngOutAudioFrame;
    audio_spk_frame NsOutAudioFrame;
    audio_spk_frame AgcOutAudioFrame;

    awi_initialize_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

    awi_initialize_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank,  BfOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank,  NsOutSynBuffer,    AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer,   AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);


    printf(" aec_part1:%.2f;\n aec_part2:%.2f;\n fbf:%.2f;\n abf_part1:%.2f;\n abf_part2:%.2f;\n aes:%.2f;\n, cng:%.2f;\n agc:%.2f;\n",
            sizeof(awi_subband_aec_part1_t)/1024.f, sizeof(awi_subband_aec_part2_t)/1024.f, sizeof(fbf)/1024.f, sizeof(abf_part1)/1024.f, sizeof(abf_part2)/1024.f, sizeof(aes)/1024.f,
           sizeof(cng)/1024.f, sizeof(agc)/1024.f);

    printf(" nse:%.2f;\n ns:%.2f;\n vad:%.2f;\n bf:%.2f;\n",
           sizeof(nse_gsc)/1024.f*2, sizeof(ns)/1024.f, sizeof(vad_gsc)/1024.f, sizeof(bf)/1024.f);


    int block_num = 0;
    int offset = AWI_FRAME_BAND_COUNT * 2;
    memset(SpkFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
    memset(MicFilterBankFrame_Prev, 0, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));
    while(1)
    {
        if(frameSize != sf_readf_float(MicSignalFile, MicAudioFrame, frameSize))
            break;

        if(frameSize != sf_readf_float(SpkSignalFile, SpkAudioFrame, frameSize))
            break;

        awi_inplace_split_channel(SpkAudioFrame, frameSize, AWI_SPK_CHANNEL);

        awi_inplace_split_channel(MicAudioFrame, frameSize, AWI_MIC_CHANNEL);

        for(int i = 0; i < AWI_MIC_CHANNEL; i++)
        {
            awi_filter_dc_notch16(MicAudioFrame + i * frameSize, notch_filter_radius, MicAudioFrame + i * frameSize, frameSize, mem + 2 * i);
        }

        if (emphasis.cfg.enable_flag)
        {
            for(int i = 0; i < AWI_MIC_CHANNEL; i++)
                awi_pre_emphasis_process(&emphasis, MicAudioFrame + i * frameSize, &(emphasis.pre_emp_first_in_chs[i]));

            awi_pre_emphasis_process(&emphasis, SpkAudioFrame, &(emphasis.pre_emp_first_in_spk));
        }

        awi_process_dft_analysis_filter_bank(&MicSignalFilterBank, anaMicBuffer, MicAudioFrame, MicFilterBankFrame);
        
        awi_process_dft_analysis_filter_bank(&SpkSignalFilterBank, anaSpkBuffer, SpkAudioFrame, SpkFilterBankFrame);

        float *mic_sp, *fbf_sp, *aec_sp, *gsc_sp, *aes_sp, *bf_sp, *cng_sp, *ns_sp, *agc_sp;
        float spk_sp[2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2];

        if (aec_part1.init_aec_flag)
        {
            memset(spk_sp, 0, 2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
        }

        for(int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;
            mic_sp = MicFilterBankFrame_Prev + pos * AWI_MIC_CHANNEL * offset;
            // memcpy(spk_sp, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));

            if (pos == 0)
            {
                memcpy(spk_sp, SpkFilterBankFrame_Prev + ( pos + 1 ) * offset, offset * sizeof(float));
                memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            } else
            {
                memcpy(spk_sp, SpkFilterBankFrame + ( pos - 1 ) * offset, offset * sizeof(float));
                memcpy(spk_sp + 1 * offset, SpkFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
            }

            fbf_sp = FbfOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aec_sp = AecOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            gsc_sp = GscOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            aes_sp = AesOutFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
            bf_sp  = BfOutFilterBankFrame  + pos * 1 * offset;
            cng_sp = CngOutFilterBankFrame + pos * 1 * offset;
            ns_sp  = NsOutFilterBankFrame  + pos * 1 * offset;
            agc_sp = AgcOutFilterBankFrame + pos * 1 * offset;

            for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            {
                awi_fbf_process_single(&fbf, mic_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
            }

            awi_subband_aec_process(&aec_part1, &aec_part2, fbf_sp, spk_sp, aec_sp);

            // int refId[AWI_ABF_REFERENCE];
            // for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
            // {
            //     for (int ref_id = 0; ref_id < AWI_ABF_REFERENCE; ref_id++ )
            //         refId[ref_id] = fbf.cfg.block_angle[beam_id * AWI_ABF_REFERENCE + ref_id];

            //     float *ref1_sp = aec_sp + refId[0] * offset;
            //     float *ref2_sp = aec_sp + refId[1] * offset;
            //     float *ref3_sp = aec_sp + refId[2] * offset;

            //     awi_abf_process(&abf_part1, &abf_part2, aec_sp + beam_id * offset, ref1_sp, ref2_sp, ref3_sp,
            //                          gsc_sp + beam_id * offset, beam_id);
            // }

            memcpy(gsc_sp, aec_sp, AWI_BEAM_CHANNEL * offset * sizeof(float));

            awi_nse_process_multichannel(&nse_gsc, gsc_sp, 0, 0);

            awi_simple_vad_process(&vad_gsc, nse_gsc.bkg_est_psd, nse_gsc.recur_block_psd, nse_gsc.avg_inst_noisy_psd);

            //for test
            vad_fullband_snr[pos] = vad_gsc.fullband_snr;
            memcpy(vad_critical_band_snr[pos], vad_gsc.critical_band_snr, sizeof(float)*AWI_FREQ_BANDS);

            awi_aes_process(&aes, &aec_part1, fbf_sp, aec_sp, gsc_sp, aes_sp, vad_gsc.bkg_est_psd_floor);

            awi_bf_process(&bf, mic_sp, aes_sp, bf_sp);

            awi_cng_process(&cng, aes.min_aes_gain, bf_sp, cng_sp, 
                            vad_gsc.bkg_est_psd_floor, vad_gsc.recur_fullband_bkg_energy);

            awi_nse_process_mono(&nse_cng, cng_sp, vad_gsc.is_speech_triggered, 1);

            awi_ns_process(&ns, cng_sp, ns_sp, nse_cng.avg_inst_noisy_psd, nse_cng.recur_block_psd, nse_cng.bkg_est_psd, cng.pure_noise_flag);

            //awi_agc_process(&agc, &aes, ns_sp, agc_sp, vad_gsc.fullband_snr, vad_gsc.critical_band_snr,
            //                vad_gsc.is_speech_triggered, cng.pure_noise_flag, ns.inst_post_snr_cng);

        }

        memcpy(SpkFilterBankFrame_Prev, SpkFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));
        memcpy(MicFilterBankFrame_Prev, MicFilterBankFrame, AWI_FRAME_PER_BLOACK * AWI_MIC_CHANNEL * offset * sizeof(float));

        awi_process_dft_synthesis_filter_bank(&AecOutSynthesisFilterBank, AecOutSynBuffer, AecOutFilterBankFrame, AecOutAudioFrame);
        awi_inplace_merge_channel(AecOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aecOutFile, AecOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&FbfOutSynthesisFilterBank, FbfOutSynBuffer, FbfOutFilterBankFrame, FbfOutAudioFrame);
        awi_inplace_merge_channel(FbfOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(fbfOutFile, FbfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&GscOutSynthesisFilterBank, GscOutSynBuffer, GscOutFilterBankFrame, GscOutAudioFrame);
        awi_inplace_merge_channel(GscOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(gscOutFile, GscOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&AesOutSynthesisFilterBank, AesOutSynBuffer, AesOutFilterBankFrame, AesOutAudioFrame);
        awi_inplace_merge_channel(AesOutAudioFrame, AWI_BAND_COUNT, AWI_MIC_CHANNEL);
        sf_writef_float(aesOutFile, AesOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&BfOutSynthesisFilterBank, BfOutSynBuffer, BfOutFilterBankFrame, BfOutAudioFrame);
        sf_writef_float(bfOutFile, BfOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&CngOutSynthesisFilterBank, CngOutSynBuffer, CngOutFilterBankFrame, CngOutAudioFrame);
        sf_writef_float(cngOutFile, CngOutAudioFrame, frameSize);

        awi_process_dft_synthesis_filter_bank(&NsOutSynthesisFilterBank, NsOutSynBuffer, NsOutFilterBankFrame, NsOutAudioFrame);
        sf_writef_float(nsOutFile, NsOutAudioFrame, frameSize);

        for (size_t i = 0; i < frameSize; i++)
        {
            NsOutAudioFrame[i] = NsOutAudioFrame[i]* 32767 *  16;
        }
    
        float vad_prob = rnnoise_process_frame_vad(rnnoise_st, NsOutAudioFrame);
        
        if (vad_file == NULL)
        {
            vad_file = fopen("vad_out.wav", "wb");
        }else{
            short tmp[256];
            for (size_t i = 0; i < frameSize; i++)
            {
                tmp[i] = vad_prob * 16384.0f;
            }
            int len = fwrite(tmp, sizeof(short), frameSize, vad_file);
            //printf("len %d, vad prob %f\n", len, vad_prob);
        }
        
        for (int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
        {
            block_num++;
            ns_sp  = NsOutFilterBankFrame  + pos * 1 * offset;
            agc_sp = AgcOutFilterBankFrame + pos * 1 * offset;
            vad_gsc.is_speech_triggered = vad_prob > 0.90f?1:0;
            cng.pure_noise_flag = vad_prob > 0.9?1:0;
            awi_agc_process(&agc, &aes, ns_sp, agc_sp, vad_fullband_snr[pos], vad_critical_band_snr[pos],
                            vad_gsc.is_speech_triggered, cng.pure_noise_flag, ns.inst_post_snr_cng);
        }
        
        awi_process_dft_synthesis_filter_bank(&AgcOutSynthesisFilterBank, AgcOutSynBuffer, AgcOutFilterBankFrame, AgcOutAudioFrame);
        sf_writef_float(agcOutFile, AgcOutAudioFrame, frameSize);

        if (emphasis.cfg.enable_flag)
        {
            awi_de_emphasis_process(&emphasis, AgcOutAudioFrame);
        }

        awi_limiter_process(&limiter, AgcOutAudioFrame, cng.pure_noise_flag && vad_gsc.is_speech_triggered);

        awi_noise_gate_process(&noise_gate, AgcOutAudioFrame);

        sf_writef_float(MonoOutFile, AgcOutAudioFrame, frameSize);

    }

    rnnoise_destroy(rnnoise_st);

    free(mem);
    sf_close(MicSignalFile);
    sf_close(SpkSignalFile);
    sf_close(aecOutFile);
    sf_close(fbfOutFile);
    sf_close(gscOutFile);
    sf_close(aesOutFile);
    sf_close(bfOutFile);
    sf_close(cngOutFile);
    sf_close(nsOutFile);
    sf_close(agcOutFile);
    sf_close(MonoOutFile);
}
#endif
