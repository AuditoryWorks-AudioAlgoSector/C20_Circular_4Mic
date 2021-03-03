#include "awi_algo.h"
#include "awi_algo_module_accu.h"
#include "awi_utils.h"
#include "awi_channel_convert.h"

#define DEBUG_VAD 0

#if DEBUG_VAD
#include <sndfile.h>
static FILE *nnvad_file = NULL;
#endif

void *awi_algo_init()
{

    awi_algo_modules_accu_t *st = (awi_algo_modules_accu_t *)malloc(sizeof(awi_algo_modules_accu_t));

    if (NULL == st)
    {
        printf("Malloc Error in function awi_algo_init()!\n\n");
        return NULL;
    }

    memset(st, 0, sizeof(awi_algo_modules_accu_t));
 
    memset(st->notch_mem, 0, 2 * AWI_MIC_CHANNEL * sizeof(float));

#ifdef WPE_USE
    printf("eigen wpe use!!!\n");
    awi_wpe_cfg_init(&st->wpe_cfg, AWI_MIC_CHANNEL);
    st->wpe = awi_eigen_wpe_create();
    awi_eigen_wpe_init(st->wpe, &st->wpe_cfg);
#endif 

    awi_fbf_cfg_init(&st->fbf.cfg);
    awi_fbf_init(&st->fbf);

    awi_subband_aec_cfg_init(&st->aec_part1.cfg);
    awi_subband_aec_init(&st->aec_part1, &st->aec_part2);

    awi_abf_cfg_init(&st->abf_part1.cfg);
    awi_abf_init( &st->abf_part1, &st->abf_part2);

    awi_aes_cfg_init(&st->aes.cfg);
    awi_aes_init(&st->aes);

    awi_nse_cfg_init(&st->nse_gsc.cfg, kGlobalWindow120);
    awi_nse_init(&st->nse_gsc);

    awi_simple_vad_cfg_init(&st->vad_gsc.cfg);
    awi_simple_vad_init(&st->vad_gsc);

    awi_bf_cfg_init(&st->bf.cfg);
    awi_bf_init(&st->bf, 0);

    awi_cng_cfg_init(&st->cng.cfg);
    awi_cng_init(&st->cng);

    awi_nse_cfg_init(&st->nse_cng.cfg, kGlobalWindow120);
    awi_nse_init(&st->nse_cng);

    awi_ns_cfg_init(&st->ns.cfg);
    awi_ns_init(&st->ns);

    awi_agc_cfg_init(&st->agc.cfg);
    awi_agc_init(&st->agc);

    awi_limiter_cfg_init(&st->voip_limiter.cfg);
    awi_limiter_init(&st->voip_limiter);

    awi_noise_gate_cfg_init(&st->noise_gate.cfg);
    awi_noise_gate_init(&st->noise_gate);

    awi_emphasis_cfg_init(&st->emphasis.cfg);
    awi_emphasis_init(&st->emphasis);

    awi_limiter_cfg_init(&st->auxi_limiter.cfg);
    st->auxi_limiter.cfg.max_amplifier   = 30;
    st->auxi_limiter.cfg.recur_amplifier = 30;
    st->auxi_limiter.cfg.amp_low_th      = 2;
    st->auxi_limiter.cfg.amp_high_th     = 30;
    awi_limiter_init(&st->auxi_limiter);

    awi_initialize_dft_analysis_filter_bank(&st->MicSignalFilterBank, st->anaMicBuffer, AWI_MIC_CHANNEL, AWI_FRAME_BAND_COUNT);
    awi_initialize_dft_analysis_filter_bank(&st->SpkSignalFilterBank, st->anaSpkBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

    awi_initialize_dft_synthesis_filter_bank(&st->voipDftSynthesisFilterBank, st->voipSynBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);

#ifdef NNVAD_USE
    awi_initialize_dft_synthesis_filter_bank(&st->CngOutSynthesisFilterBank, st->CngOutSynBuffer, AWI_SPK_CHANNEL, AWI_FRAME_BAND_COUNT);
    st->nnvad_st =  awi_nnvad_agc_create();
    awi_nnvad_agc_init(st->nnvad_st);

    memset(st->cng_fd_delay_buffer, 0, AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * 2 * NNVAD_DELAY_FRAMES * sizeof(float));

    memset(st->bkg_est_psd, 0, AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * NNVAD_DELAY_FRAMES * sizeof(float));

    memset(st->aec_on_flag, 0, AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES * sizeof(float));

    memset(st->bkg_energy, 0, AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES * sizeof(float));
#endif

    return st;

}


void awi_algo_process(void *algo, short *multi_chs_in, short *output)
{
    awi_algo_modules_accu_t *st = (awi_algo_modules_accu_t *)algo;

    awi_normal_split_channel(st->audioFrame, multi_chs_in, AWI_FRAME_LENGTH, AWI_MIC_CHANNEL + AWI_SPK_CHANNEL + 1);

    for(int i = 1; i < AWI_MIC_CHANNEL + 1; i++)
    {
        awi_filter_dc_notch16(st->audioFrame + i * AWI_FRAME_LENGTH, 0.982f, st->audioFrame + i * AWI_FRAME_LENGTH, AWI_FRAME_LENGTH, st->notch_mem + 2 * ( i - 1 ));
    }

    if (st->emphasis.cfg.enable_flag)
    {
        for(int i = 1; i < AWI_MIC_CHANNEL + 1; i++)
        {
            awi_pre_emphasis_process(&st->emphasis, st->audioFrame + i * AWI_FRAME_LENGTH, &(st->emphasis.pre_emp_first_in_chs[i-1]));
        }

        awi_pre_emphasis_process(&st->emphasis, st->audioFrame + (AWI_MIC_CHANNEL + 1) * AWI_FRAME_LENGTH, &(st->emphasis.pre_emp_first_in_spk));
    }

    awi_process_dft_analysis_filter_bank(&st->MicSignalFilterBank, st->anaMicBuffer, 
    st->audioFrame + 1 * AWI_FRAME_LENGTH, st->MicSignalFilterBankFrame);

    awi_process_dft_analysis_filter_bank(&st->SpkSignalFilterBank, st->anaSpkBuffer, 
    st->audioFrame + (AWI_MIC_CHANNEL + 1) * AWI_FRAME_LENGTH, st->SpkSignalFilterBankFrame);

    int offset = AWI_FRAME_BAND_COUNT * 2;
#ifdef WPE_USE
    float *wpe_sp;
#endif
    float *mic_sp, *fbf_sp, *aec_sp, *gsc_sp, *aes_sp, *bf_sp, *cng_sp, *ns_sp, *agc_sp;
    float spk_sp[2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2];
#ifdef NNVAD_USE
    float cng_sp_again[AWI_FRAME_BAND_COUNT * 2];
    memset(cng_sp_again, 0, AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
#endif

    if (st->aec_part1.init_aec_flag)
    {
        memset(spk_sp, 0, 2 * AWI_SPK_CHANNEL * AWI_FRAME_BAND_COUNT * 2 * sizeof(float));
    }

    for(int pos = 0; pos < AWI_FRAME_PER_BLOACK; ++pos)
    {
        mic_sp = st->MicSignalFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
        memcpy(spk_sp, st->SpkSignalFilterBankFrame + pos * offset, offset * sizeof(float));

        // if (pos == 0)
        // {
        //     memcpy(spk_sp, st->SpkSignalFilterBankFrame_Prev + ( pos + 1 ) * offset, offset * sizeof(float));
        //     memcpy(spk_sp + 1 * offset, st->SpkSignalFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
        // } else
        // {
        //     memcpy(spk_sp, st->SpkSignalFilterBankFrame + ( pos - 1 ) * offset, offset * sizeof(float));
        //     memcpy(spk_sp + 1 * offset, st->SpkSignalFilterBankFrame_Prev + pos * offset, offset * sizeof(float));
        // }
#ifdef WPE_USE
        wpe_sp = st->wpeOutSignalFilterBankFrame;
#endif
        fbf_sp = st->fbfOutSignalFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
        aec_sp = st->aecOutSignalFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
        gsc_sp = st->gscOutSignalFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
        aes_sp = st->aesOutSignalFilterBankFrame + pos * AWI_MIC_CHANNEL * offset;
        bf_sp  = st->bfOutSignalFilterBankFrame  + pos * 1 * offset;
        cng_sp = st->cngOutSignalFilterBankFrame + pos * 1 * offset;

#ifndef NNVAD_USE
        ns_sp  = st->nsOutSignalFilterBankFrame  + pos * 1 * offset;
        agc_sp = st->agcOut + pos * 1 * offset;
#endif

#ifdef WPE_USE
        awi_eigen_wpe_process_freq_dep(st->wpe, mic_sp, wpe_sp);
#endif

        for (int beam_id = 0; beam_id < AWI_BEAM_CHANNEL; beam_id++)
        {
#ifdef WPE_USE
            awi_fbf_process_single(&st->fbf, wpe_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
#else
            awi_fbf_process_single(&st->fbf, mic_sp, fbf_sp + beam_id * 2 * AWI_FRAME_BAND_COUNT, beam_id);
#endif
        }

        awi_subband_aec_process(&st->aec_part1, &st->aec_part2, fbf_sp, spk_sp, aec_sp);

        memcpy(gsc_sp, aec_sp, AWI_BEAM_CHANNEL * offset * sizeof(float));

        awi_nse_process_multichannel(&st->nse_gsc, gsc_sp, 0, 0);

        awi_simple_vad_process(&st->vad_gsc, st->nse_gsc.bkg_est_psd, st->nse_gsc.recur_block_psd, st->nse_gsc.avg_inst_noisy_psd);

        awi_aes_process(&st->aes, &st->aec_part1, fbf_sp, aec_sp, gsc_sp, aes_sp, st->vad_gsc.bkg_est_psd_floor);

#ifdef WPE_USE
        awi_bf_process(&st->bf, wpe_sp, aes_sp, bf_sp);
#else
        awi_bf_process(&st->bf, mic_sp, aes_sp, bf_sp);
#endif

        awi_cng_process(&st->cng, st->aes.min_aes_gain, bf_sp, cng_sp, 
                        st->vad_gsc.bkg_est_psd_floor, st->vad_gsc.recur_fullband_bkg_energy);

#ifndef NNVAD_USE
        awi_nse_process_mono(&st->nse_cng, cng_sp, st->vad_gsc.is_speech_triggered, 1);

        awi_ns_process(&st->ns, cng_sp, ns_sp, st->nse_cng.avg_inst_noisy_psd, st->nse_cng.bkg_est_psd, st->cng.pure_noise_flag);

        awi_agc_process(&st->agc, &st->aes, ns_sp, agc_sp, st->vad_gsc.fullband_snr, st->vad_gsc.critical_band_snr,
                        st->vad_gsc.is_speech_triggered, st->cng.pure_noise_flag, st->ns.ns_gain_seq);
#else
        st->aec_on_flag[AWI_FRAME_PER_BLOACK*(NNVAD_DELAY_FRAMES-1)+pos] = st->aec_part1.is_aec_on;

        st->bkg_energy[AWI_FRAME_PER_BLOACK*(NNVAD_DELAY_FRAMES-1)+pos] = st->vad_gsc.recur_fullband_bkg_energy;

        memcpy(st->bkg_est_psd + AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * (NNVAD_DELAY_FRAMES-1) + pos * AWI_FRAME_BAND_COUNT, 
                st->vad_gsc.bkg_est_psd_floor, AWI_FRAME_BAND_COUNT * sizeof(float));
#endif
    }


#ifdef NNVAD_USE
    awi_process_dft_synthesis_filter_bank(&st->CngOutSynthesisFilterBank, st->CngOutSynBuffer, 
                st->cngOutSignalFilterBankFrame, st->CngOutAudioFrame);

    memcpy(st->cng_fd_delay_buffer + (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * offset, 
                st->cngOutSignalFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));

    awi_limiter_process(&st->auxi_limiter, st->CngOutAudioFrame, st->cng.pure_noise_flag && st->vad_gsc.is_speech_triggered);

    float vad_prob = awi_nnvad_determine_f(st->nnvad_st, st->CngOutAudioFrame);

    #if DEBUG_VAD
    if ( vad_prob < 0 )
    {
        printf("vad_flag:%.4f\n", vad_prob);
    }
    
    if (nnvad_file == NULL)
    {
        nnvad_file = fopen("nnvad_algo_out.wav", "wb");
        short tmp[256];
        for (size_t i = 0; i < 256; i++)
        {
            tmp[i] = vad_prob * 16384.0f;
        }
        int len = fwrite(tmp, sizeof(short), 256, nnvad_file);
    }
    else
    {
        short tmp[256];
        for (size_t i = 0; i < 256; i++)
        {
            tmp[i] = vad_prob * 16384.0f;
        }
        int len = fwrite(tmp, sizeof(short), 256, nnvad_file);
    }
    #endif
    int inst_nn_vad, com_vad;
    for (int pos = 0; pos < AWI_FRAME_PER_BLOACK; pos++)
    {
        cng_sp = st->cng_fd_delay_buffer + pos * 1 * offset;
        ns_sp  = st->nsOutSignalFilterBankFrame  + pos * 1 * offset;
        agc_sp = st->agcOut + pos * 1 * offset;
        inst_nn_vad = vad_prob > 0.50f ? 1 : 0;

        if ( inst_nn_vad )
        {
            st->vad_gsc.nn_speech_counter++;
            if ( st->vad_gsc.nn_speech_counter >= st->vad_gsc.cfg.nn_speech_counter_th )
            {
                st->vad_gsc.nn_post_vad_flag = 1;
                st->vad_gsc.nn_non_speech_counter = 0;
            }
        }
        else
        {
            st->vad_gsc.nn_non_speech_counter++;
            if ( st->vad_gsc.nn_non_speech_counter >= st->vad_gsc.cfg.nn_non_speech_counter_th )
            {
                st->vad_gsc.nn_post_vad_flag = 0;
                st->vad_gsc.nn_speech_counter = 0;
            }
        }
        com_vad = st->vad_gsc.nn_post_vad_flag;

        awi_cng_process_mod(&st->cng, cng_sp, cng_sp_again,
                    com_vad, st->bkg_est_psd + pos * AWI_FRAME_BAND_COUNT, st->bkg_energy[pos]);

        awi_nse_process_mono(&st->nse_cng, cng_sp_again, com_vad, 1);

        awi_ns_process(&st->ns, cng_sp_again, ns_sp, st->nse_cng.avg_inst_noisy_psd, st->nse_cng.bkg_est_psd, st->cng.pure_noise_flag && com_vad);

        awi_agc_process(&st->agc, &st->aes, ns_sp, agc_sp, st->vad_gsc.fullband_snr, st->vad_gsc.critical_band_snr,
                        com_vad, st->cng.pure_noise_flag && com_vad, st->ns.inst_post_snr_cng);

    }

    memmove(st->cng_fd_delay_buffer, st->cng_fd_delay_buffer + AWI_FRAME_PER_BLOACK * offset, 
            (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * offset * sizeof(float));

    memmove(st->bkg_est_psd, st->bkg_est_psd + AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT, 
            (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * sizeof(float));

    memmove(st->aec_on_flag, st->aec_on_flag + AWI_FRAME_PER_BLOACK, 
            (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * sizeof(float));

    memmove(st->bkg_energy, st->bkg_energy + AWI_FRAME_PER_BLOACK, 
            (NNVAD_DELAY_FRAMES-1) * AWI_FRAME_PER_BLOACK * sizeof(float));
#endif

    memcpy(st->MicSignalFilterBankFrame_Prev, st->MicSignalFilterBankFrame, AWI_MIC_CHANNEL * AWI_FRAME_PER_BLOACK * offset * sizeof(float));
    
    memcpy(st->SpkSignalFilterBankFrame_Prev, st->SpkSignalFilterBankFrame, AWI_FRAME_PER_BLOACK * offset * sizeof(float));

    awi_process_dft_synthesis_filter_bank(&st->voipDftSynthesisFilterBank, st->voipSynBuffer, st->agcOut, st->voipOutSignalAudioFrame);

    if (st->emphasis.cfg.enable_flag)
    {
        awi_de_emphasis_process(&st->emphasis, st->voipOutSignalAudioFrame);
    }

    awi_limiter_process(&st->voip_limiter, st->voipOutSignalAudioFrame, st->cng.pure_noise_flag && st->vad_gsc.is_speech_triggered);

    awi_noise_gate_process(&st->noise_gate, st->voipOutSignalAudioFrame);

    for(int i = 0; i < AWI_FRAME_LENGTH; ++i) 
    {
        multi_chs_in[i * ( AWI_MIC_CHANNEL + AWI_SPK_CHANNEL + 1 ) + 0] = st->voipOutSignalAudioFrame[i] * 32768.f;
        output[i] = st->voipOutSignalAudioFrame[i] * 32768.f;
    }

}


int awi_algo_deinit(void *algo)
{
    awi_algo_modules_accu_t *st;
    st = (awi_algo_modules_accu_t *)algo;
#ifdef NNVAD_USE
    awi_nnvad_agc_free(st->nnvad_st);
#endif
    if ( NULL != st )
    {
        free(st);
        st = NULL;
        return 0;
    }
    else
    {
        return -1;
    }
}


#ifdef VERSION_CONTROL
#include "algo_config.h"
char* awi_algo_get_version()
{
	return AWI_ALGO_VERSION;
}

char* awi_algo_get_compile_time()
{
	return AWI_ALGO_COMPILE_TIME;
}
#endif
