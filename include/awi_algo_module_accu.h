
#ifndef _AWI_ALGO_MODULE_ACCU_H_
#define _AWI_ALGO_MODULE_ACCU_H_


#include "awi_constant.h"
#include "awi_emphasis.h"
#include "awi_dft_analysis_filter_bank.h"
#include "awi_dft_synthesis_filter_bank.h"
#include "awi_fbf.h"
#include "awi_subband_aec.h"
#include "awi_abf.h"
#include "awi_aes.h"
#include "awi_beam_fusion.h"
#include "awi_nse.h"
#include "awi_simple_vad.h"
#include "awi_cng.h"
#include "awi_ns.h"
#include "awi_agc.h"
#include "awi_limiter.h"
#include "awi_noise_gate.h"
#include "awi_wpe_cfg.h"
#include "nn_vad_agc.h"
#include "awi_algo.h"

#include "awi_eigen_wpe.h"

typedef struct awi_algo_modules_accu
{

    float notch_mem[2 * AWI_MIC_CHANNEL];

    awi_emphasis_t emphasis;

    dft_analysis_filter_bank MicSignalFilterBank;
    dft_analysis_filter_bank SpkSignalFilterBank;

    awi_wpe_cfg_t wpe_cfg; 
    awi_eigen_wpe_t* wpe;
    awi_fbf_t fbf;

    awi_subband_aec_part1_t aec_part1;
    awi_subband_aec_part2_t aec_part2;

    awi_abf_part1_t abf_part1;
    awi_abf_part2_t abf_part2;

    awi_nse_t nse_gsc;

    awi_simple_vad_t vad_gsc;

    awi_aes_t aes;

    awi_bf_t bf;

    awi_cng_t cng;

    awi_nse_t nse_cng;

    awi_ns_t ns;

    awi_agc_t agc;

    dft_synthesis_filter_bank voipDftSynthesisFilterBank;
    dft_synthesis_filter_bank CngOutSynthesisFilterBank;
    syn_spk_buffer voipSynBuffer;
    syn_spk_buffer CngOutSynBuffer;

    awi_limiter_t voip_limiter;
    awi_limiter_t auxi_limiter;

    awi_noise_gate_t noise_gate;

    ana_mic_buffer anaMicBuffer;
    ana_spk_buffer anaSpkBuffer;

    float audioFrame[(AWI_MIC_CHANNEL + AWI_SPK_CHANNEL + 1) * AWI_FRAME_LENGTH];

    float MicSignalFilterBankFrame_Prev[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_MIC_CHANNEL * 2];
    float SpkSignalFilterBankFrame_Prev[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_SPK_CHANNEL * 2];

    float MicSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_MIC_CHANNEL * 2];
    float SpkSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_SPK_CHANNEL * 2];

	float wpeOutSignalFilterBankFrame[AWI_FRAME_BAND_COUNT * AWI_MIC_CHANNEL * 2];

    float fbfOutSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * 2];

    float aecOutSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * 2];

    float gscOutSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * 2];

    float aesOutSignalFilterBankFrame[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * AWI_BEAM_CHANNEL * 2];

    float bfOutSignalFilterBankFrame[AWI_FRAME_BAND_COUNT * AWI_FRAME_PER_BLOACK * 2];

    float cngOutSignalFilterBankFrame[AWI_FRAME_BAND_COUNT * AWI_FRAME_PER_BLOACK * 2];

    float nsOutSignalFilterBankFrame[AWI_FRAME_BAND_COUNT * AWI_FRAME_PER_BLOACK * 2];

    float agcOut[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * 2];

    float voipOutSignalAudioFrame[AWI_FRAME_LENGTH];

    float CngOutAudioFrame[AWI_FRAME_LENGTH];

#ifdef NNVAD_USE
    NnvadAgcState* nnvad_st;

    float cng_fd_delay_buffer[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * 2 * NNVAD_DELAY_FRAMES];

    float bkg_est_psd[AWI_FRAME_PER_BLOACK * AWI_FRAME_BAND_COUNT * NNVAD_DELAY_FRAMES];

    int aec_on_flag[AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES];

    float bkg_energy[AWI_FRAME_PER_BLOACK * NNVAD_DELAY_FRAMES];
#endif

}awi_algo_modules_accu_t;


#endif





