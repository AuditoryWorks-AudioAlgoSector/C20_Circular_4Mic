#ifndef _AWI_SIMPLE_VAD_CFG_H_
#define _AWI_SIMPLE_VAD_CFG_H_


typedef struct awi_simple_vad_cfg
{
    int snr_counter_th;
    int speech_counter_th;
    int non_speech_counter_th;
    int nn_speech_counter_th;
    int nn_non_speech_counter_th;

    float subband_snr_ratio;
    float fullband_snr_ratio;
    float alpha_valid_bkg_fast;
    float alpha_valid_bkg_slow;
    float alpha_valid_bkg;
    float alpha_local_snr;
    float alpha_global_snr;
    float bkg_energy_high_th;

} awi_simple_vad_cfg_t;


void awi_simple_vad_cfg_init(awi_simple_vad_cfg_t *cfg);
#endif
