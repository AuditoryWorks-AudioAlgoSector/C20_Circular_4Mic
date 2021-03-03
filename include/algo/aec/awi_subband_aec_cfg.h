#ifndef _AWI_SUBBAND_AEC_CFG_H_
#define _AWI_SUBBAND_AEC_CFG_H_
#include <stdbool.h>


typedef struct awi_subband_aec_cfg
{
    int taps;                          /* number of aec taps               */
    int bandCnts;                      /* fft length or frame length: 256  */
    int micChannels;                   /* number of aec channels: 8        */
    int spkChannels;                   /* number of spk channels: 1        */
    int frameBandCnts;                 /* frequency bins: 129              */
    int vssLowBand;
    int vssHighBand;
    float gamma;
    float coh_reg_factor;
    float stat_ff;
    float divergT;                     /* threshold for adaptive filter divergence judgement */
    float alpha_far_at;
    float alpha_far_ar;
    float spk_energy_th;               /* threshold to enable aec                            */
    float sig_th;
    float coh_th;

    float factor_fast;
    float mufb_fast;
    float err_step_size_scale;
    float dtd_reg_factor_basis;
    float reg_factor_basis;
    float dtd_reg_factor_basis_head;
    float reg_factor_basis_head;
    float copy_fast_energy_th;
    float copy_slow_energy_th;

    float alpha_eo_energy;
    float spk_energy_head_th;
    int header_counter_th;

} awi_subband_aec_cfg_t;


void awi_subband_aec_cfg_init(awi_subband_aec_cfg_t *cfg);



#endif
