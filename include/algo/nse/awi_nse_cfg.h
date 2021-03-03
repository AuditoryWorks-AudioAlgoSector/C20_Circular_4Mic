#ifndef AWI_NSE_CFG_H
#define AWI_NSE_CFG_H

enum {
    kGlobalWindow60  = 0,
    kGlobalWindow120 = 1,
};



typedef struct awi_nse_cfg
{
    int nis_blocks;
    float epsth;
    float alpha_init_bkg;
    float alpha_init_speech;
    float alpha_correction_factor;
    float alpha_opt_max;
    float alpha_opt_min_max;
    float beta_max;
    float snr_power;
    int D;
    float MD;
    float HD;
    int V;
    float MV;
    float HV;
    int U;
    int sub_win_counter;
    int local_ms_counter;
    float alpha_v;
    int global_win_len_flag;           /* value must be within {0, 1} */
    float alpha_bkg_inc;
    float alpha_bkg_dec;
    float bkg_update_th;
    int high_freq_id_th;
} awi_nse_cfg_t;


void awi_nse_cfg_init(awi_nse_cfg_t *p, int global_win_flag);


#endif
