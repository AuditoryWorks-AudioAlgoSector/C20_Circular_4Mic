#include "awi_nse_cfg.h"
#include "awi_constant.h"

void awi_nse_cfg_init(awi_nse_cfg_t *p, int global_win_flag)
{
    p->nis_blocks                = 30;
    p->epsth                     = 1e-35f;
    p->alpha_init_bkg            = 0.95f;
    p->alpha_init_speech         = 0.90f;
    p->alpha_correction_factor   = 0.9f;
    p->alpha_opt_max             = 0.98f;
    p->alpha_opt_min_max         = 0.7f;
    p->beta_max                  = 0.8f;
    p->snr_power                 = -0.125f;
    p->alpha_v                   = 2.12f;
    p->V                         = 10;
    p->MV                        = 0.61f;
    p->HV                        = 0.98f;

    p->global_win_len_flag       = global_win_flag;
    p->alpha_bkg_inc             = 0.98;
    p->alpha_bkg_dec             = 0.90;
    p->bkg_update_th             = 1.2589f;
    p->high_freq_id_th           = 120;

    if ( p-> global_win_len_flag == kGlobalWindow60 )
    {
        p->D                     = p->V * AWI_BKG_U_SMALL;
        p->MD                    = 0.841f;
        p->HD                    = 2.9f;
    }

    if ( p->global_win_len_flag == kGlobalWindow120 )
    {
        p->D                     = p->V * AWI_BKG_U_LARGE;
        p->MD                    = 0.89f;
        p->HD                    = 4.0f;
    }

    p->U                         = p->D / p->V;
    p->sub_win_counter           = 0;
    p->local_ms_counter          = p->V - 1;
}
