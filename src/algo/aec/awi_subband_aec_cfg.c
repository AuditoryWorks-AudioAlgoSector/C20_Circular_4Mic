#include <awi_constant.h>
#include <awi_filterbankparams.h>
#include "awi_subband_aec_cfg.h"
#include "awi_constant.h"
#include "awi_filterbankparams.h"

void awi_subband_aec_cfg_init(awi_subband_aec_cfg_t *cfg)
{
    cfg->micChannels                = AWI_BEAM_CHANNEL;
    cfg->spkChannels                = AWI_SPK_CHANNEL;
    cfg->bandCnts                   = AWI_BAND_COUNT;
    cfg->frameBandCnts              = AWI_FRAME_BAND_COUNT;
    cfg->taps                       = AWI_AEC_TAP;

    cfg->vssLowBand                 = 4;
    cfg->vssHighBand                = 68;
    cfg->alpha_far_ar               = 0.80f;
    cfg->alpha_far_at               = 0.10f;
    cfg->gamma                      = 1 - 1.f / cfg->taps;
    cfg->stat_ff                    = 1 - 1.f / 2.f / cfg->taps;
    cfg->coh_reg_factor             = 1e-20f;
    cfg->divergT                    = 3.5f;
    cfg->spk_energy_th              = 1.e-7f;
    cfg->spk_energy_head_th         = 2.e-5f;
    cfg->sig_th                     = 3.e-11f;
    cfg->coh_th                     = 3.e-10f;

    cfg->factor_fast                    = 1.5f;
    cfg->mufb_fast                      = 1.f / cfg->factor_fast;
    cfg->err_step_size_scale            = 1.5f;
    cfg->dtd_reg_factor_basis           = 1e-06f;
    cfg->reg_factor_basis               = 3e-06f;
    cfg->dtd_reg_factor_basis_head      = 8e-07f;
    cfg->reg_factor_basis_head          = 9e-07f;
    cfg->copy_fast_energy_th            = 2.0;
    cfg->copy_slow_energy_th            = 1.1;
    cfg->alpha_eo_energy                = 0.95;
    cfg->header_counter_th              = 150;

}
