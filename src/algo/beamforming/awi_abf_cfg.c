#include <stdlib.h>
#include "awi_abf_cfg.h"
#include "awi_utils.h"

void awi_abf_cfg_init(awi_abf_cfg_t *p)
{
    p->aicTap             = AWI_ABF_TAP;
    p->references         = AWI_ABF_REFERENCE;
    p->AFlen              = p->aicTap * p->references;
    p->vss_div_factor     = 1e-20f / ( 25 * 25 );
    p->vss_stat_factor    = 0.75;
    p->alp                = 0.45;
    p->mufb               = 1;
    p->midPoint           = 1.2;
    p->slope              = 10;

    p->minSNR             = 0.03;
    p->maxSNR             = 33;
    p->leaking_factor     = 0.99;
    p->power_reg_factor   = 5e-7f;
    p->divergence_th      = 3.5;
    p->alpha_snr          = 0.2;

}
