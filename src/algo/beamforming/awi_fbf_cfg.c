#include "awi_fbf_cfg.h"
#include "awi_constant.h"
#include "awi_filterbankparams.h"

void awi_fbf_cfg_init(awi_fbf_cfg_t *cfg)
{
    cfg->mvdr_fbf_weight = FBF_PARAMS_256;

    cfg->block_angle[0]  = 2;
    cfg->block_angle[1]  = 3;
    cfg->block_angle[2]  = 4;


    cfg->block_angle[3]  = 3;
    cfg->block_angle[4]  = 4;
    cfg->block_angle[5]  = 1;


    cfg->block_angle[6]  = 4;
    cfg->block_angle[7]  = 1;
    cfg->block_angle[8]  = 2;

    cfg->block_angle[9]  = 1;
    cfg->block_angle[10] = 2;
    cfg->block_angle[11] = 3;

    for(int i = 0; i < AWI_BEAM_CHANNEL * AWI_ABF_REFERENCE; i++)
        cfg->block_angle[i] -= 1;

}