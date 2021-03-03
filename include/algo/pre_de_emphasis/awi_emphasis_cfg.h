//
// Created by xzc on 19-7-26.
//

#ifndef ASP_AWI_EMPHASIS_CFG_H

#include <stdbool.h>

typedef struct awi_emphasis_cfg
{
    float alpha_pre;
    float alpha_de;
    bool enable_flag;
}awi_emphasis_cfg_t;

void awi_emphasis_cfg_init(awi_emphasis_cfg_t *p);

#define ASP_AWI_EMPHASIS_CFG_H

#endif //ASP_AWI_EMPHASIS_CFG_H
