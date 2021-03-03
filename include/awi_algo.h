#ifndef _SRC_ALGO_AWI_ALGO_H_
#define _SRC_ALGO_AWI_ALGO_H_

#ifdef __cplusplus
extern "C"{
#endif


// #define  NNVAD_USE
// #define  WPE_USE

void *awi_algo_init();

void awi_algo_process(void *algo, short *multi_chs_in, short *output);

int awi_algo_deinit(void *algo);

#ifdef VERSION_CONTROL
char* awi_algo_get_version(void);
char* awi_algo_get_compile_time(void);
#endif

#ifdef __cplusplus
};
#endif
#endif

