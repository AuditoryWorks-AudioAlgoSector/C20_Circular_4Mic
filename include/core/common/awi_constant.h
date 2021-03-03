#ifndef _AWI_CONSTANT_H_
#define _AWI_CONSTANT_H_

#define ANA_FILTER_LENGTH_32    128
#define ANA_FILTER_LENGTH_64    256
#define ANA_FILTER_LENGTH_128   256
#define ANA_FILTER_LENGTH_256   512
#define ANA_FILTER_LENGTH_512   512

#define FRAME_LENGTH_32 		256
#define BAND_COUNTS_32 			32
#define DECI_RATE_32 			16
#define FRAME_BAND_COUNTS_32 	17

#define FRAME_LENGTH_64  		128
#define BAND_COUNTS_64  		64
#define DECI_RATE_64  			32
#define FRAME_BAND_COUNTS_64  	33

#define FRAME_LENGTH_128  		128
#define BAND_COUNTS_128  		128
#define DECI_RATE_128  			64
#define FRAME_BAND_COUNTS_128  	65

#define FRAME_LENGTH_256 		256
#define BAND_COUNTS_256 		256
#define DECI_RATE_256 			128
#define FRAME_BAND_COUNTS_256 	129

#define FRAME_LENGTH_512 		256
#define BAND_COUNTS_512 		512
#define DECI_RATE_512 			256
#define FRAME_BAND_COUNTS_512 	257

#define AWI_FRAME_LENGTH  		FRAME_LENGTH_512

#define AWI_MIC_CHANNEL   		4
#define AWI_SPK_CHANNEL   		1
#define AWI_BEAM_CHANNEL        AWI_MIC_CHANNEL
#define AWI_AEC_CHANNELS        AWI_BEAM_CHANNEL
#define AWI_FBF_MIC_NUM         4


#define AWI_BAND_COUNT   		256
#define AWI_AEC_TAP		  		25
#define AWI_AEC_TAP_LOWEST      8

#define AWI_ABF_TAP		1
#define AWI_ABF_REFERENCE  3

#define AWI_FREQ_BANDS 			24
#define AWI_FREQ_BANDS_VALID    22

#define AWI_BKG_V				10
#define AWI_BKG_U_SMALL         6
#define AWI_BKG_U_LARGE         12
#define AWI_BKG_U				AWI_BKG_U_LARGE

#define AWI_DOA_BIN_NUM         AWI_MIC_CHANNEL
#define AWI_DOA_EST_WIN_LEN     13
#define AWI_BKG_EST_WIN_LEN     1250


#define CNG_RAN_NUM				5

#define AWI_EPS					2.2204e-16f

#define TOTAL_TAP_NUM           2150

#define NNVAD_DELAY_FRAMES      4

#define NS_LOCAL_WIN_LEN        3

#define NS_GLOBAL_WIN_LEN       31

#endif
