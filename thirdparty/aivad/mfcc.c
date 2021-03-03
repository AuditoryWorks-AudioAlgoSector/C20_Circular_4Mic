#include "mfcc.h"


void mfcc(float *mfcc_result,float *powspectrum,float *data, int Sr, int n_filt, int n_cep, int n_fft,
          int cep_lifter, int append_Energy) {

    float PI = 3.14159265358979323846264338327f;
    float FLT_EPSILON =1.192092896e-07f;

    //初始化
    fftwf_plan p;
    fftwf_complex temp_out[n_fft / 2 + 1];
    p = fftwf_plan_dft_r2c_1d(n_fft, data, temp_out, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);

    // Compute power spectra
    for (int i = 0; i < n_fft / 2 + 1; i++){
        powspectrum[i] = (1.0f / n_fft) * (powf(temp_out[i][0], 2) + powf(temp_out[i][1], 2));
    }
    // Compute the filterbank parameters
    float **filterbank = (float **)malloc(n_filt * sizeof(float *));
	for (int i = 0; i < n_filt; i++)
		filterbank[i] = (float *)malloc((n_fft / 2 + 1) * sizeof(float ));

    get_filterbank_parameters(filterbank, n_filt, Sr, n_fft);

    // 计算能量
    float specenergy = 0.0f;
    for (int i = 0; i < n_fft / 2 + 1; i++)
        specenergy += powspectrum[i];
    if (specenergy <= 0.0f)
        specenergy = FLT_EPSILON;

    // Get filter bank output
    float feat[n_filt];
    for (int i = 0; i < n_filt; i++) {
        feat[i] = 0.0f;
        for (int j = 0; j < n_fft / 2 + 1; j++)
            feat[i] += powspectrum[j] * filterbank[i][j];

        if (feat[i] > 0.0)
            feat[i] = logf(feat[i]);
        else
            feat[i] = FLT_EPSILON;
    }
    for (int k = 0; k <  n_filt; ++k) {
        free(filterbank[k]);
    }
    free(filterbank);

    // DCT - II of filter bank output
    for (int i = 0; i < n_cep; i++)
    {
        mfcc_result[i] = 0.0f;
        for (int j = 0; j < n_filt; j++)
            mfcc_result[i] += feat[j] * cosf((i * PI / n_filt) * (j + 0.5));
        // Orthogonalization of DCT output
        if (i == 0)
            mfcc_result[i] *= sqrtf(1.0 / n_filt);
        else
            mfcc_result[i] *= sqrtf(2.0 / n_filt);

        // Ceplifter
        if (cep_lifter != 0)
            mfcc_result[i] *= 1.0f + (cep_lifter / 2.0) * sinf(PI * i / cep_lifter);
    }
    // Append Energy
    if (append_Energy == 1)
        mfcc_result[0] = logf(specenergy);
}

void get_filterbank_parameters(float **fbank, int n_filt, int sampleRate, int n_fft) {

    float lowmel = hztomel(0.0f);
    float highmel = hztomel(sampleRate / 2.0f);

    // Generate n_filt center frequencies linearly spaced in the mel scale
//    float bin[n_filt + 2];
    float bin[n_filt + 2];
    for (int i = 0; i <= n_filt + 1; i++){
        bin[i] = floorf(meltohz(i * (highmel - lowmel) / (n_filt + 1) + lowmel) * (n_fft + 1) / sampleRate);
    }
    // Triangular Filter Banks
    for (int i = 0; i < n_filt; i++)
	{
		memset(fbank[i], 0, (n_fft / 2 + 1)*sizeof(float ));
		for (int j = (int)bin[i]; j < (int)bin[i + 1]; j++)
			fbank[i][j] = (j - bin[i]) / (bin[i + 1] - bin[i]);
		for (int j = (int)bin[i + 1]; j < (int)bin[i + 2]; j++)
			fbank[i][j] = (bin[i + 2] - j) / (bin[i + 2] - bin[i + 1]);
	}
}

float hztomel(float hz) {
    return 2595 * log10f(1 + hz / 700.0f);
}

float meltohz(float mel) {
    return 700 * (powf(10, mel / 2595.0f) - 1);
}
