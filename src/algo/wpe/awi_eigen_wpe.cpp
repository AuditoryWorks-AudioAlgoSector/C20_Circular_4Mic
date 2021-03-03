#include "awi_eigen_wpe.h"
#include "awi_constant.h"
#include "awi_filterbankparams.h"
#include "awi_filterbank.h"
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <iostream>

using namespace Eigen;

struct eigen_wpe {

    awi_wpe_cfg_t *cfg;
    bool real_denominator_flag;
    bool alternate_flag;

    MatrixXcf sp_buffer;
    MatrixXcf filter_taps;
    MatrixXcf inv_cov;

    VectorXcf window;
    RowVectorXcf row_window;
    VectorXcf kalman_gain;
    VectorXcf nominator;
    VectorXcf denominator;

    VectorXf power;

    MatrixXcf bandMat_Input_sp;
    MatrixXcf bandMat_Out_sp;
    VectorXcf bandMat_Out_sp_col;
};

awi_eigen_wpe_t* awi_eigen_wpe_create(){
    awi_eigen_wpe_t* st = (awi_eigen_wpe_t*)malloc(sizeof(awi_eigen_wpe_t));
    return st;
}

void awi_eigen_wpe_init(awi_eigen_wpe_t *wpe, awi_wpe_cfg_t *wpe_cfg)
{
    wpe->cfg = wpe_cfg;
    wpe->real_denominator_flag = 0;
    wpe->alternate_flag        = 1;
    unsigned int tapMultiCh = wpe_cfg->taps * wpe_cfg->channels; 
    unsigned int tapPlusDelayPlus1 = wpe_cfg->taps + wpe_cfg->predict_delay + 1;

    wpe->sp_buffer   = MatrixXcf::Zero(tapPlusDelayPlus1 * wpe_cfg->channels, AWI_FRAME_BAND_COUNT);
    wpe->filter_taps = MatrixXcf::Zero(tapMultiCh * wpe_cfg->channels, AWI_FRAME_BAND_COUNT);
    wpe->inv_cov     = MatrixXcf::Zero(tapMultiCh *tapMultiCh, AWI_FRAME_BAND_COUNT);

    for (unsigned int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
    {
        Map<MatrixXcf>bandMat_Inv_cov(wpe->inv_cov.col(band).data(), tapMultiCh, tapMultiCh);
        bandMat_Inv_cov.setIdentity();
    }

    wpe->window        = VectorXcf::Zero(tapMultiCh);
    wpe->row_window    = RowVectorXcf::Zero(tapMultiCh);
    wpe->kalman_gain   = VectorXcf::Zero( tapMultiCh);
    wpe->nominator     = VectorXcf::Zero(tapMultiCh);
    wpe->denominator   = VectorXf::Zero(1);
    wpe->power         = VectorXf::Zero(AWI_FRAME_BAND_COUNT);

    wpe->bandMat_Input_sp    = MatrixXcf::Zero(AWI_FRAME_BAND_COUNT, wpe->cfg->channels);
    wpe->bandMat_Out_sp      = MatrixXcf::Zero( wpe->cfg->channels, AWI_FRAME_BAND_COUNT);
    wpe->bandMat_Out_sp_col  = VectorXf::Zero(wpe->cfg->channels);

    for (int i = 0; i < AWI_FRAME_BAND_COUNT; i++){
        WPE_FREQ_TAPS[i] = (int)(WPE_FREQ_TAPS[i]/10.0f * wpe_cfg->taps + 0.5f);
        printf("%d",WPE_FREQ_TAPS[i]);
        if (i%10  == 0)
        {
            printf("\n");
        }
    }
    printf("\n");
}



void awi_eigen_wpe_process(awi_eigen_wpe_t *wpe, float *input_sp, float *out_sp)
{

    unsigned int tapMultiCh = wpe->cfg->taps * wpe->cfg->channels; 
    unsigned int buf_vec_len = ( wpe->cfg->taps + wpe->cfg->predict_delay + 1 ) * wpe->cfg->channels; 
    unsigned int shift_buf_len = buf_vec_len - wpe->cfg->channels; 

    float alpha_psd_inv = 1.f / wpe->cfg->alpha_psd;

    memcpy(wpe->bandMat_Input_sp.data(), input_sp, AWI_FRAME_BAND_COUNT * wpe->cfg->channels * 2 * sizeof(float));

    if (wpe->real_denominator_flag)
    {
        float real_denominator;
        for (unsigned int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            // predict
            Map<VectorXcf>sp_buffer_col(wpe->sp_buffer.col(band).data(), buf_vec_len);


            wpe->window = sp_buffer_col.head(tapMultiCh);
            wpe->row_window = wpe->window.adjoint();
            Map<MatrixXcf>bandMat_Filter_tap(wpe->filter_taps.col(band).data(), wpe->cfg->channels, tapMultiCh);
            wpe->bandMat_Out_sp_col = sp_buffer_col.tail(wpe->cfg->channels) - bandMat_Filter_tap.conjugate() * wpe->window;
            wpe->bandMat_Out_sp.col(band) = wpe->bandMat_Out_sp_col;
            

            // update buffer
            sp_buffer_col.head(shift_buf_len) = sp_buffer_col.tail(shift_buf_len);
            sp_buffer_col.tail(wpe->cfg->channels) = wpe->bandMat_Input_sp.row(band);

            if (wpe->alternate_flag)
            {
                // update power
                wpe->power(band) = sp_buffer_col.cwiseAbs2().mean();

                // compute kalman gain
                Map<MatrixXcf>bandMat_Inv_cov(wpe->inv_cov.col(band).data(), tapMultiCh, tapMultiCh);
                wpe->nominator = bandMat_Inv_cov * wpe->window;
                real_denominator = wpe->cfg->alpha_psd * wpe->power(band) + (wpe->row_window * wpe->nominator)(0).real() + AWI_EPS;
                wpe->kalman_gain = wpe->nominator / real_denominator;

                // update inv_cov
                bandMat_Inv_cov -= wpe->kalman_gain * ( wpe->row_window * bandMat_Inv_cov);
                bandMat_Inv_cov *= alpha_psd_inv;

                // update filter coefficients
                bandMat_Filter_tap += wpe->bandMat_Out_sp_col.conjugate() * wpe->kalman_gain.transpose();
            }
        }
    }
    else
    {
        float sp_re, sp_im;
        float denominator_square_am;
        for (unsigned int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            // predict
            Map<VectorXcf>sp_buffer_col(wpe->sp_buffer.col(band).data(), buf_vec_len);
            wpe->window = sp_buffer_col.head(tapMultiCh);
            wpe->row_window = wpe->window.adjoint();
            Map<MatrixXcf>bandMat_Filter_tap(wpe->filter_taps.col(band).data(), wpe->cfg->channels, tapMultiCh);
            wpe->bandMat_Out_sp_col = sp_buffer_col.tail(wpe->cfg->channels) - bandMat_Filter_tap.conjugate() * wpe->window;
            wpe->bandMat_Out_sp.col(band) = wpe->bandMat_Out_sp_col;
            

            // update buffer
            sp_buffer_col.head(shift_buf_len) = sp_buffer_col.tail(shift_buf_len);
            sp_buffer_col.tail(wpe->cfg->channels) = wpe->bandMat_Input_sp.row(band);
            
            if (wpe->alternate_flag)
            {
                // update power
                wpe->power(band) = sp_buffer_col.cwiseAbs2().mean();

                // compute kalman gain
                Map<MatrixXcf>bandMat_Inv_cov(wpe->inv_cov.col(band).data(), tapMultiCh, tapMultiCh);
                wpe->nominator = bandMat_Inv_cov * wpe->window;
                wpe->denominator(0).real(wpe->cfg->alpha_psd * wpe->power(band));
                wpe->denominator(0).imag(0);
                wpe->denominator += wpe->row_window * wpe->nominator;
                sp_re = wpe->denominator(0).real();
                sp_im = wpe->denominator(0).imag();
                denominator_square_am = sp_re * sp_re + sp_im * sp_im + AWI_EPS;
                wpe->denominator(0).imag(-sp_im);
                wpe->kalman_gain = wpe->nominator * wpe->denominator / denominator_square_am;

                // update inv_cov
                bandMat_Inv_cov -= wpe->kalman_gain * ( wpe->row_window * bandMat_Inv_cov);
                bandMat_Inv_cov *= alpha_psd_inv;

                // update filter coefficients
                bandMat_Filter_tap += wpe->bandMat_Out_sp_col.conjugate() * wpe->kalman_gain.transpose();
            }  
        }
    }
    
    wpe->bandMat_Input_sp = wpe->bandMat_Out_sp.transpose();
    memcpy(out_sp, wpe->bandMat_Input_sp.data(), AWI_FRAME_BAND_COUNT * wpe->cfg->channels * 2 * sizeof(float));

    wpe->alternate_flag = !wpe->alternate_flag;

}

void awi_eigen_wpe_process_freq_dep(awi_eigen_wpe_t *wpe, float *input_sp, float *out_sp)
{

    unsigned int tapMultiCh;
    unsigned int buf_vec_len;
    unsigned int shift_buf_len; 
    unsigned int valid_buf_vec_len;
    unsigned int valid_power_vec_len;
    unsigned int valid_buf_vec_offset;
    unsigned int max_taps = wpe->cfg->taps;

    float alpha_psd_inv = 1.f / wpe->cfg->alpha_psd;

    VectorXcf valid_window;
    RowVectorXcf valid_row_window;
    MatrixXcf sub_bandMat_Filter_tap;
    MatrixXcf sub_bandMat_Inv_cov;
    VectorXcf valid_nominator;
    VectorXcf valid_kalman_gain;

    memcpy(wpe->bandMat_Input_sp.data(), input_sp, AWI_FRAME_BAND_COUNT * wpe->cfg->channels * 2 * sizeof(float));

    if (wpe->real_denominator_flag)
    {
        float real_denominator;
        for (unsigned int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            tapMultiCh = wpe->cfg->taps * wpe->cfg->channels; 
            buf_vec_len = ( wpe->cfg->taps + wpe->cfg->predict_delay + 1 ) * wpe->cfg->channels; 
            shift_buf_len = buf_vec_len - wpe->cfg->channels; 
            valid_buf_vec_offset = (max_taps - WPE_FREQ_TAPS[band]) * wpe->cfg->channels; 
            valid_buf_vec_len = tapMultiCh - valid_buf_vec_offset;
            valid_power_vec_len = ( WPE_FREQ_TAPS[band] + wpe->cfg->predict_delay + 1 ) * wpe->cfg->channels; 

            // predict
            Map<VectorXcf>sp_buffer_col(wpe->sp_buffer.col(band).data(), buf_vec_len);
            valid_window = sp_buffer_col.segment(valid_buf_vec_offset, valid_buf_vec_len);
            valid_row_window = valid_window.adjoint();
            Map<MatrixXcf>bandMat_Filter_tap(wpe->filter_taps.col(band).data(), wpe->cfg->channels, tapMultiCh);
            sub_bandMat_Filter_tap = bandMat_Filter_tap.rightCols(valid_buf_vec_len);
            wpe->bandMat_Out_sp_col = sp_buffer_col.tail(wpe->cfg->channels) - sub_bandMat_Filter_tap.conjugate() * valid_window;
            wpe->bandMat_Out_sp.col(band) = wpe->bandMat_Out_sp_col;
            

            // update buffer
            sp_buffer_col.head(shift_buf_len) = sp_buffer_col.tail(shift_buf_len);
            sp_buffer_col.tail(wpe->cfg->channels) = wpe->bandMat_Input_sp.row(band);

            if (wpe->alternate_flag)
            {
                // update power
                wpe->power(band) = sp_buffer_col.tail(valid_power_vec_len).cwiseAbs2().mean();

                // compute kalman gain
                Map<MatrixXcf>bandMat_Inv_cov(wpe->inv_cov.col(band).data(), tapMultiCh, tapMultiCh);
                sub_bandMat_Inv_cov = bandMat_Inv_cov.block(valid_buf_vec_offset, valid_buf_vec_offset, valid_buf_vec_len, valid_buf_vec_len);
                valid_nominator = sub_bandMat_Inv_cov * valid_window;
                real_denominator = wpe->cfg->alpha_psd * wpe->power(band) + (valid_row_window * valid_nominator)(0).real() + AWI_EPS;
                valid_kalman_gain = valid_nominator / real_denominator;

                // update inv_cov
                sub_bandMat_Inv_cov -= valid_kalman_gain * ( valid_row_window * sub_bandMat_Inv_cov);
                sub_bandMat_Inv_cov *= alpha_psd_inv;
                bandMat_Inv_cov.block(valid_buf_vec_offset, valid_buf_vec_offset, valid_buf_vec_len, valid_buf_vec_len) = sub_bandMat_Inv_cov;

                // update filter coefficients
                sub_bandMat_Filter_tap += wpe->bandMat_Out_sp_col.conjugate() * valid_kalman_gain.transpose();
                bandMat_Filter_tap.rightCols(valid_buf_vec_len) = sub_bandMat_Filter_tap;
            }
        }
    }
    else
    {
        float sp_re, sp_im;
        float denominator_square_am;
        for (unsigned int band = 0; band < AWI_FRAME_BAND_COUNT; band++)
        {
            tapMultiCh = wpe->cfg->taps * wpe->cfg->channels; 
            buf_vec_len = ( wpe->cfg->taps + wpe->cfg->predict_delay + 1 ) * wpe->cfg->channels; 
            shift_buf_len = buf_vec_len - wpe->cfg->channels; 
            valid_buf_vec_offset = (max_taps -  WPE_FREQ_TAPS[band]) * wpe->cfg->channels; 
            valid_buf_vec_len = tapMultiCh - valid_buf_vec_offset;
            valid_power_vec_len = ( WPE_FREQ_TAPS[band] + wpe->cfg->predict_delay + 1 ) * wpe->cfg->channels; 

            // predict
            Map<VectorXcf>sp_buffer_col(wpe->sp_buffer.col(band).data(), buf_vec_len);
            valid_window = sp_buffer_col.segment(valid_buf_vec_offset, valid_buf_vec_len);
            valid_row_window = valid_window.adjoint();
            Map<MatrixXcf>bandMat_Filter_tap(wpe->filter_taps.col(band).data(), wpe->cfg->channels, tapMultiCh);
            sub_bandMat_Filter_tap = bandMat_Filter_tap.rightCols(valid_buf_vec_len);
            // std::cout << sub_bandMat_Filter_tap - bandMat_Filter_tap << std::endl;
           
            wpe->bandMat_Out_sp_col = sp_buffer_col.tail(wpe->cfg->channels) - sub_bandMat_Filter_tap.conjugate() * valid_window;
            wpe->bandMat_Out_sp.col(band) = wpe->bandMat_Out_sp_col;
            

            // update buffer
            sp_buffer_col.head(shift_buf_len) = sp_buffer_col.tail(shift_buf_len);
            sp_buffer_col.tail(wpe->cfg->channels) = wpe->bandMat_Input_sp.row(band);
            
            if (wpe->alternate_flag)
            {
                // update power
                wpe->power(band) = sp_buffer_col.tail(valid_power_vec_len).cwiseAbs2().mean();

                // compute kalman gain
                Map<MatrixXcf>bandMat_Inv_cov(wpe->inv_cov.col(band).data(), tapMultiCh, tapMultiCh);
                sub_bandMat_Inv_cov = bandMat_Inv_cov.block(valid_buf_vec_offset, valid_buf_vec_offset, valid_buf_vec_len, valid_buf_vec_len);
                valid_nominator = sub_bandMat_Inv_cov * valid_window;
                wpe->denominator(0).real(wpe->cfg->alpha_psd * wpe->power(band));
                wpe->denominator(0).imag(0);
                wpe->denominator += valid_row_window * valid_nominator;
                sp_re = wpe->denominator(0).real();
                sp_im = wpe->denominator(0).imag();
                denominator_square_am = sp_re * sp_re + sp_im * sp_im + AWI_EPS;
                wpe->denominator(0).imag(-sp_im);
                valid_kalman_gain = valid_nominator * wpe->denominator / denominator_square_am;

                // update inv_cov
                sub_bandMat_Inv_cov -= valid_kalman_gain * ( valid_row_window * sub_bandMat_Inv_cov);
                sub_bandMat_Inv_cov *= alpha_psd_inv;
                bandMat_Inv_cov.block(valid_buf_vec_offset, valid_buf_vec_offset, valid_buf_vec_len, valid_buf_vec_len) = sub_bandMat_Inv_cov;

                // update filter coefficients
                sub_bandMat_Filter_tap += wpe->bandMat_Out_sp_col.conjugate() * valid_kalman_gain.transpose();
                bandMat_Filter_tap.rightCols(valid_buf_vec_len) = sub_bandMat_Filter_tap;
            }  
        }
    }
    
    wpe->bandMat_Input_sp = wpe->bandMat_Out_sp.transpose();
    memcpy(out_sp, wpe->bandMat_Input_sp.data(), AWI_FRAME_BAND_COUNT * wpe->cfg->channels * 2 * sizeof(float));

    wpe->alternate_flag = !wpe->alternate_flag;

}
