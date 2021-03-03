#include "awi_utils.h"

float awi_calculate_vss(float corr)
{
    if (corr >= 0.6f)
        return 1.;
    else if ( corr >= 0.4f && corr < 0.6f)
        return corr * corr / 0.36f;
    else
        return 0.01f;
}

float awi_sigmoid(float corr, int midPoint, int slope)
{
    float step = (corr - midPoint) * slope;
    step = 1 - .5 * (1 + step / (1 + fabs(step)));
    return fmax(fmin(step, 1), 0);
}

float awi_db2pow(float ydB)
{
    return pow(10, ydB / 10);
}

float awi_get_radius(int sampling_rate)
{
    float res;
    if (sampling_rate < 12000)
        res = 0.9;
    else if (sampling_rate < 24000)
        res = 0.92;
    else
        res = 0.992;

    return res;
}

void awi_filter_dc_notch16(const float *in, float radius, float *out, int len, float *mem)
{
    int i;
    float den2;

    den2 = radius * radius + .7f * (1 - radius) * (1 - radius);
    for (i = 0; i < len; i++)
    {
        float vin = in[i];
        float vout = mem[0] + vin;
        mem[0] = mem[1] + 2 * (-vin + radius * vout);
        mem[1] = vin - den2 * vout;
        out[i] = radius * vout;
    }
}

void awi_cross_fadef(float *output, float *newChannel, float *oldChannel, int n, int id, int totalFrame)
{
    float alph;
    for (int i = 0; i < n; ++i)
    {
        alph = cos((id * n + i) * PI / (2 * (totalFrame * n - 1)));
        output[i] = alph * oldChannel[i] + (1 - alph) * newChannel[i];
    }
}

int awi_is_active(float *arr, int length, float threshold)
{
    float sum = 0;
    for (int i = 0; i < length; ++i)
        sum += arr[i] * arr[i];

    return sum / length > threshold;
}

void awi_normalize(float *dst, short *src, int length)
{
    float factor = 1.0f / 32768;
    for (int i = 0; i < length; ++i)
    {
        int value = src[i];
        dst[i] = value * factor;
    }
}

void awi_unnormalize(short *dst, float *src, int length)
{
    for(int i = 0; i < length; ++i)
        dst[i] = src[i] * 32768;
}

int awi_get_max_index_float(float *arr, int length)
{
    float max = -FLT_MAX;
    int index = 0;
    for (int i = 0; i < length; i++)
        if (max < arr[i])
        {
            max = arr[i];
            index = i;
        }
    return index;
}

float awi_get_max_float(float *arr, int length)
{
    float max = -FLT_MAX;
    for (int i = 0; i < length; i++)
        if (max < arr[i])
            max = arr[i];
    return max;
}


float awi_get_mean_float(float *arr, int length)
{
    float mean = 0;
    for (int i = 0; i < length; ++i)
        mean += arr[i];
    return mean / length;
}

float awi_get_min_float(float *arr, int length)
{
    float min = FLT_MAX;
    for (int i = 0; i < length; i++)
        if (min > arr[i])
        {
            min = arr[i];
        }
    return min;
}

int awi_get_min_index_float(float *arr, int length)
{
    float min = FLT_MAX;
    int index;

    for (int i = 0; i < length; i++)
        if (min > arr[i])
        {
            min = arr[i];
            index = i;
        }
    return index;
}

int awi_get_max_index_int(int *arr, int length)
{
    int max = INT_MIN;
    int index = 0;
    for (int i = 0; i < length; i++)
        if (max < arr[i])
        {
            max = arr[i];
            index = i;
        }
    return index;
}

float awi_get_median_float(float *a, int n, int i)
{
    if (n == 1)
    {
        return a[0];
    }

    int n_meds = 0;
    for (int j = 0; j < n; j += 5)
    {
        int l = AWI_MIN(5, n - j);
        awi_insertion_sort_float(a + j, l);
        float tmp = a[j / 5];
        a[j / 5] = a[j + l / 2];
        a[j + l / 2] = tmp;
        n_meds++;
    }

    float median_of_medians;
    if (n_meds > 1)
    {
        median_of_medians = awi_get_median_float(a, n_meds, n_meds / 2);
    }
    else
    {
        median_of_medians = a[0];
    }

    int k = awi_median_partition_float(a, n, median_of_medians);

    if (k == i)
    {
        return median_of_medians;
    }
    else if (i < k)
    {
        return awi_get_median_float(a, k, i);
    }
    else
    {
        return awi_get_median_float(a + k, n - k, i - k);
    }
}


void awi_insertion_sort_float(float *a, int n)
{
    for (int j = 1; j < n; j++)
    {
        float key = a[j];
        int i = j - 1;
        while ((i >= 0) && (a[i] > key))
        {
            a[i + 1] = a[i];
            i--;
        }
        a[i + 1] = key;
    }
}


int awi_median_partition_float(float *a, int n, float x)
{
    for (int i = 0; i < n; i++)
    {
        if(a[i] == x)
        {
            a[i] = a[n - 1];
            a[n - 1] = x;
        }
    }

    int id = 0;
    for (int j = 0; j < (n - 1); j++)
    {
        if (a[j] <= x)
        {
            float tmp = a[j];
            a[j] = a[id];
            a[id] = tmp;
            id++;
        }
    }

    a[n - 1] = a[id];
    a[id] = x;

    return id;
}

void awi_correct_delay(float *buffer, float *data, int dataSize, int channel, int delaySize)
{
    for(int i = 0; channel; i++)
    {
        memmove(buffer, buffer + i * (dataSize + delaySize) + dataSize, sizeof(float) * delaySize);
        memcpy(buffer + i * (dataSize + delaySize) + delaySize,  data + i * dataSize, sizeof(float) * dataSize);
    }
}

float NormC(float *z)
{
    return z[0] * z[0] + z[1] * z[1];
}


void CXC(float *z, float *x, float *y)
{
    float re = x[0] * y[0] - x[1] * y[1];
    float im = x[1] * y[0] + x[0] * y[1];
    z[0] = re;
    z[1] = im;
}

void CPlusC(float *z, float *x, float *y)
{
    z[0] = x[0] + y[0];
    z[1] = x[1] + y[1];
}

void CSubtractC(float *z, float *x, float *y)
{
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
}

void awi_init_step_size(float *arr, int bins)
{
    float mid = 1;
    float slope = 5;
    float step;

    for(int i = 0; i < bins; i++)
    {
        step = 1.0 * i / (bins - 1);
        step = (step - mid) * slope;
        arr[i] = 1 - 0.5 * (1 + step / (1 + fabs(step)));
    }

    float normalizedFactor = 1.f / awi_get_max_float(arr, bins);

    for(int i = 0; i < bins; i++)
        arr[i] *= normalizedFactor;
}

void microphone_position_init(float *mic, int r, int channels)
{
    float angle = 2 * PI / channels;
    for(int i = 0; i < channels; i++)
    {
        mic[3 * i] = r * cos(angle * i);
        mic[3 * i + 1] = r * sin(angle * i);
        mic[3 * i + 2] = 0.;
    }
}

float fminf_local(float a, float b)
{
	return (a < b) ? a : b;
}

float fmaxf_local(float a, float b)
{
	return (a > b) ? a : b;
}


void awi_smoothing_1D_vector(float *dst, const float *src, const float *window, int data_length, int window_length)
{
    int half_length = ( window_length - 1 ) / 2;
    float tmp_src;
    for (int i = 0; i < data_length; i++)
    {
        dst[i] = 0;

        for (int j = i-half_length; j <= i+half_length; j++)
        {
            if ( j <= 0) 
            {
                tmp_src = src[0];
            }
            else if ( j >= data_length - 1 )
            {
                tmp_src = src[data_length - 1];
            }
            else
            {
                tmp_src = src[j];
            }

            dst[i] += tmp_src * window[j-i+half_length];
        }
    }
}