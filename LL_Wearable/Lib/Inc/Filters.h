
#ifndef FILTERS_H
#define FILTERS_H

#include <stdint.h>

#define SCALE_FACTOR 65536 //2^16 16bit a fixed point
#define FILTER_ORDER 2
#define SECTIONS 2
/*#define FILTER_TAP_NUM 100
#define SAMPLE_SIZE 50
#define FILTER_TAP_NUM 100
#define HPFILTER_TAP_NUM 101*/

int32_t BWHPF(int32_t input);
int32_t BWLPF(int32_t input);

/*void FIRF_Init(void);
float FIRF_Process(float input);
float MAF(float new_sample);
float EWMAF(float new_measurement, float prev_ewma, float alpha);

void HighPassFilter_Init(void);
float HighPassFilter_Process(float input);
float applyLowPassFilter(float input);

typedef struct {
    float estimate;
    float error_estimate;
    float error_measure;
    float kalman_gain;
} KMF;

void KMF_Init(KMF *kf, float init_estimate, float init_error_estimate, float error_measure);
void KMF_Update(KMF *kf, float measurement);*/

#endif /* FILTERS_H */
