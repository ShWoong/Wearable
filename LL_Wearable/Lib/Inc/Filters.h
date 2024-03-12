
#ifndef FILTERS_H
#define FILTERS_H

#include <stdint.h>

#define FILTER_ORDER 2
#define SECTIONS 2
/*#define FILTER_TAP_NUM 100
#define SAMPLE_SIZE 50
#define FILTER_TAP_NUM 100
#define HPFILTER_TAP_NUM 101*/

float BWHPF(float input);
float BWLPF(float input);
float NEURAL_ACTIVATION(float emg);
float MUSCLE_ACTIVATION(float neural_activation);
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
