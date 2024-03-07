
#include "Filters.h"

float hpf_x_buffer[SECTIONS][2] = {0};
float hpf_y_buffer[SECTIONS][2] = {0};
float lpf_x_buffer[SECTIONS][2] = {0};
float lpf_y_buffer[SECTIONS][2] = {0};


float hpf_sos[SECTIONS][6] = {
    /*{0.78136727, -1.56273453,  0.78136727, 1.0, -1.67466095, 0.70485868},
    {1.0, -2.0,  1.0, 1.0, -1.83312526, 0.86618045}*/ //30Hz
	/*{0.8484753, -1.69695059, 0.8484753, 1.0, -1.77831349, 0.79244747},
	{1.0, -2.0, 1.0, 1.0, -1.8934156, 0.90846441}*/ //20Hz
	{0.88409204, -1.76818409, 0.88409204, 1.0, -1.83185386, 0.84001994},
	{1.0, -2.0, 1.0, 1.0, -1.92190889, 0.93047642} //15Hz
	/*{0.99182421, -1.98364842, 0.99182421, 1.0, -1.98841802, 0.98845727},
	{1.0, -2.0, 1.0, 1.0, -1.99516324, 0.99520262}*/ //1Hz
};

float lpf_sos[SECTIONS][6] = {
	/*{0.06290094, 0.12580189, 0.06290094, 1.0, -0.1964664, 0.0484845},
	{1.0, 2.0, 1.0, 1.0, -0.27237536, 0.45358867}*/ //220Hz
    {1.20231162e-07, 2.40462324e-07, 1.20231162e-07, 1.0, -1.93132782, 0.932701053},
    {1.0, 2.0, 1.0, 1.0, -1.97016249, 0.971563337} //6Hz
    /*{1.32937289e-05, 2.65874578e-05, 1.32937289e-05, 1.00000000, -1.77831349, 0.79244747},
	{1.00000000, 2.00000000, 1.00000000, 1.00000000, -1.89341560, 0.90846441}*/ //20Hz
    /*{7.69909891e-09, 1.53981978e-08, 7.69909891e-09, 1.00000000, -1.96541950, 0.96576872},
	{1.00000000, 2.00000000, 1.00000000, 1.00000000, -1.98532459, 0.98567734}*/ //3Hz
	/*{1.53324552e-09, 3.06649104e-09, 1.53324552e-09, 1.00000000, -1.97689135, 0.97704745},
	{1.00000000, 2.00000000, 1.00000000, 1.00000000, -1.99027124, 0.99042840}*/ //2Hz
    /*{9.66139663e-11, 1.93227933e-10, 9.66139663e-11, 1.00000000, -1.98841802, 0.98845727},
	{1.00000000, 2.00000000, 1.00000000, 1.00000000, -1.99516324, 0.99520262}*/ //1Hz
    /*{9.73291699e-15, 1.94658340e-14, 9.73291699e-15, 1.0, -1.99883930, 0.99883969},
	{1.0, 2.0, 1.0, 1.0, -1.99951883, 0.99951922}*/ //0.1Hz
};

float BWHPF(float input) {
    float output = 0.0;
    for (int i = 0; i < SECTIONS; i++) {
        float xn = (i == 0) ? input : output;

        output = hpf_sos[i][0] * xn + hpf_sos[i][1] * hpf_x_buffer[i][0] + hpf_sos[i][2] * hpf_x_buffer[i][1]
                 - hpf_sos[i][4] * hpf_y_buffer[i][0] - hpf_sos[i][5] * hpf_y_buffer[i][1];

        hpf_x_buffer[i][1] = hpf_x_buffer[i][0];
        hpf_x_buffer[i][0] = xn;
        hpf_y_buffer[i][1] = hpf_y_buffer[i][0];
        hpf_y_buffer[i][0] = output;
    }
    return output;
}

float BWLPF(float input) {
    float output = 0.0;
    for (int i = 0; i < SECTIONS; i++) {
        float xn = (i == 0) ? input : output;

        output = lpf_sos[i][0] * xn + lpf_sos[i][1] * lpf_x_buffer[i][0] + lpf_sos[i][2] * lpf_x_buffer[i][1]
                 - lpf_sos[i][4] * lpf_y_buffer[i][0] - lpf_sos[i][5] * lpf_y_buffer[i][1];

        lpf_x_buffer[i][1] = lpf_x_buffer[i][0];
        lpf_x_buffer[i][0] = xn;
        lpf_y_buffer[i][1] = lpf_y_buffer[i][0];
        lpf_y_buffer[i][0] = output;
    }
    return output;
}


#define LPF_SECTION_COUNT 1  // 2차 필터는 단일 섹션으로 구성됩니다.

// 필터 계수 정의
const float lpf_coefficients[LPF_SECTION_COUNT][6] = {
    {8.76555488e-05, 1.75311098e-04, 8.76555488e-05, 1.0, -1.97334425, 0.973694872}
};

// 입력 및 출력 버퍼 정의
float lpf_input_buffer[LPF_SECTION_COUNT][2] = {0};
float lpf_output_buffer[LPF_SECTION_COUNT][2] = {0};

float applyLowPassFilter(float input) {
    float output = input;

    for (int i = 0; i < LPF_SECTION_COUNT; i++) {
        float b0 = lpf_coefficients[i][0];
        float b1 = lpf_coefficients[i][1];
        float b2 = lpf_coefficients[i][2];
        float a1 = lpf_coefficients[i][4];
        float a2 = lpf_coefficients[i][5];

        // 현재 섹션의 출력 계산
        output = b0 * input + b1 * lpf_input_buffer[i][0] + b2 * lpf_input_buffer[i][1]
                 - a1 * lpf_output_buffer[i][0] - a2 * lpf_output_buffer[i][1];

        // 버퍼 업데이트
        lpf_input_buffer[i][1] = lpf_input_buffer[i][0];
        lpf_input_buffer[i][0] = input;
        lpf_output_buffer[i][1] = lpf_output_buffer[i][0];
        lpf_output_buffer[i][0] = output;

        // 다음 섹션의 입력은 현재 섹션의 출력
        input = output;
    }

    return output;
}

static float firFilterCoeffs[FILTER_TAP_NUM] = {
    // 여기에 계산된 계수 넣기
    0.00149401, 0.00151131, 0.00156314, 0.00164928, 0.00176939,
    0.00192298, 0.00210944, 0.00232802, 0.00257783, 0.00285788,
    0.00316703, 0.00350404, 0.00386755, 0.0042561 , 0.00466813,
    0.00510197, 0.00555589, 0.00602804, 0.00651654, 0.00701941,
    0.00753463, 0.00806013, 0.00859378, 0.00913344, 0.00967695,
    0.0102221 , 0.0107667 , 0.01130856, 0.01184551, 0.01237536,
    0.012896  , 0.01340533, 0.01390129, 0.01438189, 0.01484519,
    0.01528933, 0.01571251, 0.01611304, 0.0164893 , 0.01683977,
    0.01716305, 0.01745783, 0.01772293, 0.01795727, 0.01815992,
    0.01833005, 0.01846699, 0.01857018, 0.0186392 , 0.01867379,
    0.01867379, 0.0186392 , 0.01857018, 0.01846699, 0.01833005,
    0.01815992, 0.01795727, 0.01772293, 0.01745783, 0.01716305,
    0.01683977, 0.0164893 , 0.01611304, 0.01571251, 0.01528933,
    0.01484519, 0.01438189, 0.01390129, 0.01340533, 0.012896  ,
    0.01237536, 0.01184551, 0.01130856, 0.0107667 , 0.0102221 ,
    0.00967695, 0.00913344, 0.00859378, 0.00806013, 0.00753463,
    0.00701941, 0.00651654, 0.00602804, 0.00555589, 0.00510197,
    0.00466813, 0.0042561 , 0.00386755, 0.00350404, 0.00316703,
    0.00285788, 0.00257783, 0.00232802, 0.00210944, 0.00192298,
    0.00176939, 0.00164928, 0.00156314, 0.00151131, 0.00149401
};

float firFilterState[FILTER_TAP_NUM] = {0};

// FIR 필터 초기화 함수
void FIRF_Init(void) {
    // 필터 상태 초기화
    for (int i = 0; i < FILTER_TAP_NUM; ++i) {
        firFilterState[i] = 0;
    }
}

// FIR 필터 처리 함수
float FIRF_Process(float input) {

    // 상태 배열 업데이트 (가장 오래된 값을 버리고, 새로운 입력을 배열의 처음에 추가)
    for (int i = FILTER_TAP_NUM - 1; i > 0; --i) {
        firFilterState[i] = firFilterState[i - 1];
    }
    firFilterState[0] = input;

    // FIR 필터 적용 (계수와 상태 배열의 곱의 합 계산)
    float output = 0;
    for (int i = 0; i < FILTER_TAP_NUM; ++i) {
        output += firFilterState[i] * firFilterCoeffs[i];
    }

    return output;
}

static float highPassFilterCoeffs[FILTER_TAP_NUM] = {
    -3.19998918e-06, -3.23629810e-06, -3.34507992e-06, -3.52590537e-06, -3.77806088e-06,
    -4.10055137e-06, -4.49210415e-06, -4.95117398e-06, -5.47594919e-06, -6.06435875e-06,
    -6.71408052e-06, -7.42255037e-06, -8.18697234e-06, -9.00432962e-06, -9.87139651e-06,
    -1.07847511e-05, -1.17407888e-05, -1.27357366e-05, -1.37656680e-05, -1.48265182e-05,
    -1.59141005e-05, -1.70241229e-05, -1.81522044e-05, -1.92938932e-05, -2.04446835e-05,
    -2.16000336e-05, -2.27553839e-05, -2.39061747e-05, -2.50478644e-05, -2.61759473e-05,
    -2.72859713e-05, -2.83735556e-05, -2.94344081e-05, -3.04643419e-05, -3.14592925e-05,
    -3.24153332e-05, -3.33286909e-05, -3.41957610e-05, -3.50131215e-05, -3.57775468e-05,
    -3.64860198e-05, -3.71357447e-05, -3.77241573e-05, -3.82489352e-05, -3.87080075e-05,
    -3.90995624e-05, -3.94220547e-05, -3.96742114e-05, -3.98550376e-05, -3.99638196e-05,
     9.99963200e-01, -3.99638196e-05, -3.98550376e-05, -3.96742114e-05, -3.94220547e-05,
    -3.90995624e-05, -3.87080075e-05, -3.82489352e-05, -3.77241573e-05, -3.71357447e-05,
    -3.64860198e-05, -3.57775468e-05, -3.50131215e-05, -3.41957610e-05, -3.33286909e-05,
    -3.24153332e-05, -3.14592925e-05, -3.04643419e-05, -2.94344081e-05, -2.83735556e-05,
    -2.72859713e-05, -2.61759473e-05, -2.50478644e-05, -2.39061747e-05, -2.27553839e-05,
    -2.16000336e-05, -2.04446835e-05, -1.92938932e-05, -1.81522044e-05, -1.70241229e-05,
    -1.59141005e-05, -1.48265182e-05, -1.37656680e-05, -1.27357366e-05, -1.17407888e-05,
    -1.07847511e-05, -9.87139651e-06, -9.00432962e-06, -8.18697234e-06, -7.42255037e-06,
    -6.71408052e-06, -6.06435875e-06, -5.47594919e-06, -4.95117398e-06, -4.49210415e-06,
    -4.10055137e-06, -3.77806088e-06, -3.52590537e-06, -3.34507992e-06, -3.23629810e-06,
    -3.19998918e-06
};

// 필터 상태
static float highPassFilterState[HPFILTER_TAP_NUM] = {0};

// FIR 고역 통과 필터 초기화
void HighPassFilter_Init(void) {
    for (int i = 0; i < HPFILTER_TAP_NUM; ++i) {
        highPassFilterState[i] = 0;
    }
}

// FIR 고역 통과 필터 처리 함수
float HighPassFilter_Process(float input) {
    for (int i = HPFILTER_TAP_NUM - 1; i > 0; --i) {
        highPassFilterState[i] = highPassFilterState[i - 1];
    }
    highPassFilterState[0] = input;

    float output = 0;
    for (int i = 0; i < HPFILTER_TAP_NUM; ++i) {
        output += highPassFilterState[i] * highPassFilterCoeffs[i];
    }

    return output;
}

float MAF(float new_sample) {
    static float samples[SAMPLE_SIZE] = {0};
    static int index = 0;
    static float sum = 0;
    float average = 0;

    // 이전 합계에서 가장 오래된 샘플 제거
    sum -= samples[index];
    // 새 샘플 추가
    samples[index] = new_sample;
    sum += new_sample;
    // 다음 샘플을 위한 인덱스 업데이트
    index = (index + 1) % SAMPLE_SIZE;

    // 평균 계산
    average = sum / SAMPLE_SIZE;

    return average;
}

float EWMAF(float new_measurement, float prev_ewma, float alpha) {
    return alpha * new_measurement + (1 - alpha) * prev_ewma;
}

void KMF_Init(KMF *kf, float init_estimate, float init_error_estimate, float error_measure) {
    kf->estimate = init_estimate;
    kf->error_estimate = init_error_estimate;
    kf->error_measure = error_measure;
}

void KMF_Update(KMF *kf, float measurement) {
    kf->kalman_gain = kf->error_estimate / (kf->error_estimate + kf->error_measure);
    kf->estimate = kf->estimate + kf->kalman_gain * (measurement - kf->estimate);
    kf->error_estimate = (1 - kf->kalman_gain) * kf->error_estimate;
}
