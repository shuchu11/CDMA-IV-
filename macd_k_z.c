#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>


#define N 16       // 簽名波形長度（16 個晶片）
#define k 1        // user index  k若為正累積相位趨勢為正，反之趨勢為負
#define MAX_ITER 1000 // 最大迭代次數
#define alpha 0.1   // 根升餘弦濾波器的滾降係數
#define Tc 1.0     // 晶片週期
#define p 0.333333    // 馬可夫機率
#define oversampling (32)    // 采樣率(Hz)  [時間單位Tc]
#define FILTER_SIZE (oversampling*N+1)  // 根升餘弦濾波器陣列大小
#define PI 3.141592653589793
//////////////////////////////// 公式(6)子函式 //////////////////////////////////////////////////////////////////////////////////////
// 半正弦波形 psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
		return cos(PI * t / Tc);
    }
    return 0.0;
}
//////////////////////////////// 公式(6)主函式 ///////////////////////////////////////////
// 初始化簽名波形 c_k^(0)(t)
void initialize_signature_waveform(double complex *c_k, double *s_k, int user_index, int length) {
    int extended_length = 2 * length; // 2N
    double step = Tc / (oversampling); // 確保為浮點運算
    int total_samples = N * oversampling+1;          // 更動( extended_length * oversampling )

	// 初始狀態：隨機選擇 +1 或 -1
    s_k[0] = (rand() % 2 ? 1 : -1);
    // 生成馬可夫隨機序列

/*  // 方法一 : 在此涵式內直接生成序列
    for (int i = 1; i < extended_length; i++) {
        double rand_prob = (double)(rand() % 100) / 100; // 隨機機率 [0, 1)
        if (c_k[i - 1] == 1) {
            c_k[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 的機率 p
        }
        else {
            c_k[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 的機率 p
        }
    }
*/

    // 方法二 : 由 Mk_chain.m 生成序列，再導入此函式
    FILE *file = fopen("chain_data.txt", "r");

    if (file == NULL) {
        printf("Error: Unable to open file 'chain_data.txt'.\n");
    }
    // 從檔案中讀取數據並存入 s_k[] 陣列
    for (int i = 0; i < 2*N; i++) {
        if (fscanf(file, "%lf", &s_k[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
        }
    }
    fclose(file);

	// 生成 C_k^(0)(t)
    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // 當前時間
        double complex value = 0.0+0.0*I;
        // 根據公式 (6) 計算波形c_k
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * PI / 2.0; // (-1)^k n 的相位偏移
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
			double psi_value = t - n * Tc / 2.0 ;

			value += phase * s_k[n] * psi(psi_value);
        }
        c_k[i] = value; // 保存波形
    }
}
//////////////////////////////// 公式(8)子函式 //////////////////////////////////////////////////////////////////////////////////////
// 振幅正規化   [ 公式(9) ]
void normalize(double complex *c_k, int length) {
	double max ; // 儲存當前最大振幅
	double c_k_amp[length] ; // Ck 振幅
	double c_k_arg[length];  // Ck 相位

	// 提取 Ck 振幅和相位
    for (int i = 0; i < length; i++) {
        c_k_amp[i] = sqrt(creal(c_k[i]) * creal(c_k[i]) + cimag(c_k[i]) * cimag(c_k[i]));     //存取 c_k振幅
		c_k_arg[i] = atan2(cimag(c_k[i]), creal(c_k[i]));        //存取 c_k相位
    }
	max = c_k_amp[0]; // 將最大幅值定為第一個元素
	// 尋找最大幅值
	for (int i = 0; i < length; i++){
		if(c_k_amp[i] >= max )
			max = c_k_amp[i];
	}
	// 將已正規化的 c_k 取代原陣列
	for (int i = 0; i < length; i++){
		c_k_amp[i] = c_k_amp[i] / max;
		c_k[i] = c_k_amp[i] * cexp(I * c_k_arg[i]);
    }
}

// 根升餘弦濾波器生成
void generate_rrc_filter(double *filter, int length) {
    FILE *file = fopen("filter_data.txt", "r");
    // 確認欲讀取檔案存在
    if (file == NULL) {
        printf("Error: Unable to open file 'filter_data.txt'.\n");
    }
    // 從檔案中讀取數據並存入 filter[] 陣列
    for (int i = 0; i < FILTER_SIZE; i++) {
        if (fscanf(file, "%lf", &filter[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
        }
    }
    fclose(file);

    // 彌補濾波器存入矩陣後對稱於y軸特性消失，將濾波器前半段和後半段互換
    double filter_copy[length];
    for(int i=0;i<length;i++){
        filter_copy[i] = filter[i];  //複製原濾波器陣列
    }
    //濾波器前半段和後半段互換
    for(int i=0;i<N*N;i++){
        filter_copy[i] = filter[i+N*N+1];
        filter_copy[i+N*N+1] = filter[i];
    }
    // 完成濾波器修改
    for(int i=0;i<length;i++) filter[i] = filter_copy[i];
}

// 函數：計算累積相位
void calculate_unwrapped_phase(double complex *c_k, double *unwrapped_phase, int length) {
    double real[length];
    double imag[length];
    double prev_phase = 0.0; // 初始相位
    double cumulative_phase = 0.0; //當前已累積相位

    // 分離累積相位公式需要的值
    for (int i = 0; i < length; i++) {
        real[i] = creal(c_k[i]);
        imag[i] = cimag(c_k[i]);
    }

    for (int i = 0; i < length; i++) {
        // 計算當前點的瞬時相位
        double current_phase = atan2(imag[i], real[i]);

        if (i > 0) {
            // 計算相鄰點的相位差
            double phase_diff = current_phase - prev_phase;

            // 修正相位跳變
            if (phase_diff > PI)
                phase_diff -= 2.0 * PI;
            else if (phase_diff < -PI)
                phase_diff += 2.0 * PI;
            // 更新累積相位
            cumulative_phase += phase_diff;
        }
		else
            cumulative_phase = current_phase; // 對於第一點，累積相位等於瞬時相位
        // 存儲累積相位
        unwrapped_phase[i] = cumulative_phase;
        prev_phase = current_phase;
    }
    //unwrapped_phase[length-1] = unwrapped_phase[length-2];
}

// 計算循環捲積
void circular_convolution(double complex* x1, double* x2, int len_x1 ,int len_x2 , int n, double complex* y) {
    int len_conv = len_x1 + len_x2 - 1;  // 線性捲積長度
    double complex linear_result[len_conv]; // 線性捲積值

    // 先計算線性捲積
    for (int i = 0; i < len_x1; i++) {
        for (int j = 0; j < len_x2; j++) {
            linear_result[i + j] += x1[i] * x2[j];
        }
    }
    // 計算循環捲積
    for (int i = 0; i < n; i++) {
        if (i < len_conv - n) {
            y[i] = linear_result[n + i];
        }
        y[i] += linear_result[i];
    }
}
//////////////////////////////// 公式(8)主函式 //////////////////////////////////////////////////////////////////////////////////////
void generate_ce_waveform(double complex *c_k, int length, const char *output_file ) {
    double filter[length];  // 濾波器
    generate_rrc_filter(filter, length ); // 生成根升餘弦濾波器

    // 開啟將要存入累積相位的.csv檔
    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Iteration,Time,Phase\n"); // 寫入表頭

    int save_iters[] = {0, 1, 10, 100, 1000 }; // 要保存的指定迭代
    int save_index = 0; //指定迭代已執行次數
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]); // 指定迭代個數

    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // iter = 0 不進行循環捲積
        if(iter!=0){
            double complex y[length]; // 循環捲積結果值
            int n = length; // 指定循環捲積的長度

            // 調用 circular_convolution 函式
            circular_convolution(c_k , filter, length , FILTER_SIZE , n , y);

            // 新C_k為循環捲積結果值
            for(int i=0;i<length;i++)
                c_k[i] = y[i];

            // 正規化幅值
            normalize(c_k, length);
        }

        // 保存數據（僅在指定迭代次數）
        if (save_index < num_saves && iter == save_iters[save_index]) {

			double unwrapped_phase[length]; // 存儲累積相位

			// 計算累積相位
			calculate_unwrapped_phase(c_k, unwrapped_phase, length);

			// 將累積相位數據輸出，將交由 ce_fig1_.m 繪圖
            for (int i = 0; i < length; i++) {
                fprintf(file, "%d,%f,%f\n", iter, i*(Tc/ oversampling), unwrapped_phase[i]/PI ) ; // 計算相位
            }
            save_index++;
        }
    }
    fclose(file);
}

//////////////////////////////// 主程式 //////////////////////////////////////////////////////////////////////////////////////

int main() {
    int total_samples = N * oversampling + 1; // 總采樣點數
    double s_k[2 * N]; // 初始晶片序列
    double complex c_k[total_samples]; // 簽名波形

    srand((unsigned int)time(NULL)); // 初始化隨機數種子

    initialize_signature_waveform(c_k,s_k, k ,N) ;// 初始化波形（馬可夫過程）
    generate_ce_waveform(c_k, total_samples,"ce_waveform_data.csv"); // 保存數據

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
