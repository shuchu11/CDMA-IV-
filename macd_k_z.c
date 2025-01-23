#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>


#define N 16       // ñ�W�i�Ϊ��ס]16 �Ӵ����^
#define k 1        // user index  k�Y�����ֿn�ۦ��Ͷլ����A�Ϥ��Ͷլ��t
#define MAX_ITER 1000 // �̤j���N����
#define alpha 0.1   // �ڤɾl���o�i�����u���Y��
#define Tc 1.0     // �����g��
#define p 0.333333    // ���i�Ҿ��v
#define oversampling (32)    // ���˲v(Hz)  [�ɶ����Tc]
#define FILTER_SIZE (oversampling*N+1)  // �ڤɾl���o�i���}�C�j�p
#define PI 3.141592653589793
//////////////////////////////// ����(6)�l�禡 //////////////////////////////////////////////////////////////////////////////////////
// �b�����i�� psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
		return cos(PI * t / Tc);
    }
    return 0.0;
}
//////////////////////////////// ����(6)�D�禡 ///////////////////////////////////////////
// ��l��ñ�W�i�� c_k^(0)(t)
void initialize_signature_waveform(double complex *c_k, double *s_k, int user_index, int length) {
    int extended_length = 2 * length; // 2N
    double step = Tc / (oversampling); // �T�O���B�I�B��
    int total_samples = N * oversampling+1;          // ���( extended_length * oversampling )

	// ��l���A�G�H����� +1 �� -1
    s_k[0] = (rand() % 2 ? 1 : -1);
    // �ͦ����i���H���ǦC

/*  // ��k�@ : �b���[���������ͦ��ǦC
    for (int i = 1; i < extended_length; i++) {
        double rand_prob = (double)(rand() % 100) / 100; // �H�����v [0, 1)
        if (c_k[i - 1] == 1) {
            c_k[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 �����v p
        }
        else {
            c_k[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 �����v p
        }
    }
*/

    // ��k�G : �� Mk_chain.m �ͦ��ǦC�A�A�ɤJ���禡
    FILE *file = fopen("chain_data.txt", "r");

    if (file == NULL) {
        printf("Error: Unable to open file 'chain_data.txt'.\n");
    }
    // �q�ɮפ�Ū���ƾڨæs�J s_k[] �}�C
    for (int i = 0; i < 2*N; i++) {
        if (fscanf(file, "%lf", &s_k[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
        }
    }
    fclose(file);

	// �ͦ� C_k^(0)(t)
    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // ��e�ɶ�
        double complex value = 0.0+0.0*I;
        // �ھڤ��� (6) �p��i��c_k
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * PI / 2.0; // (-1)^k n ���ۦ찾��
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
			double psi_value = t - n * Tc / 2.0 ;

			value += phase * s_k[n] * psi(psi_value);
        }
        c_k[i] = value; // �O�s�i��
    }
}
//////////////////////////////// ����(8)�l�禡 //////////////////////////////////////////////////////////////////////////////////////
// ���T���W��   [ ����(9) ]
void normalize(double complex *c_k, int length) {
	double max ; // �x�s��e�̤j���T
	double c_k_amp[length] ; // Ck ���T
	double c_k_arg[length];  // Ck �ۦ�

	// ���� Ck ���T�M�ۦ�
    for (int i = 0; i < length; i++) {
        c_k_amp[i] = sqrt(creal(c_k[i]) * creal(c_k[i]) + cimag(c_k[i]) * cimag(c_k[i]));     //�s�� c_k���T
		c_k_arg[i] = atan2(cimag(c_k[i]), creal(c_k[i]));        //�s�� c_k�ۦ�
    }
	max = c_k_amp[0]; // �N�̤j�T�ȩw���Ĥ@�Ӥ���
	// �M��̤j�T��
	for (int i = 0; i < length; i++){
		if(c_k_amp[i] >= max )
			max = c_k_amp[i];
	}
	// �N�w���W�ƪ� c_k ���N��}�C
	for (int i = 0; i < length; i++){
		c_k_amp[i] = c_k_amp[i] / max;
		c_k[i] = c_k_amp[i] * cexp(I * c_k_arg[i]);
    }
}

// �ڤɾl���o�i���ͦ�
void generate_rrc_filter(double *filter, int length) {
    FILE *file = fopen("filter_data.txt", "r");
    // �T�{��Ū���ɮצs�b
    if (file == NULL) {
        printf("Error: Unable to open file 'filter_data.txt'.\n");
    }
    // �q�ɮפ�Ū���ƾڨæs�J filter[] �}�C
    for (int i = 0; i < FILTER_SIZE; i++) {
        if (fscanf(file, "%lf", &filter[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
        }
    }
    fclose(file);

    // �����o�i���s�J�x�}���٩�y�b�S�ʮ����A�N�o�i���e�b�q�M��b�q����
    double filter_copy[length];
    for(int i=0;i<length;i++){
        filter_copy[i] = filter[i];  //�ƻs���o�i���}�C
    }
    //�o�i���e�b�q�M��b�q����
    for(int i=0;i<N*N;i++){
        filter_copy[i] = filter[i+N*N+1];
        filter_copy[i+N*N+1] = filter[i];
    }
    // �����o�i���ק�
    for(int i=0;i<length;i++) filter[i] = filter_copy[i];
}

// ��ơG�p��ֿn�ۦ�
void calculate_unwrapped_phase(double complex *c_k, double *unwrapped_phase, int length) {
    double real[length];
    double imag[length];
    double prev_phase = 0.0; // ��l�ۦ�
    double cumulative_phase = 0.0; //��e�w�ֿn�ۦ�

    // �����ֿn�ۦ줽���ݭn����
    for (int i = 0; i < length; i++) {
        real[i] = creal(c_k[i]);
        imag[i] = cimag(c_k[i]);
    }

    for (int i = 0; i < length; i++) {
        // �p���e�I�����ɬۦ�
        double current_phase = atan2(imag[i], real[i]);

        if (i > 0) {
            // �p��۾F�I���ۦ�t
            double phase_diff = current_phase - prev_phase;

            // �ץ��ۦ����
            if (phase_diff > PI)
                phase_diff -= 2.0 * PI;
            else if (phase_diff < -PI)
                phase_diff += 2.0 * PI;
            // ��s�ֿn�ۦ�
            cumulative_phase += phase_diff;
        }
		else
            cumulative_phase = current_phase; // ���Ĥ@�I�A�ֿn�ۦ쵥�����ɬۦ�
        // �s�x�ֿn�ۦ�
        unwrapped_phase[i] = cumulative_phase;
        prev_phase = current_phase;
    }
    //unwrapped_phase[length-1] = unwrapped_phase[length-2];
}

// �p��`�����n
void circular_convolution(double complex* x1, double* x2, int len_x1 ,int len_x2 , int n, double complex* y) {
    int len_conv = len_x1 + len_x2 - 1;  // �u�ʱ��n����
    double complex linear_result[len_conv]; // �u�ʱ��n��

    // ���p��u�ʱ��n
    for (int i = 0; i < len_x1; i++) {
        for (int j = 0; j < len_x2; j++) {
            linear_result[i + j] += x1[i] * x2[j];
        }
    }
    // �p��`�����n
    for (int i = 0; i < n; i++) {
        if (i < len_conv - n) {
            y[i] = linear_result[n + i];
        }
        y[i] += linear_result[i];
    }
}
//////////////////////////////// ����(8)�D�禡 //////////////////////////////////////////////////////////////////////////////////////
void generate_ce_waveform(double complex *c_k, int length, const char *output_file ) {
    double filter[length];  // �o�i��
    generate_rrc_filter(filter, length ); // �ͦ��ڤɾl���o�i��

    // �}�ұN�n�s�J�ֿn�ۦ쪺.csv��
    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Iteration,Time,Phase\n"); // �g�J���Y

    int save_iters[] = {0, 1, 10, 100, 1000 }; // �n�O�s�����w���N
    int save_index = 0; //���w���N�w���榸��
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]); // ���w���N�Ӽ�

    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // iter = 0 ���i��`�����n
        if(iter!=0){
            double complex y[length]; // �`�����n���G��
            int n = length; // ���w�`�����n������

            // �ե� circular_convolution �禡
            circular_convolution(c_k , filter, length , FILTER_SIZE , n , y);

            // �sC_k���`�����n���G��
            for(int i=0;i<length;i++)
                c_k[i] = y[i];

            // ���W�ƴT��
            normalize(c_k, length);
        }

        // �O�s�ƾڡ]�Ȧb���w���N���ơ^
        if (save_index < num_saves && iter == save_iters[save_index]) {

			double unwrapped_phase[length]; // �s�x�ֿn�ۦ�

			// �p��ֿn�ۦ�
			calculate_unwrapped_phase(c_k, unwrapped_phase, length);

			// �N�ֿn�ۦ�ƾڿ�X�A�N��� ce_fig1_.m ø��
            for (int i = 0; i < length; i++) {
                fprintf(file, "%d,%f,%f\n", iter, i*(Tc/ oversampling), unwrapped_phase[i]/PI ) ; // �p��ۦ�
            }
            save_index++;
        }
    }
    fclose(file);
}

//////////////////////////////// �D�{�� //////////////////////////////////////////////////////////////////////////////////////

int main() {
    int total_samples = N * oversampling + 1; // �`�����I��
    double s_k[2 * N]; // ��l�����ǦC
    double complex c_k[total_samples]; // ñ�W�i��

    srand((unsigned int)time(NULL)); // ��l���H���ƺؤl

    initialize_signature_waveform(c_k,s_k, k ,N) ;// ��l�ƪi�Ρ]���i�ҹL�{�^
    generate_ce_waveform(c_k, total_samples,"ce_waveform_data.csv"); // �O�s�ƾ�

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
