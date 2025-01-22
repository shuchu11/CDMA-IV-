#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>


#define N 16       // ñ�W�i�Ϊ��ס]16 �Ӵ����^
#define MAX_ITER 1000 // ���N����
#define alpha 0.1   // �ڤɾl���o�i�����u���Y��
#define Tc 1.0     // �����g��
#define p 0.333333    // ���i�Ҿ��v
#define oversampling (2*N)     // ���˲v�]�C�Ӵ��������q�ơ^
#define FILTER_SIZE (16*32+1)  // �ھ� MATLAB �p�⪺�Y�Ƽƶq�]�m�j�p

#define PI 3.141592653589793



// �b�����i�� psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
		return cos(PI * t / Tc);
    }
    return 0.0;
}
// ��l��ñ�W�i�� c_k^(0)(t)
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length) {
    int extended_length = 2 * length; // 2N
    double step = Tc / (oversampling); // �T�O���B�I�B��
    int total_samples = N * oversampling+1;          // ���( extended_length * oversampling )
	//printf("Tc :%f , oversampling : %d , step : %f ",Tc,oversampling,step);
	// ��l���A�G�H����� +1 �� -1
    s_k[0] = (rand() % 2 ? 1 : -1);
    // �ͦ����i�ҧǦC
/*
    for (int i = 1; i < extended_length; i++) {
        double rand_prob = (double)(rand() % 100) / 100; // �H�����v [0, 1)
        if (s_k[i - 1] == 1) {
            s_k[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 �����v p
        }
        else {
            s_k[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 �����v p
        }
    }
*/

FILE *file = fopen("chain_data.txt", "r");

    if (file == NULL) {
        printf("Error: Unable to open file 'chain_data.txt'.\n");
        return ;
    }

    // �q�ɮפ�Ū���ƾڨæs�J s_k[] �}�C
    for (int i = 0; i < 2*N; i++) {
        if (fscanf(file, "%lf", &s_k[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
            return ;
        }
    }
    fclose(file);


	// c_k^(0)(t)
    for (int i = 0; i < total_samples; i++) {
        double t = i * step; // ��e�ɶ�
        double complex value = 0.0+0.0*I;

        // �ھڤ��� (6) �p��i��c_k
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * PI / 2.0; // (-1)^k n ���ۦ찾��
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
			double psi_value = t - n * Tc / 2.0 ;

            //value += phase * s_k[n] * psi_value;
			value += phase * s_k[n] * psi(psi_value);
            //printf("%f \n ",psi(psi_value));
			//printf(" t: %f , nTc/2 : %f , theda : %f , psi : %f \n", t , n * Tc / 2.0 , psi_value , psi(psi_value) );//�d��phase���Ҧ������ƭ�
			//printf(" phase : %.1f %.1fj ,s_k[%d] : %.1f , psi: %.3f \n",creal(phase),cimag(phase),n,s_k[n],psi_value);
        }

        k_z[i] = value; // �O�s�i��
    }
    /*for (int i=0;i<total_samples;i++)
        printf("%f %fi\n",creal(k_z[i]),cimag(k_z[i]));*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ���T���W��
void normalize(double complex *k_z, int length) {
	double max ;
	double k_z_amp[length] ;
	double k_z_arg[length];
    for (int i = 0; i < length; i++) {
        k_z_amp[i] = sqrt(creal(k_z[i]) * creal(k_z[i]) + cimag(k_z[i]) * cimag(k_z[i]));     //�s�� c_k���T
		k_z_arg[i] = atan2(cimag(k_z[i]), creal(k_z[i]));        //�s�� c_k�ۦ�
    }
	max = k_z_amp[0];
	for (int i = 0; i < length; i++){
		if(k_z_amp[i] >= max )
			max = k_z_amp[i];
	}
	for (int i = 0; i < length; i++){
		k_z_amp[i] = k_z_amp[i] / max;
		k_z[i] = k_z_amp[i] * cexp(I * k_z_arg[i]);
    }

}

// �ڤɾl���o�i���ͦ�
void generate_rrc_filter(double *filter, int length) {
    FILE *file = fopen("filter_data.txt", "r");

    if (file == NULL) {
        printf("Error: Unable to open file 'filter_data.txt'.\n");
        return ;
    }

    // �q�ɮפ�Ū���ƾڨæs�J filter[] �}�C
    for (int i = 0; i < FILTER_SIZE; i++) {
        if (fscanf(file, "%lf", &filter[i]) != 1) {
            printf("Error: Failed to read data at index %d.\n", i);
            fclose(file);
            return ;
        }
    }
    fclose(file);

    double filter_copy[length];
    for(int i=0;i<length;i++){
        filter_copy[i] = filter[i];
    }
    for(int i=0;i<N*N;i++){
        filter_copy[i] = filter[i+N*N+1];
        filter_copy[i+N*N+1] = filter[i];
    }
    for(int i=0;i<length;i++) filter[i] = filter_copy[i];


/*
    // ���տ�X����Ū�����ƾ�
    printf("Filter data (first 10 values):\n");
    for (int i = 0; i < 10; i++) {
        printf("filter[%d] = %.8e\n", i, filter[i]);
    }
    */
}

// ��ơG�p��ֿn�ۦ�
void calculate_unwrapped_phase(const double real[], const double imag[], double unwrapped_phase[], int length) {
    double prev_phase = 0.0; // ��l�ۦ�
    double cumulative_phase = 0.0;

    for (int i = 0; i < length; i++) {
        // �p���e�I�����ɬۦ�
        double current_phase = atan2(imag[i], real[i]);

        if (i > 0) {
            // �p��۾F�I���ۦ�t
            double phase_diff = current_phase - prev_phase;

            // �ץ��ۦ����
            if (phase_diff > PI) {
                phase_diff -= 2.0 * PI;
            } else if (phase_diff < -PI) {
                phase_diff += 2.0 * PI;
            }

            // ��s�ֿn�ۦ�
            cumulative_phase += phase_diff;
        }
		else {
            // ���Ĥ@�I�A�ֿn�ۦ쵥�����ɬۦ�
            cumulative_phase = current_phase;
        }

        // �s�x�ֿn�ۦ�
        unwrapped_phase[i] = cumulative_phase;
        prev_phase = current_phase;
    }
    //unwrapped_phase[length-1] = unwrapped_phase[length-2];
}

void circular_convolution(double complex* x1, double* x2, int len_x1, int len_x2, int n, double complex* y) {
    int len_conv = len_x1 + len_x2 - 1;
    double complex linear_result[len_conv];

    // ��l�� linear_result
    for (int i = 0; i < len_conv; i++) {
        linear_result[i] = 0.0 + 0.0 * I;
    }

    // �p��u�ʱ��n
    for (int i = 0; i < len_x1; i++) {
        for (int j = 0; j < len_x2; j++) {
            linear_result[i + j] += x1[i] * x2[j];
        }
    }

    // ��l�� y
    for (int i = 0; i < n; i++) {
        y[i] = 0.0 + 0.0 * I;
    }

    // �p��`�����n
    for (int i = 0; i < n; i++) {
        if (i < len_conv - n) {
            y[i] = linear_result[n + i];
        }
        y[i] += linear_result[i];
    }
}
// �D�y�{�A�ëO�s�S�w���N���ƪ��ƾ�
void generate_ce_waveform(double complex *k_z, int length, const char *output_file ) {
    double filter[length];
    generate_rrc_filter(filter, length ); // �ͦ��ڤɧE���o�i��

    //for(int i=0; i < length ;i++)
    //    printf("%f\n",filter[i]);

    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Iteration,Time,Phase\n"); // �g�J���Y

    int save_index = 0;
    int save_iters[] = {0, 1, 10, 100, 1000 }; // �n�O�s�����N����
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]);
    for (int iter = 0; iter <= MAX_ITER; iter++) {
        if(iter!=0){

		// �o�i�]²�ƪ��`�����n�^
		//for (int i = 0; i < length; i++) {
          /*  double complex x1[] = {1 , 2 , 3 , 4 };
            double x2[] = {4, 3, 2, 1};
          */
            double complex y[length]; // �אּ double complex ���A�����G�}�C
            int len_x1 = length;
            int len_x2 = sizeof(filter) / sizeof(filter[0]);
            int n = length; // ���w�`�����n������
            // �ե� circular_convolution �禡
            circular_convolution(k_z, filter, len_x1, len_x2, n, y);

      /*      // ��X���G
            printf("Circular Convolution Result:\n");
            for (int i = 0; i < n; i++) {
                printf("%.2f ", y[i]);
            }
            printf("\n");
      */
            for(int i=0;i<length;i++)
                k_z[i] = y[i];

	/*		double complex temp = 0.0;
			for (int j = 0; j < length; j++) {
				temp += k_z[j] * filter[(i - j + length) % length];
			}
			k_z[i] = temp;
		}
    */
		normalize(k_z, length); // ���W�ƴT��
        }
        //for (int i=0;i< N * oversampling;i++)
        //printf("%f %fi\n",creal(k_z[i]),cimag(k_z[i]));
		//printf("%d",num_saves);
        // �O�s�ƾڡ]�Ȧb���w���N���ơ^
        if (save_index < num_saves && iter == save_iters[save_index]) {
			//printf("save_index : %d , iter : %d\n",save_index,iter);       //���� && iter == save_iters[save_index]
			double real[length];
			double imag[length];
			double unwrapped_phase[length]; // �s�x�ֿn�ۦ�
			for (int i = 0; i < length; i++) {
				real[i] = creal(k_z[i]);
				imag[i] = cimag(k_z[i]);
			}
			// �p��ֿn�ۦ�
			calculate_unwrapped_phase(real, imag, unwrapped_phase, length);

            //printf("This is length : %d \n",length);
            for (int i = 0; i < length; i++) {
                //printf("This is (i * Tc) /oversampling : %f \n" , i*(Tc/ oversampling) );
                fprintf(file, "%d,%f,%f\n", iter, i*(Tc/ oversampling), unwrapped_phase[i]/PI ) ; // �p��ۦ�
            }
            save_index++;
        }
    }
    fclose(file);
}


int main() {
    int total_samples = N * oversampling + 1; // �`�����I��       // ���( 2 * N * oversampling )
    double s_k[2 * N]; // ��l�����ǦC
    double complex k_z[total_samples]; // ñ�W�i��


    srand((unsigned int)time(NULL)); // ��l���H���ƺؤl

    initialize_signature_waveform(k_z,s_k, 3 ,N) ;// ��l�ƪi�Ρ]���i�ҹL�{�^
    generate_ce_waveform(k_z, total_samples,"ce_waveform_data.csv"); // �O�s�ƾ�

    printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
