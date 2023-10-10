#include <stdio.h>
#include "math.h"
#include <complex.h>
#include <fftw3.h>
#include "SSCEq_2D_full.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

void gethostname();

time_t s_time;
time_t f_time;


double Random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}

void CheckFile(FILE* file_name)
{
	if (file_name == NULL)
	{
		printf("Cannot open or find the file \n");
		exit(0);
	}
}

void Set_new_data()
{
	// C (k, t) 

/*	for (int kx = 0; kx < Nx; kx++)
	{
		for (int ky = 0; ky < Ny; ky++)
		{
			C1_k[kx][ky] = (1e-12 + I * 0.0) * cexp(I * Random(0.0, 2.0 * M_PI));
			C2_k[kx][ky] = (1e-12 + I * 0.0) * cexp(I * Random(0.0, 2.0 * M_PI));
		}
	}*/

	//C1_k[10][1] = 0.001 + I * 0.0;
	//C1_k[25][2] = 0.001 + I * 0.0;
	//C1_k[50][3] = 0.01 + I * 0.0;
	//C1_k[50][Ny - 1] = 0.01 + I * 0.0;
	C1_k[100][0] = 0.0025;// * (1.0 + I);
	C1_k[130][3] = 0.0001;// * (1.0 + I);
	C1_k[70][Ny - 3] = 0.0001;// * (1.0 + I);
	//C2_k[Nx - 10][1] = 0.01 + I * 0.0;
	//C2_k[Nx - 25][2] = 0.01 + I * 0.0;
	//C2_k[Nx - 30][0] = 0.00275 * (1.0 + I);
	//C2_k[Nx - 50][0] = 0.05 + I * 0.0;
	//C2_k[Nx - 110][Ny - 3] = 0.01 + I * 0.0;

/*	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = Nx/2; kx < Nx; kx++)
		{
			C2_k[kx][ky] = C1_k[kx - Nx/2][ky];
		}
	}*/

	printf("%s\n", "New data was successfully set");
}

void Load_data()
{
	int i,j;

	pfile = fopen ("./Restart_Files/0_restart","rb");
	CheckFile(pfile);
	binary_number = fread(&b_k_restart[0], sizeof(fftw_complex), Nx * Ny/2, pfile);
	fclose(pfile);

	for (int ky = 0; ky < Ny/4; ky++)
	{
		for (int kx = 0; kx < Nx/4; kx++)
		{
			C1_k[kx][ky] = sqrt(K[kx][ky]) * b_k_restart[kx + ky * Nx];
		}
	}

	for (int ky = Ny - Ny/4 + 1; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/4; kx++)
		{
			C1_k[kx][ky] = sqrt(K[kx][ky]) * b_k_restart[kx + (ky - Ny/2) * Nx];
		}
	}

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = Nx/2; kx < Nx; kx++)
		{
			C1_k[kx][ky] = 0.0 + I * 0.0;
		}
	}

/*	pfile = fopen ("./Restart_Files/C1_k_data_7","rb");
	CheckFile(pfile);
	binary_number = fread(&C1_k[0][0], sizeof(fftw_complex), Nx * Ny, pfile);
	fclose(pfile);

	pfile = fopen ("./Restart_Files/C2_k_data_7","rb");
	CheckFile(pfile);
	binary_number = fread(&C2_k[0][0], sizeof(fftw_complex), Nx * Ny, pfile);
	fclose(pfile);*/

/*		for (int kx = 0; kx < Nx; kx++)
		{
			for (int ky = 0; ky < Ny; ky++)
			{
				C1_kxky[kx][ky] = C1_kykx[ky][kx];
				C1_kxky[kx][ky] = C1_kxky[kx][ky] * cexp(I * kx * Dk_x * 2500);	
			}
		}

		for ( int kx = 1; kx < Nx/2; kx++)
		{
			C2_kxky[Nx - kx][0] = C1_kxky[kx][0];				
		} 
*/
/*		for ( int kx = 0; kx < Nx/2; kx++)
		{
			C1_kxky[kx][0] = 0.0 + I * 0.0;
		} */


	printf("%s\n", "The data loaded successfully");
}

void Set_data()
{
	//  W(kx,ky) massive initializing

	t = 0.0;
	hx = (bx - ax)/Nx;
	hy = (by - ay)/Ny;
	Dk_x = 2 * M_PI/ (bx - ax);
	Dk_y = 2 * M_PI/ (by - ay);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			C1_k[kx][ky] = 0.0 + I * 0.0;
			C2_k[kx][ky] = 0.0 + I * 0.0;
		} 
	}

	for (int kx = 0; kx < Nx/2; kx++)
	{
		K[kx][Ny/2] = 0.0 + I * 0.0;
		for (int ky = 0; ky < Ny/2; ky++)
		{
			K[kx][ky] = sqrt(kx * kx * Dk_x * Dk_x + ky * ky * Dk_y * Dk_y);
			W[kx][ky] = sqrt(g * K[kx][ky]);
		}

		for (int ky = Ny/2 + 1; ky < Ny; ky++)
		{
			K[kx][ky] = sqrt(kx * kx * Dk_x * Dk_x + (ky - Ny) * (ky - Ny) * Dk_y * Dk_y);
			W[kx][ky] = sqrt(g * K[kx][ky]);
		}
	}

	for (int kx = Nx/2 + 1; kx < Nx; kx++)
	{
		K[kx][Ny/2] = 0.0 + I * 0.0;
		for (int ky = 0; ky < Ny/2; ky++)
		{
			K[0][ky] = 0.0 + I * 0.0;
			K[Nx/2][ky] = 0.0 + I * 0.0;
			K[kx][ky] = sqrt((kx - Nx) * (kx - Nx) * Dk_x * Dk_x + ky * ky * Dk_y * Dk_y);
			W[kx][ky] = sqrt(g * K[kx][ky]);
		}

		for (int ky = Ny/2; ky < Ny; ky++)
		{
			K[0][ky] = 0.0 + I * 0.0;
			K[Nx/2][ky] = 0.0 + I * 0.0;
			K[kx][ky] = sqrt((kx - Nx) * (kx - Nx) * Dk_x * Dk_x + (ky - Ny) * (ky - Ny) * Dk_y * Dk_y);
			W[kx][ky] = sqrt(g * K[kx][ky]);
		}
	}

	//Set_new_data();
	
	Load_data();

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			C2_k[Nx/2 + kx][ky] = C1_k[Nx/2 - kx][ky];
		}
	}

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
			{
				C1_k_old[kx][ky] = C1_k[kx][ky];
				C2_k_old[kx][ky] = C2_k[kx][ky];
			} 
	}
}


void Plans_Creating()
{

	C1k_to_C1x_2D = fftw_plan_dft_2d(Nx,Ny,&C1_k[0][0], &C1_x[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	C2k_to_C2x_2D = fftw_plan_dft_2d(Nx,Ny,&C2_k[0][0], &C2_x[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	dC1_k_to_dC1_x_2D = fftw_plan_dft_2d(Nx,Ny,&dC1_k[0][0], &dC1_x[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	dC2_k_to_dC2_x_2D = fftw_plan_dft_2d(Nx,Ny,&dC2_k[0][0], &dC2_x[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	Eta_1_forward_2D = fftw_plan_dft_2d(Nx,Ny,&Eta_1[0][0], &Eta_1[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	Eta_1_backward_2D = fftw_plan_dft_2d(Nx,Ny,&Eta_1[0][0], &Eta_1[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	Eta_2_to_Eta_2_forward_2D = fftw_plan_dft_2d(Nx,Ny,&Eta_2[0][0], &Eta_2[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	Eta_2_to_Eta_2_backward_2D = fftw_plan_dft_2d(Nx,Ny,&Eta_2[0][0], &Eta_2[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	NL_1p_to_NL_1p_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_1p[0][0], &NL_1p[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	NL_1m_to_NL_1m_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_1m[0][0], &NL_1m[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

	NL_2p_to_NL_2p_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_2p[0][0], &NL_2p[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	NL_2m_to_NL_2m_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_2m[0][0], &NL_2m[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

	NL_2p_to_NL_2p_backward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_2p[0][0], &NL_2p[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	NL_2m_to_NL_2m_backward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_2m[0][0], &NL_2m[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	NL_3p_to_NL_3p_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_3p[0][0], &NL_3p[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	NL_3m_to_NL_3m_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_3m[0][0], &NL_3m[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

	NL_4p_to_NL_4p_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_4p[0][0], &NL_4p[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	NL_4m_to_NL_4m_forward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_4m[0][0], &NL_4m[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

	NL_4p_to_NL_4p_backward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_4p[0][0], &NL_4p[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	NL_4m_to_NL_4m_backward_2D = fftw_plan_dft_2d(Nx, Ny, &NL_4m[0][0], &NL_4m[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	K_abs_C1x_forward_2D = fftw_plan_dft_2d(Nx, Ny, &K_abs_C1[0][0], &K_abs_C1[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
	K_abs_C2x_forward_2D = fftw_plan_dft_2d(Nx, Ny, &K_abs_C2[0][0], &K_abs_C2[0][0], FFTW_FORWARD, FFTW_ESTIMATE);

	K_abs_C1x_backward_2D = fftw_plan_dft_2d(Nx, Ny, &K_abs_C1[0][0], &K_abs_C1[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	K_abs_C2x_backward_2D = fftw_plan_dft_2d(Nx, Ny, &K_abs_C2[0][0], &K_abs_C2[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

	dEta_dx_backward_2D = fftw_plan_dft_2d(Nx, Ny, &dEta_dx[0][0], &dEta_dx[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
	dEta_dy_backward_2D = fftw_plan_dft_2d(Nx, Ny, &dEta_dy[0][0], &dEta_dy[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
}

void Plans_Destroying()
{

	fftw_destroy_plan(C1k_to_C1x_2D);
	fftw_destroy_plan(C2k_to_C2x_2D);

	fftw_destroy_plan(dC1_k_to_dC1_x_2D);
	fftw_destroy_plan(dC2_k_to_dC2_x_2D);

	fftw_destroy_plan(Eta_1_forward_2D);
	fftw_destroy_plan(Eta_1_backward_2D);

	fftw_destroy_plan(Eta_2_to_Eta_2_forward_2D);
	fftw_destroy_plan(Eta_2_to_Eta_2_backward_2D);

	fftw_destroy_plan(Mu_2_to_Mu_2_backward_2D);

	fftw_destroy_plan(NL_1p_to_NL_1p_forward_2D);
	fftw_destroy_plan(NL_1m_to_NL_1m_forward_2D);

	fftw_destroy_plan(NL_2p_to_NL_2p_forward_2D);
	fftw_destroy_plan(NL_2m_to_NL_2m_forward_2D);

	fftw_destroy_plan(NL_2p_to_NL_2p_backward_2D);
	fftw_destroy_plan(NL_2m_to_NL_2m_backward_2D);

	fftw_destroy_plan(NL_3p_to_NL_3p_forward_2D);
	fftw_destroy_plan(NL_3m_to_NL_3m_forward_2D);

	fftw_destroy_plan(NL_4p_to_NL_4p_forward_2D);
	fftw_destroy_plan(NL_4m_to_NL_4m_forward_2D);

	fftw_destroy_plan(NL_4p_to_NL_4p_backward_2D);
	fftw_destroy_plan(NL_4m_to_NL_4m_backward_2D);

	fftw_destroy_plan(K_abs_C1x_forward_2D);
	fftw_destroy_plan(K_abs_C2x_forward_2D);

	fftw_destroy_plan(K_abs_C1x_backward_2D);
	fftw_destroy_plan(K_abs_C2x_backward_2D);

	fftw_destroy_plan(dEta_dx_backward_2D);
	fftw_destroy_plan(dEta_dy_backward_2D);

}

double MaxFinder(fftw_complex Ck[][Ny])
{
	double Max = 0.0;

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			if (cabs(Ck[kx][ky]) > Max)
			{
				Max = cabs(Ck[kx][ky]);
			}
		}

	}

	return Max;
}

void NonLin_1()
{

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_1p[kx][ky] = I * (C1_x[kx][ky] * conj(C1_x[kx][ky]) - C2_x[kx][ky] * conj(C2_x[kx][ky])) * dC1_x[kx][ky];
			NL_1m[kx][ky] = I * (C2_x[kx][ky] * conj(C2_x[kx][ky]) - C1_x[kx][ky] * conj(C1_x[kx][ky])) * dC2_x[kx][ky];
		}
	}

	fftw_execute(NL_1p_to_NL_1p_forward_2D); // Финальное преобразование в тот же массив
	fftw_execute(NL_1m_to_NL_1m_forward_2D);


}

void NonLin_2()
{
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_2p[kx][ky] = (C1_x[kx][ky] * conj(C1_x[kx][ky]) - C2_x[kx][ky] * conj(C2_x[kx][ky]));
			NL_2m[kx][ky] = (C2_x[kx][ky] * conj(C2_x[kx][ky]) - C1_x[kx][ky] * conj(C1_x[kx][ky]));
		}
	}

	fftw_execute(NL_2p_to_NL_2p_forward_2D); // Фурье преобразование записывается в тот же самый массив!
	fftw_execute(NL_2m_to_NL_2m_forward_2D);// Фурье преобразование записывается в тот же самый массив!

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_2p[kx][ky] = K[kx][ky] * NL_2p[kx][ky]/(Nx * Ny); // Действуем оператором k, нормируем двумерное Фурье, записываем результат в тот же массив!
			NL_2m[kx][ky] = K[kx][ky] * NL_2m[kx][ky]/(Nx * Ny);
		}
	}

	fftw_execute(NL_2p_to_NL_2p_backward_2D); // Обратное Фурье преобразование в тот же самый массив!
	fftw_execute(NL_2m_to_NL_2m_backward_2D); // Обратное Фурье преобразование в тот же самый массив!

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_2p[kx][ky] = C1_x[kx][ky] * NL_2p[kx][ky]; // Получаем слагаемое С+ k (|C+|^2 - |C-|^2) и записываем его в тот же массив
			NL_2m[kx][ky] = C2_x[kx][ky] * NL_2m[kx][ky]; // Получаем слагаемое С- k (|C-|^2 - |C+|^2) и записываем его в тот же массив

		}
	}

	fftw_execute(NL_2p_to_NL_2p_forward_2D); // Финальное преобразование в тот же массив
	fftw_execute(NL_2m_to_NL_2m_forward_2D);
}

void NonLin_3()
{
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_3p[kx][ky] = I * C1_x[kx][ky] * C2_x[kx][ky] * conj(dC2_x[kx][ky]);
			NL_3m[kx][ky] = I * C2_x[kx][ky] * C1_x[kx][ky] * conj(dC1_x[kx][ky]);
		}
	}

	fftw_execute(NL_3p_to_NL_3p_forward_2D); // Фурье преобразование в тот же массив
	fftw_execute(NL_3m_to_NL_3m_forward_2D);
}

void NonLin_4()
{
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_4p[kx][ky] = C1_x[kx][ky] * C2_x[kx][ky];
			NL_4m[kx][ky] = C1_x[kx][ky] * C2_x[kx][ky];

			K_abs_C1[kx][ky] = C1_x[kx][ky] * conj(C1_x[kx][ky]);
			K_abs_C2[kx][ky] = C2_x[kx][ky] * conj(C2_x[kx][ky]);
		}
	}

	fftw_execute(NL_4p_to_NL_4p_forward_2D); // Фурье преобразование в тот же массив
	fftw_execute(NL_4m_to_NL_4m_forward_2D);

	fftw_execute(K_abs_C1x_forward_2D); // для расчёта энергии
	fftw_execute(K_abs_C2x_forward_2D);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_4p[kx][ky] = K[kx][ky] * NL_4p[kx][ky]/(Nx * Ny); // Действуем оператором k и записываем в тот же массив
			NL_4m[kx][ky] = K[kx][ky] * NL_4m[kx][ky]/(Nx * Ny);

			K_abs_C1[kx][ky] = K[kx][ky] * K_abs_C1[kx][ky]/(Nx * Ny); // слагаемые для энергии k (|C+|^2) и  k (|C-|^2) 
			K_abs_C2[kx][ky] = K[kx][ky] * K_abs_C2[kx][ky]/(Nx * Ny);
		}
	}

	fftw_execute(NL_4p_to_NL_4p_backward_2D); // Фурье преобразование в тот же массив
	fftw_execute(NL_4m_to_NL_4m_backward_2D);

	fftw_execute(K_abs_C1x_backward_2D); // слагаемые для энергии k (|C+|^2) и  k (|C-|^2) в х-пространстве!
	fftw_execute(K_abs_C2x_backward_2D);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			NL_4p[kx][ky] = conj(C2_x[kx][ky]) * NL_4p[kx][ky]; // Получаем слагаемое С*- k ( C+ C-) и записываем его в тот же массив
			NL_4m[kx][ky] = conj(C1_x[kx][ky]) * NL_4m[kx][ky]; // Получаем слагаемое С*+ k ( C+ C-) и записываем его в тот же массив

			E_int2[kx][ky] = conj(C1_x[kx][ky]) * NL_4p[kx][ky]; // c+* c-* k (c+ c-) для энергии
		}
	}

	fftw_execute(NL_4p_to_NL_4p_forward_2D); // Фурье преобразование в тот же массив
	fftw_execute(NL_4m_to_NL_4m_forward_2D);

}

void RHS(fftw_complex RHS[][Ny])
{
	#pragma omp parallel for num_threads(all_threads)
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			dC1_k[kx][ky] = I * K[kx][ky] * C1_k[kx][ky];
			dC2_k[kx][ky] = 0.0 + I * 0.0;
		}

		for (int kx = Nx/2; kx < Nx; kx++)
		{
			dC1_k[kx][ky] = 0.0 + I * 0.0;
			dC2_k[kx][ky] = -I * K[kx][ky] * C2_k[kx][ky];
		}
	}


	#pragma omp parallel private(t_id) num_threads(all_threads)
	{
		t_id = omp_get_thread_num(); // номер потока

		if (t_id == 0)
		{
			fftw_execute(dC1_k_to_dC1_x_2D);
		}

		if (t_id == 1)
		{
			fftw_execute(dC2_k_to_dC2_x_2D);
		}

		if (t_id == 2)
		{
			fftw_execute(C1k_to_C1x_2D);
		}

		if (t_id == 3)
		{
			fftw_execute(C2k_to_C2x_2D);
		}
	}

	#pragma omp parallel private(t_id) num_threads(all_threads)
	{
		t_id = omp_get_thread_num(); // номер потока

		if(t_id == 0)
		{
			NonLin_1();	
		}
		if(t_id == 1)
		{
			NonLin_2();	
		}
		if(t_id == 2)
		{
			NonLin_3();	
		}
		if(t_id == 3)
		{
			NonLin_4();	
		}
	}					

	#pragma omp parallel for num_threads(all_threads)
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			RHS[kx][ky] = -I * W[kx][ky] * C1_k[kx][ky] +  I * K[kx][ky]/(Nx * Ny) * (NL_1p[kx][ky] + NL_2p[kx][ky] - NL_3p[kx][ky] - NL_4p[kx][ky]);
		}

		for (int kx = Nx/2; kx < Nx; kx++)
		{
			RHS[kx][ky] = -I * W[kx][ky] * C2_k[kx][ky] - I * K[kx][ky]/(Nx * Ny) * (NL_1m[kx][ky] - NL_2m[kx][ky] - NL_3m[kx][ky] + NL_4m[kx][ky]);
		}

	}	

}

void Shift(fftw_complex RHS[][Ny], double step)
{
	#pragma omp parallel for num_threads(all_threads)
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			C1_k[kx][ky] = C1_k_old[kx][ky] + step * RHS[kx][ky];
			C2_k[kx][ky] = 0.0 + I * 0.0;
		}

		for (int kx = Nx/2; kx < Nx; kx++)
		{
			C1_k[kx][ky] = 0.0 + I * 0.0;
			C2_k[kx][ky] = C2_k_old[kx][ky] + step * RHS[kx][ky];
		}
		
	}
}

int RK4_step()
{
	#pragma omp parallel for num_threads(all_threads)
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			C1_k_new[kx][ky] = C1_k_old[kx][ky] + (tau/6.0) * (RHS_1[kx][ky] + 2.0 * RHS_2[kx][ky] + 2.0 * RHS_3[kx][ky] + RHS_4[kx][ky]);
			C2_k_new[kx][ky] = 0.0 + I * 0.0;
		}

		for (int kx = Nx/2; kx < Nx; kx++)
		{
			C2_k_new[kx][ky] = C2_k_old[kx][ky] + (tau/6.0) * (RHS_1[kx][ky] + 2.0 * RHS_2[kx][ky] + 2.0 * RHS_3[kx][ky] + RHS_4[kx][ky]);
			C1_k_new[kx][ky] = 0.0 + I * 0.0;
		}
	}

	#pragma omp parallel for num_threads(all_threads)
	for (int ky = Ny/4; ky < 3 * Ny/4; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			C1_k_new[kx][ky] = 0.0 + I * 0.0;
			C2_k_new[kx][ky] = 0.0 + I * 0.0;
		}
	}

	#pragma omp parallel for num_threads(all_threads)
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			C1_k_old[kx][ky] = C1_k_new[kx][ky];
			C2_k_old[kx][ky] = C2_k_new[kx][ky];
			C1_k[kx][ky] = C1_k_new[kx][ky];
			C2_k[kx][ky] = C2_k_new[kx][ky];
		}
	}

}

double aver_ampl_x()
{
	double result = 0.0;

	for (int ky = 0; ky < ky_damping; ky++)
	{
		for (int kx = Nx/2 - kx_high; kx < Nx/2; kx++)
		{
			result = result + cabs(C1_k[kx][ky]);
			result = result + cabs(C1_k[kx][Ny - ky]);
		}

		for (int kx = Nx/2 + 1; kx < Nx/2 + kx_high; kx++)
		{
			result = result + cabs(C2_k[kx][ky]);
			result = result + cabs(C2_k[kx][Ny - ky]);
		}
	}

	result = result/((2.0 * kx_high - 1) * (2.0 * ky_damping + 1));

	return result;
}

double aver_ampl_y()
{
	double result = 0.0;

	for (int kx = 0; kx < kx_damping; kx++)
	{
		for (int ky = Ny/4 - ky_high; ky < Ny/4; ky++)
		{
			result = result + cabs(C1_k[kx][ky]);
			result = result + cabs(C2_k[Nx - kx][ky]);
		}

		for (int ky = Ny - Ny/4; ky < Ny - Ny/4 + ky_high; ky++)
		{
			result = result + cabs(C1_k[kx][ky]);
			result = result + cabs(C2_k[Nx - kx][ky]);
		}
	}

	result = result/(2.0 * ky_high * (2.0 * kx_damping + 1));

	return result;
}

void Damping_procedure()
{

	// damping in kx
/*	if (aver_ampl_x() > Threshold)
	{*/
		#pragma omp parallel for num_threads(all_threads)
		for(int i = 0; i < Ny; i++)
		{
			for(int j = kx_damping; j < Nx/2; j++)
			{
				C1_k[j][i] = exp( -Dx * sqrt( log( cosh( pow( Alpha_x * (j - kx_damping) * Dk_x, 2.0))))) * C1_k[j][i];
			}

			for(int j = Nx/2 + 1; j < Nx - kx_damping; j++)
			{
				C2_k[j][i] = exp( -Dx * sqrt( log( cosh( pow( Alpha_x * (j - Nx + kx_damping) * Dk_x, 2.0 ))))) * C2_k[j][i];
			}
		}
	//}


	// damping in ky
/*	if (aver_ampl_y() > Threshold)
	{*/
		#pragma omp parallel for num_threads(all_threads)
		for(int j = 0; j < Nx/2; j++)
		{
			for(int i = ky_damping; i < Ny/4; i++)
			{
				C1_k[j][i] = exp( -Dy * sqrt( log( cosh( pow( Alpha_y * (i - ky_damping) * Dk_y, 2.0 ))))) * C1_k[j][i];	
			}

			for(int i = Ny - Ny/4; i <= Ny - ky_damping; i++)
			{
				C1_k[j][i] = exp( -Dy * sqrt( log( cosh( pow( Alpha_y * (i - Ny + ky_damping) * Dk_y, 2.0 )))) ) * C1_k[j][i];	
			}
		}

		#pragma omp parallel for num_threads(all_threads)
		for(int j = Nx/2; j < Nx; j++)
		{
			for(int i = ky_damping; i < Ny/4; i++)
			{
				C2_k[j][i] = exp( -Dy * sqrt( log( cosh( pow( Alpha_y * (i - ky_damping) * Dk_y, 2.0 ))))) * C2_k[j][i];	
			}

			for(int i = Ny - Ny/4; i <= Ny - ky_damping; i++)
			{
				C2_k[j][i] = exp( -Dy * sqrt( log( cosh( pow( Alpha_y * (i - Ny + ky_damping) * Dk_y, 2.0 )))) ) * C2_k[j][i];	
			}
		}
	//}
}

void Eta_Mu_calculating()
{
	int kx, ky;

	for (ky = 0; ky < Ny; ky++)
	{
		for (kx = 0; kx < Nx; kx++)
		{
			Eta_1[kx][ky] = (C1_x[kx][ky] + C2_x[kx][ky]) + conj(C1_x[kx][ky] + C2_x[kx][ky]);
		}
	}

	fftw_execute(Eta_1_forward_2D);

	for (ky = 0; ky < Ny/2; ky++)
	{
		for (kx = 1; kx < Nx/2; kx++)
		{
			Eta_1[kx][ky] = 1.0/(sqrt(2.0) * pow(g, 0.25)) * pow(K[kx][ky], -0.25) * Eta_1[kx][ky]/(Nx * Ny);
		}

		for (kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Eta_1[kx][ky] = 1.0/(sqrt(2.0) * pow(g, 0.25)) * pow(K[kx][ky], -0.25) * Eta_1[kx][ky]/(Nx * Ny);
		}
	}

	for (ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (kx = 1; kx < Nx/2; kx++)
		{
			Eta_1[kx][ky] = 1.0/(sqrt(2.0) * pow(g, 0.25)) * pow(K[kx][ky], -0.25) * Eta_1[kx][ky]/(Nx * Ny);
		}

		for (kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Eta_1[kx][ky] = 1.0/(sqrt(2.0) * pow(g, 0.25)) * pow(K[kx][ky], -0.25) * Eta_1[kx][ky]/(Nx * Ny);
		}
	}

	for (kx = 0; kx < Nx; kx++)
	{
		Eta_1[kx][Ny/2] = 0.0 + I * 0.0;
	}

	for (ky = 0; ky < Ny; ky++)
	{
		Eta_1[0][ky] = 0.0 + I * 0.0;
		Eta_1[Nx/2][ky] = 0.0 + I * 0.0;
	}

	for (ky = 0; ky < Ny; ky++)
	{
		for (kx = 0; kx < Nx/2; kx++)
		{
			dEta_dx[kx][ky] = I * kx * Dk_x * Eta_1[kx][ky];
		}

		for (kx = Nx/2; kx < Nx; kx++)
		{
			dEta_dx[kx][ky] = I * (kx - Nx) * Dk_x * Eta_1[kx][ky];
		}
	}

	for (kx = 0; kx < Nx; kx++)
	{
		for (ky = 0; ky < Ny/2; ky++)
		{
			dEta_dy[kx][ky] = I * ky * Dk_y * Eta_1[kx][ky];
		}

		for (ky = Ny/2; ky < Ny; ky++)
		{
			dEta_dy[kx][ky] = I * (ky - Ny) * Dk_y * Eta_1[kx][ky];
		}
	}

	fftw_execute(dEta_dx_backward_2D);
	fftw_execute(dEta_dy_backward_2D);

	fftw_execute(Eta_1_backward_2D);

	for (ky = 0; ky < Ny; ky++)
	{
		for(kx = 0; kx < Nx; kx++)
		{
			Mu_1[kx][ky] = sqrt(dEta_dx[kx][ky] * dEta_dx[kx][ky] + dEta_dy[kx][ky] * dEta_dy[kx][ky]);
		}
	}
}

void Eta_Mu_2_calculating()
{
	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			Eta_2[kx][ky] = (C1_x[kx][ky] - conj(C1_x[kx][ky]) - C2_x[kx][ky] + conj(C2_x[kx][ky]));
		}
	}

	fftw_execute(Eta_2_to_Eta_2_forward_2D);

	for (int ky = 0; ky < Ny/2; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Eta_2[kx][ky] = pow(K[kx][ky], -0.25) * Eta_2[kx][ky]/(Nx * Ny);
		}

		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Eta_2[kx][ky] = pow(K[kx][ky], -0.25) * Eta_2[kx][ky]/(Nx * Ny);
		}
	}

	for (int ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Eta_2[kx][ky] = pow(K[kx][ky], -0.25) * Eta_2[kx][ky]/(Nx * Ny);
		}

		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Eta_2[kx][ky] = pow(K[kx][ky], -0.25) * Eta_2[kx][ky]/(Nx * Ny);
		}
	}

	for (int kx = 0; kx < Nx; kx++)
	{
		Eta_2[kx][Ny/2] = 0.0 + I * 0.0;
	}

	for (int ky = 0; ky < Ny; ky++)
	{
		Eta_2[0][ky] = 0.0 + I * 0.0;
		Eta_2[Nx/2][ky] = 0.0 + I * 0.0;
	}

	fftw_execute(Eta_2_to_Eta_2_backward_2D);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			Eta_2[kx][ky] = Eta_2[kx][ky] * Eta_2[kx][ky];
		}
	}

	fftw_execute(Eta_2_to_Eta_2_forward_2D);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx/2; kx++)
		{
			Eta_2[kx][ky] = K[kx][ky]/(4.0 * sqrt(g)) * Eta_2[kx][ky]/(Nx * Ny);
			Mu_2[kx][ky] = I * kx * Dk_x * Eta_2[kx][ky];
		}

		for (int kx = Nx/2; kx < Nx; kx++)
		{
			Eta_2[kx][ky] = K[kx][ky]/(4.0 * sqrt(g)) * Eta_2[kx][ky]/(Nx * Ny);
			Mu_2[kx][ky] = I * (kx - Nx) * Dk_x * Eta_2[kx][ky];
		}
	}

	for (int kx = 0; kx < Nx; kx++)
	{
		Mu_2[kx][Ny/2] = 0.0 + I * 0.0;
	}

	for (int ky = 0; ky < Ny; ky++)
	{
		Mu_2[0][ky] = 0.0 + I * 0.0;
		Mu_2[Nx/2][ky] = 0.0 + I * 0.0;
	}

	fftw_execute(Eta_2_to_Eta_2_backward_2D);
	fftw_execute(Mu_2_to_Mu_2_backward_2D);
}

void Print_C_kxky(int index)
{
	int i, j;
	double kx, ky;

	sprintf(Name, "./Ck//C_kxky_%d.dat", index);

	pfile = fopen(Name,"w");

	for (i = Ny - Ny/4 + 1; i < Ny; i++)
	{
		ky = Dk_y * ( (double)(i - Ny) );

		for (j = Nx - Nx/2 + 1; j < Nx; j++)
		{
			kx = Dk_x * ( (double)(j - Nx) );
			fprintf (pfile, "%e   %e   %e\n", ky, kx, cabs(C2_k[j][i]));
		}

		for(j = 0; j < Nx/2; j++)
		{
			kx = Dk_x * ( (double)j );
			fprintf (pfile, "%e   %e   %e\n", ky, kx, cabs(C1_k[j][i]));
		}

		fprintf (pfile, "\n");
	}

	for(i = 0; i < Ny/4; i++)
	{
		ky = Dk_y * ( (double)i );
		for(j = Nx - Nx/2 + 1; j < Nx; j++)
		{
			kx = Dk_x * ( (double)(j - Nx) );
			fprintf (pfile, "%e   %e   %e\n", ky, kx, cabs(C2_k[j][i]));
		}

		for(j = 0; j < Nx/2; j++)
		{
			kx = Dk_x * ( (double)j );
			fprintf (pfile, "%e   %e   %e\n", ky, kx, cabs(C1_k[j][i]));
		}
		fprintf (pfile, "\n");
	}

	fclose(pfile);

}

void Print_b_1D(int index)
{

	for (int kx = 1; kx < Nx/2; kx++)
	{
		Bk_x[kx] = C1_k[kx][0]/sqrt(K[kx][0]);
	}

	for (int kx = Nx/2 + 1; kx < Nx; kx++)
	{
		Bk_x[kx] = C2_k[kx][0]/sqrt(K[kx][0]);
	}

	Bk_x[0] = 0.0 + I * 0.0;
	Bk_x[Nx/2] = 0.0 + I * 0.0;

	sprintf(Name, "./b_k//bk_1d_%d.dat", index);

	pfile = fopen(Name,"w");

	for(int j = Nx/2; j < Nx; j++)
	{
		fprintf (pfile, "%e\t%e\n", (j - Nx) * Dk_x, cabs(Bk_x[j]));
	}

	for(int j = 0; j < Nx/2; j++)
	{
		fprintf (pfile, "%e\t%e\n", j * Dk_x, cabs(Bk_x[j]));
	}

	fclose(pfile);
}


void Print_Eta_Mu(int index)
{
	Eta_Mu_calculating();
	//Eta_Mu_2_calculating();

	sprintf(Name, "./Eta_Mu//Eta_Mu_%d.dat", index);
	pfile = fopen(Name, "w"); 
	CheckFile(pfile);

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			fprintf(pfile, "%e\t%e\t%e\t%e\n", ax + kx * hx, ay + ky * hy, creal(Eta_1[kx][ky]), creal(Mu_1[kx][ky]));
		}
		fprintf (pfile, "\n");
	}

	fclose(pfile);


}

double N_plus()
{
	double Num = 0.0;

	for (int ky = 0; ky < Ny/2; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Num = Num + creal((C1_k[kx][ky] * conj(C1_k[kx][ky]))/K[kx][ky]);
		}
	}

	for (int ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Num = Num + creal((C1_k[kx][ky] * conj(C1_k[kx][ky]))/K[kx][ky]);
		}
	}


	return Num;
}

double N_minus()
{
	double Num = 0.0;

	for (int ky = 0; ky < Ny/2; ky++)
	{
		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Num = Num + creal((C2_k[kx][ky] * conj(C2_k[kx][ky]))/K[kx][ky]);
		}
	}

	for (int ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Num = Num + creal((C2_k[kx][ky] * conj(C2_k[kx][ky]))/K[kx][ky]);
		}
	}

	return Num;
}

double P_x()
{
	double Mom = 0.0;

	for (int ky = 0; ky < Ny/2; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Mom = Mom + (kx * Dk_x)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}

		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Mom = Mom + ((kx - Nx) * Dk_x)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}
	}

	for (int ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Mom = Mom + (kx * Dk_x)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}

		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Mom = Mom + ((kx - Nx) * Dk_x)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}
	}

	return Mom;
}


double P_y()
{
	double Mom = 0.0;

	for (int ky = 0; ky < Ny/2; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Mom = Mom + (ky * Dk_y)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}

		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Mom = Mom + (ky * Dk_y)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}
	}

	for (int ky = Ny/2 + 1; ky < Ny; ky++)
	{
		for (int kx = 1; kx < Nx/2; kx++)
		{
			Mom = Mom + ((ky - Ny) * Dk_y)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}
		for (int kx = Nx/2 + 1; kx < Nx; kx++)
		{
			Mom = Mom + ((ky - Ny) * Dk_y)/K[kx][ky] * (C1_k[kx][ky] * conj(C1_k[kx][ky]) + C2_k[kx][ky] * conj(C2_k[kx][ky]));
		}
	}

	return Mom;
}

double E_lin()
{
	double E_lin = 0.0;

	for (int kx = 1; kx < Nx/2; kx++)
	{
		for (int ky = 0; ky < Ny/2; ky++)
		{
			E_lin = E_lin + creal(conj(C1_k[kx][ky]) * W[kx][ky]/K[kx][ky] * C1_k[kx][ky]);
		}

		for (int ky = Ny/2 + 1; ky < Ny; ky++)
		{
			E_lin = E_lin + creal(conj(C1_k[kx][ky]) * W[kx][ky]/K[kx][ky] * C1_k[kx][ky]);
		}
	}

	for (int kx = Nx/2 + 1; kx < Nx; kx++)
	{
		for (int ky = 0; ky < Ny/2; ky++)
		{
			E_lin = E_lin + creal(conj(C2_k[kx][ky]) * W[kx][ky]/K[kx][ky] * C2_k[kx][ky]);
		}

		for (int ky = Ny/2 + 1; ky < Ny; ky++)
		{
			E_lin = E_lin + creal(conj(C2_k[kx][ky]) * W[kx][ky]/K[kx][ky] * C2_k[kx][ky]);
		}
	}

	return E_lin/g;
}


double E_Nonlin()
{
	double E_nonlin = 0.0;

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			E_nonlin = E_nonlin + creal(C1_x[kx][ky] * conj(C1_x[kx][ky])/2.0 * (cimag(conj(C1_x[kx][ky]) * dC1_x[kx][ky]) - K_abs_C1[kx][ky]));	
		}
	}

	for (int ky = 0; ky < Ny; ky++)
	{
		for (int kx = 0; kx < Nx; kx++)
		{
			E_nonlin = E_nonlin + creal(C2_x[kx][ky] * conj(C2_x[kx][ky])/2.0 * (cimag(C2_x[kx][ky] * conj(dC2_x[kx][ky])) - K_abs_C2[kx][ky]));
		}
	}

	return E_nonlin/(Nx * Ny * g);
}

double E_interact()
{
	double E_int = 0.0;

	for (int kx = 0; kx < Nx; kx++)
	{
		for (int ky = 0; ky < Ny; ky++)
		{
			E_int = E_int + creal(C1_x[kx][ky] * conj(C1_x[kx][ky]) * K_abs_C2[kx][ky]);
			E_int = E_int + creal(E_int2[kx][ky]);
			E_int = E_int + creal(I * conj(C1_x[kx][ky]) * C2_x[kx][ky] * (dC1_x[kx][ky] * conj(C2_x[kx][ky]) + C1_x[kx][ky] * conj(dC2_x[kx][ky])));
		}
	}

	return E_int/(Nx * Ny * g);
}


void Print_NPE()
{
	pfile = fopen("NPE.txt", "a");
	CheckFile(pfile);

	if (nstep == 0)
	{
		fprintf (pfile, "%s\t\t%s\t\t\t%s\t\t\t%s\t%s\t%s\n", " Time [t] ", " N_plus [N+] ", " N_minus [N-]", " Total Momentum [Px] ", " Total Momentum [Py] ", " Total Energy [E] ");
	}
	
	fprintf (pfile, "%e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", t,  N_plus(), N_minus(), P_x(), P_y(), E_lin() + E_Nonlin() + E_interact());

	fclose(pfile);
}

void Increment()
{
	pfile = fopen("Gamma.txt", "a");
	CheckFile(pfile);

	fprintf (pfile, "%e\t%e\t%e\n", t, cabs(C1_k[130][3]),cabs(C1_k[70][Ny - 3]));

	fclose(pfile);
}


void Cleaning()
{
	remove("NPE.txt");

	system("exec rm -r ./Ck/*");
	system("exec rm -r ./Eta_Mu/*");
	system("exec rm -r ./CheckPoints/*");
}

void CheckPoint(int index)
{
	int binary_number;

	sprintf(Name, "./CheckPoints//C1_k_data_%d", index);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);

	binary_number = fwrite(&C1_k[0][0], sizeof(fftw_complex), Nx * Ny, pfile);

	fclose(pfile);

	//sprintf(Name, "./CheckPoints//C2_k_data_%d", index);
	sprintf(Name, "./CheckPoints//C2_k_data_%d", index);
	pfile = fopen (Name,"w+b");
	CheckFile(pfile);

	binary_number = fwrite(&C2_k[0][0], sizeof(fftw_complex), Nx * Ny, pfile);
	fclose(pfile);
}

void Evolution(double Time, int N_rec, int N_En)
{

	printf("%s\n", "Calculation has been started");
	
	while (t < Time + tau)
	{	

		RHS(RHS_1);

		Shift(RHS_1, 0.5 * tau);

		RHS(RHS_2);

		Shift(RHS_2, 0.5 * tau);

		RHS(RHS_3);

		Shift(RHS_3, tau);

		RHS(RHS_4);

		RK4_step();

		Damping_procedure();

		if (!(nstep%(N_rec)))
		{
			Print_C_kxky(nstep/N_rec);
			Print_b_1D(nstep/N_rec);
			//Increment();

			
			Print_Eta_Mu(nstep/N_rec);

			printf("t = %e\t steps = %d\t mu_max = %e\n", t, nstep, MaxFinder(Mu_1));
		} 

		if (!(nstep%(N_En)))
		{
			Print_NPE();
			CheckPoint(nstep/N_rec);
		} 


		t = t + tau;
		nstep++;
		
	}
}
void Step_Check(double tau)
{
		double taumax_Nq_x = M_PI/(4.0 * sqrt(g * Dk_x * Nx/2.0)); // Критерий по частоте Найквиста
		double taumax_Nq_y = M_PI/(4.0 * sqrt(g * Dk_y * Ny/2.0)); // Критерий по частоте Найквиста
		double taumax_Crnt_x = 2.0 * (bx - ax) * sqrt(Dk_x)/(Nx * sqrt(g)); // Критерий Куранта
		double taumax_Crnt_y = 2.0 * (by - ay) * sqrt(Dk_y)/(Ny * sqrt(g)); // Критерий Куранта

		double Least_Nq = 0.0;
		double Least_Crnt = 0.0;

		if (taumax_Nq_x < taumax_Nq_y)
		{
			Least_Nq = taumax_Nq_x;
		}
		else
		{
			Least_Nq = taumax_Nq_y;
		}

		if (taumax_Crnt_x < taumax_Crnt_y)
		{
			Least_Crnt = taumax_Crnt_x;
		}
		else
		{
			Least_Crnt = taumax_Crnt_y;
		}


		if (tau > Least_Nq)
		{
			printf("Warning! The time step is too high. Please decrease the step to avoid possible calculation instability. \n Tau = %e\t  Tau_Nyquist_least = %e\t  Tau_Courant_least = %e\n", tau, Least_Nq, Least_Crnt);
			exit(0);
		}
		else
		{
			printf("Time step is adequate. \n Tau = %e\t Tau_Nq = %e\t Tau_Crnt = %e\n\n", tau, Least_Nq, Least_Crnt);
		}
}

void Multi_threads_info()
{
	int nthreads, tid;

	pfile = fopen("MP_info.txt", "w");
	CheckFile(pfile);

	#pragma omp parallel private(nthreads, tid, Name)
    {
        tid = omp_get_thread_num();
        gethostname(Name, 50);
        fprintf(pfile, "Thread %d started on %s and says 'Hello!' \n ", tid, Name);

        if (tid == 0) 
        {
        	nthreads = omp_get_max_threads();
            fprintf(pfile, "[Thread 0 say: there are %d threads] \n", nthreads);
        }

    }

    fclose(pfile);

}




///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////// 				MAIN
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////









int main(void)
{
	s_time = time(NULL); // запуск таймера
	srand(time(NULL)); // Обнуление генератора функции задания случайных чисел
		
	Set_data(); // Задание всех необходимых массивов и параметров

	Cleaning(); // Очистка старых файлов

	all_threads = 32;
	fftw_threads = 8;

	omp_set_num_threads(all_threads);
	fftw_init_threads();
	fftw_plan_with_nthreads(fftw_threads);

	//int fftw_planner_nthreads(void);

	Plans_Creating(); // Создание планов Фурье

	Multi_threads_info(); // Информация о количестве потоков

	Print_C_kxky(0);
	Print_b_1D(0);

	fftw_execute(C1k_to_C1x_2D);
	fftw_execute(C2k_to_C2x_2D);

	Print_Eta_Mu(0);

	t = 0.0; // начальное время
	nstep = 0; // число шагов
	tau = 0.01; // шаг по времени
	Step_Check(tau);

	double fin_time = 3000.0;

	int N_rec = 1000; // Запись эволюции
	int N_En = 1000; // Запись энергии
	printf("%s\n", "Calculations will start with next initial conditions:");
	printf(" Lx= %e\t Ly = %e\t Eta_max_ini = %e\t Mu_max_ini = %e\t finTime = %e\n\n", bx - ax, by - ay, MaxFinder(Eta_1), MaxFinder(Mu_1), fin_time);

	Evolution(fin_time, N_rec, N_En);

	Plans_Destroying(); // Освобождение памяти
	fftw_cleanup_threads();

	f_time = time(NULL); // остановка таймера

	printf("It required %ld seconds for the program to complete \n", f_time - s_time);


	return 0;
}
