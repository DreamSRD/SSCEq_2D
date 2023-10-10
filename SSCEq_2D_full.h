
	#define M_PI 3.1415926535

	#define Nx 256
	#define Ny 256

	double ax = 0.0;
	double ay = 0.0;
	double bx = 300;
	double by = 1000;

	int k0 = 0;
	int w0 = 0;
	double V = 6.24920;

	double V_br = 0.0; // Скорость бризера

	double g = 9.81;
	double Dk_x, Dk_y,hx,hy,t,tau, T0, x1,x2, T;
	int nstep = 0;

	double Dx = 10.0;
	double Dy = 200.0;
	double Alpha_x = 1.0;
	double Alpha_y = 2.0;

	char Name[100];

	int kx_damping = 64;
	int ky_damping = 32;

	int kx_high = 24;
	int ky_high = 12;

	int t_id;
	int Max_int = 0;

	int all_threads;
	int fftw_threads;

	double Threshold = 1e-10;

	int binary_number;

	double K[Nx][Ny];
	double W[Nx][Ny];

	fftw_complex b_k_restart[256 * 128];
	fftw_complex Bk_x[Nx];

	fftw_complex C1_k[Nx][Ny];
	fftw_complex C2_k[Nx][Ny];

	fftw_complex C1_x[Nx][Ny];
	fftw_complex C2_x[Nx][Ny];

	fftw_complex C1_k_old[Nx][Ny];
	fftw_complex C2_k_old[Nx][Ny];
	fftw_complex C1_k_new[Nx][Ny];
	fftw_complex C2_k_new[Nx][Ny];

	fftw_complex dC1_x[Nx][Ny];
	fftw_complex dC2_x[Nx][Ny];

	fftw_complex dC1_k[Nx][Ny];
	fftw_complex dC2_k[Nx][Ny];

	fftw_complex RHS_1[Nx][Ny];
	fftw_complex RHS_2[Nx][Ny];
	fftw_complex RHS_3[Nx][Ny];
	fftw_complex RHS_4[Nx][Ny];


	fftw_complex NL_1p[Nx][Ny];
	fftw_complex NL_1m[Nx][Ny];

	fftw_complex NL_2p[Nx][Ny];
	fftw_complex NL_2m[Nx][Ny];

	fftw_complex NL_3p[Nx][Ny];
	fftw_complex NL_3m[Nx][Ny];

	fftw_complex NL_4p[Nx][Ny];
	fftw_complex NL_4m[Nx][Ny];

	fftw_complex Eta_1[Nx][Ny];
	fftw_complex Eta_2[Nx][Ny];

	fftw_complex dEta_dx[Nx][Ny];
	fftw_complex dEta_dy[Nx][Ny];

	fftw_complex Mu_1[Nx][Ny];
	fftw_complex Mu_2[Nx][Ny];

	fftw_complex K_abs_C1[Nx][Ny];
	fftw_complex K_abs_C2[Nx][Ny];

	fftw_complex E_int2[Nx][Ny];

	fftw_plan C1k_to_C1x_2D;
	fftw_plan C2k_to_C2x_2D;

	fftw_plan dC1_k_to_dC1_x_2D;
	fftw_plan dC2_k_to_dC2_x_2D;

	fftw_plan Eta_1_forward_2D;
	fftw_plan Eta_2_to_Eta_2_forward_2D;

	fftw_plan Eta_1_backward_2D;
	fftw_plan Eta_2_to_Eta_2_backward_2D;

	fftw_plan Mu_2_to_Mu_2_backward_2D;

	fftw_plan NL_1p_to_NL_1p_forward_2D;
	fftw_plan NL_1m_to_NL_1m_forward_2D;

	fftw_plan NL_2p_to_NL_2p_forward_2D;
	fftw_plan NL_2m_to_NL_2m_forward_2D;

	fftw_plan NL_2p_to_NL_2p_backward_2D;
	fftw_plan NL_2m_to_NL_2m_backward_2D;

	fftw_plan NL_3p_to_NL_3p_forward_2D;
	fftw_plan NL_3m_to_NL_3m_forward_2D;

	fftw_plan NL_4p_to_NL_4p_forward_2D;
	fftw_plan NL_4m_to_NL_4m_forward_2D;

	fftw_plan NL_4p_to_NL_4p_backward_2D;
	fftw_plan NL_4m_to_NL_4m_backward_2D;

	fftw_plan K_abs_C1x_forward_2D;
	fftw_plan K_abs_C1x_backward_2D;
	fftw_plan K_abs_C2x_forward_2D;
	fftw_plan K_abs_C2x_backward_2D;

	fftw_plan dEta_dx_backward_2D;
	fftw_plan dEta_dy_backward_2D;

	
	
	FILE *pfile;