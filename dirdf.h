/* Dynamic Precision Linear Solver using XBLAS
JunKyu Lee
date : july. 2017	
This iterative refinement produces original (double) precision accuracy in forward error 
using adaptive precision 
*/

#include <stdio.h> 
#include <typeinfo>
#include "blas_extended.h" //xblas library
#include <mkl.h> 
#define eps_s 0.0000000596 //single prec machine epsilon
#define eps 0.000000000000000111 //double prec machine epsilon
#define standpoint 0.000000001863 // minimal required accuracy when double precision step 1 is saturated. 
#define L_THRESHOLD 10 //latency gap threshold between dbl and dbl-dbl for activation of inner loop

using namespace std;

double sdnormi (float* sp, int N) {	
	int i_indx;
	i_indx = cblas_isamax(N, sp, 1);
	return (double) fabs(sp[i_indx]);}

double ddnormi (double * dp, int N) {	
	int i_indx;
	i_indx = cblas_idamax(N, dp, 1);
	return fabs(dp[i_indx]);}
	
template <class OP, class LP>
int refine_s_ad(OP A, OP b, OP x, LP LU, OP dp, int *ipiv, int N, int ir_info[2]) {	
	int nrhs; nrhs = 1; //the number of right handside coefficient for LAPACK 
	char trans; trans = 'n'; // transpose info for LAPACK
	int i, j;
	double z_norm, del_z_norm, z_conv, cur_x_norm, x_conv; 
	static double pre_z_norm;
	static int ddon=0; //when convergence saturate, ddon = 1
	static int tt2_done=0; //when tt2 inner loop done, tt2_done=1;
	static int init=1;
	int info[1];
	
	double* resid_temp = (double *) malloc(N*sizeof(double));
	double* resid = (double *) malloc(N*sizeof(double));
	float* sp = (float *) malloc(N*sizeof(float));
	
	
	cur_x_norm = ddnormi(x, N);
	if(init) {pre_z_norm=cur_x_norm; init=0;}
	cblas_dcopy(N, b, 1, resid, 1); //copy b to residual
	
	//step 1: residual is sought using double precision for early iteration ans switch it to dbl-dbl 
	if (!ddon) cblas_dgemv(CblasColMajor,CblasNoTrans, N, N, -1, A, N, x,1,1, resid, 1);
	else BLAS_dgemv_x(blas_colmajor,blas_no_trans, N, N, -1, A, N, x,1,1, resid, 1, blas_prec_extra);
	//resid has the current residual 	

	for(i=0;i<N;i++) sp[i] = (float) resid[i]; //typecasting for residual for step 2
	
	sgetrs_(&trans, &N, &nrhs, LU, &N, ipiv, sp, &N, info); 
	//step2: seek the error in the prev sol. sp contains the approx of prev error
	for(i=0;i<N;i++) dp[i] = (double)sp[i]; //typecasting back to double for step 3 correction
	z_norm = ddnormi(dp, N); //infinity norm of the error	

	x_conv = z_norm/cur_x_norm; 
	if(x_conv < eps) {if(ir_info[0]) printf("Accuracy: %g\n",x_conv); return 0;} 
	
	//transprecision technique 2
	if (ddon&&(ir_info[1] > L_THRESHOLD)) 
		{  //printf("\ninner solver from ad...\n");
			for (i = 0; i < 30; i++)
				{	
					cblas_dcopy(N, resid, 1, resid_temp, 1);
					cblas_dgemv(CblasColMajor,CblasNoTrans, N, N, -1, A, N, dp,1, 1, resid_temp, 1);
					for(j=0;j<N;j++) {sp[j] = (float)resid_temp[j];}; 
					sgetrs_(&trans, &N, &nrhs, LU, &N, ipiv, sp, &N, info);
					del_z_norm = sdnormi(sp, N);
					BLAS_daxpby_s(N, 1, sp, 1, 1, dp, 1);
					if ((del_z_norm/z_norm) < eps_s) {tt2_done=1;  break;}
					z_norm = ddnormi(dp, N);
				}
			if (i==30-1) printf("TT2 Failure\n"); 	
		}
	
	z_conv = z_norm/pre_z_norm;//{printf("zconv:%g\n",z_conv);}
	
	if (z_conv > 0.5) //no convergence 
	{
		if (x_conv < standpoint) ddon = 1; //apply higher prec for step1 if no convergence due to lower prec for step1
		else {printf("\nFail by Slow Convergence\n"); return -1;}
	}
	
	BLAS_daxpby_x(N, 1, dp, 1, 1, x, 1, blas_prec_extra); //step3: sol update using dbl-dbl arith but using dbl data
	pre_z_norm = z_norm;
	if(tt2_done && (!ir_info[0])) { ddon = 0; tt2_done = 0; return 0;} //TT3 : no need to check accuracy 
		
	return 1; //it should be more iterated. 
};


int refine_d_dd_t(double *A, double *b, double *x, double *LU, double *p, int *ipiv, int N, int ir_info[2]) {	
		int nrhs, info[1]; nrhs = 1;
		char trans; trans = 'n';
		static int init = 1;
		double cur_z_norm, cur_x_norm, x_conv, z_conv;
		static double pre_z_norm; 
		
		cur_x_norm = ddnormi(x, N);
		if(init) {pre_z_norm = cur_x_norm; init=0;}
		
		cblas_dcopy(N, b, 1, p, 1);
		BLAS_dgemv_x(blas_colmajor,blas_no_trans, N, N, -1, A, N, x,1,1, p, 1, blas_prec_extra);
		dgetrs_(&trans, &N, &nrhs, LU, &N, ipiv, p, &N, info);
		cur_z_norm = ddnormi(p, N);
		x_conv = cur_z_norm/cur_x_norm; 
		if (x_conv < eps) {if(ir_info[0]) printf("Accuracy: %g\n", x_conv); return 0;} //successful accuracy
		
		z_conv = cur_z_norm/pre_z_norm;
		if (z_conv > 0.5) {printf("\nFail by Slow Convergence\n"); return -1;}
		BLAS_daxpby_x(N, 1, p, 1, 1, x, 1, blas_prec_extra);
		pre_z_norm = cur_z_norm;
		return 1;	//more iterations required
};

int refine_s_dd_t(double *A, double *b, double *x, float *LU, double *dp, int *ipiv, int N, int ir_info[2])
{	
	int i, nrhs; nrhs = 1;
	char trans; trans = 'n'; // Transpose info for LAPACK use
	int info[1];
	static int init = 1;
	
	double cur_z_norm, cur_x_norm, x_conv, z_conv;
	static double pre_z_norm;
	
	float* sp = (float*)malloc(N*sizeof(float));
	
	cur_x_norm = ddnormi(x, N);
	if(init) {pre_z_norm = cur_x_norm; init=0;}
	
	cblas_dcopy(N, b, 1, dp, 1);
	BLAS_dgemv_x(blas_colmajor,blas_no_trans, N, N, -1, A, N, x,1,1, dp, 1, blas_prec_extra);

	for(i=0;i<N;i++) sp[i] = (float)dp[i];
	
	sgetrs_(&trans, &N, &nrhs, LU, &N, ipiv, sp, &N, info);

	for(i=0;i<N;i++) dp[i] = (double)sp[i];
	cur_z_norm = ddnormi(dp, N);
	x_conv = cur_z_norm/cur_x_norm; 
	
	if (x_conv < eps) {if(ir_info[0]) printf("Accuracy: %g\n", x_conv); return 0;} //successful accuracy
		
	z_conv = cur_z_norm/pre_z_norm;
	if (z_conv > 0.5) {printf("\nFail by Slow Convergence\n"); return -1;}
	BLAS_daxpby_x(N, 1, dp, 1, 1, x, 1, blas_prec_extra);
	pre_z_norm = cur_z_norm;
	return 1;	
};	
 
template <class OP, class LP> 
int dirdf(OP A, OP b, OP x, LP sLU, OP dLU, int *ipiv, int N, int info[3]) {
//info[0]: irtype, info[1]: accuracy_chk on or off, info[2]: latency_dd/lantency_d for current architecture for transIR	
	double* p = (double*)malloc(N*sizeof(double));
	int tmp_info[2];
	int flag;
	tmp_info[0]=info[1]; tmp_info[1]=info[2]; flag = 1;
	
	switch(info[0]) {
		case 0:
		flag = refine_s_ad(A, b, x, sLU, p, ipiv, N, tmp_info);
		break;

		case 2:
		flag = refine_d_dd_t(A, b, x, dLU, p, ipiv, N, tmp_info);
		break;
		
		case 1:
		flag = refine_s_dd_t(A, b, x, sLU, p, ipiv, N, tmp_info);;
		break;
		
		default : printf("\nWrong IR type.\n"); flag = 0;
		};
	
	return flag;
}
