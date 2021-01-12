/* Adaptive Precision Linear Solver using XBLAS
Author and Copyright: JunKyu Lee (Queens University Belfast) 

This iterative refinement produces original (double) precision accuracy in forward error 
using adaptive precision 
*/ 
 
#include <stdio.h>
#include <string.h>
//#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include "dirdf.h"
#define LDD_LD 20 //integer latency gap between dbl-dbl and dble
#define MAX_ITER 30

double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;
  
  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
  
  return cur_time;
};


int main(int argc, char* argv[]) {	
		
	int mode, iternum, nrhs, N, n2, i, flag, ir_info[3], info[1];
	double btime, etime;
	char trans, *ir_type;
	const char *transir, *mixir, *uniir;
	
	if(argc!=4) {fprintf(stderr,"e.g., ./dirdf trans(mix or uni) matrix_size(1024) accuracy_chk\n"); return 1; }
	N= atoi(argv[2]); n2 =N*N; 
	// Memory allocations for BLAS operations
	double* A = (double *)malloc(N*N*sizeof(double));	
	double* b = (double *)malloc(N*sizeof(double));	
	float* sLU = (float *)malloc(N*N*sizeof(float));	
	double* dLU = (double *)malloc(N*N*sizeof(double));	
	double* xhat = (double *)malloc(N*sizeof(double));
	float* sxhat = (float *)malloc(N*sizeof(float));
	int* ipiv = (int *)malloc(N*sizeof(int));
	
	transir = "trans";
	mixir = "mix";
	uniir = "uni";
	ir_type = argv[1];
	if(!strcmp(transir,ir_type)) {mode = 0;}
	else if(!strcmp(mixir,ir_type)) {mode = 1;} 
	else if(!strcmp(uniir,ir_type)) {mode = 2;} 
	else {printf("your refinement type is not defined."); return -1;}
	
	//Print accuracies. 
	//Notice that if 1, the runtime for one additional dbl-dbl residual calculation will be included for transIR.
	int ACCURACY_CHECK; 
	ACCURACY_CHECK = atoi(argv[3]);
	
	nrhs = 1;  trans = 'n'; 		
	ir_info[0] = mode; ir_info[1] = ACCURACY_CHECK; ir_info[2] = LDD_LD;
	
	printf ("\n%s-ir (%dx%d matrix): \n",ir_type,N,N);
	
	//a dense matrix generation  
	for(i=0; i<N; i++) b[i] = drand48();
	for(i= 0; i< n2; i++) A[i] = drand48(); 
		

	omp_set_num_threads(1);  //set the number of threads for OpenMP for BLAS

	btime = get_cur_time(); //start timing for runtime measure 
	
	if (mode == 2) { //for uni-ir approximator (LU factor) generation 		
		cblas_dcopy(n2, A, 1, dLU, 1); //copy A to dLU(LU factors) for uni-IR case
		dgetrf_(&N, &N, dLU, &N, ipiv, info); 
		cblas_dcopy(N, b, 1, xhat, 1);
		dgetrs_(&trans, &N, &nrhs, dLU, &N, ipiv, xhat, &N, info);
	}
	else { //for mix-ir or trans-ir approximator (LU factor) generation 
		for(i=0; i<n2; i++) sLU[i] =(float) A[i]; //copy A to sLU(LU factors) using typecasting
		sgetrf_(&N, &N, sLU, &N, ipiv, info);
		for(i=0;i<N;i++) sxhat[i] = (float)b[i];
		sgetrs_(&trans, &N, &nrhs, sLU, &N, ipiv, sxhat, &N, info);
		for(i=0;i<N;i++) xhat[i] = (double)sxhat[i];
	}

	flag = 1;
	
	for(iternum = 0; iternum < MAX_ITER; iternum++)
	{ 	
		flag = dirdf(A, b, xhat, sLU, dLU, ipiv, N, ir_info);
		if(!flag) {
			printf("Success\n"); 
			etime = get_cur_time();
			printf("Runtime: %f secs\n",etime-btime);
			return 0;}
		
		if (iternum == (MAX_ITER-1)) { 
			printf("\nExceed MAX number of iterations!!\n"); return -1;
		}
	}
	
	free(A); free(b); free(xhat);
	
	return 0;
};
		

	

