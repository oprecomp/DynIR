# Dynamic precision iterative refinement v1.0
Latest Update: 25 September 2020

C Codes for Iterative Refinement Using Dynamic Precision out of Single, Double, and Double-Double Precision.  

Reference: J. Lee, H. Vandierendonck, M. Arif, G. D. Peterson and D. S. Nikolopoulos, "Energy-Efficient Iterative Refinement Using Dynamic Precision," in IEEE Journal on Emerging and Selected Topics in Circuits and Systems, vol. 8, no. 4, pp. 722-735, Dec. 2018, doi: 10.1109/JETCAS.2018.2850665. 

# Description:
dirdf.h file (Dynamic Precision Iterative Refinement for Double Precision Accuracy Forward-Error) contains dirdf() function which takes either single or double LUPP approximator and produces double precision solution accuracy.

The command line is: 

./approximator refinement_type matrix_size accuracy_check (e.g., ./lupp trans(mix or uni) 1024 1)

*lupp: LU decomposition with Partial Pivoting
- uni-ir (double lupp approximator + double-double refinement)
- mix-ir (single lupp approximator + double-double refinement)
- trans-ir (single lupp approximator + dynamic precision refinement)

# Installation 
- Install MKL and XBLAS and put XBLAS path to XBLASPATH in Makefile.
- Edit ${XBLAS}/src/blas_extended_proto.h for cplus extension:

   #ifdef __cplusplus 
   extern "C" {
   #endif
     
     functions list

   #ifdef __cplusplus
   }
   #endif

- Add source commands into your ~/.bashrc file 

e.g., source /home/user/intel/bin/compilervars.sh intel64 
      
      source /home/user/intel/mkl/bin/mklvars.sh intel64

# Execution examples: 

Type make in /DynIR/

./lupp uni 1024 1

Output:
  
  uni-ir (1024x1024 matrix): 
  
  Accuracy: 9.27044e-17
  
  Success
  
  Runtime: 0.142067 secs

This command returns the excution time and runtime for double lupp approximator generation + uni-ir refinement.
Likewise,

./lupp trans 1024 1
  
Output:

  Accuracy: 9.27169e-17
  
  Success
  
  Runtime: 0.117737 secs

./lupp trans 1024 0
  
  trans-ir (1024x1024 matrix): 
  
  Success
  
  Runtime: 0.085209 secs


In trans-IR (e.g., Dynamic Precision Iterative Refinement), one additional iteration due to accuracy check will increase the execution time in trans-ir execution, but this accuracy check is not needed since the accuracy has been guaranteed by numerical analysis (refer to the paper below).

Citations: 
J. Lee, H. Vandierendonck, M. Arif, G. D. Peterson and D. S. Nikolopoulos, "Energy-Efficient Iterative Refinement Using Dynamic Precision," in IEEE Journal on Emerging and Selected Topics in Circuits and Systems, vol. 8, no. 4, pp. 722-735, Dec. 2018, doi: 10.1109/JETCAS.2018.2850665.

Contact: JunKyu Lee (eecejk@gmail.com)
