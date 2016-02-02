# raytracing_xray
This pack includes the code to investigate different APIs for CPU parallel computing in forward projection step for X-Ray CT. 

Here, we looked at two possible methods for forward projection:
  - Ray tracing: This is an online method that is flexible with different source-detector coordinates and is popular in medical imaging community. For the pseudocode, please refer to document url provided at the end.
  - Sparse matrix vector multiplication (SMVM): This is an offline method where each element of the system matrix is pre-stored and used. 

In both cases, essentially the operation is H*x, where x is the image domain variables and H is the system matrix.

The ray-tracing forward projection here assumes point source and detectors where the corresponding x-y coordinates for each are read through a .txt file initially.

Here is a short description about the files included:
  - ray_tracing_base.cpp: A basic case with no parallelism used. A good start to explore the algorithm.
  - ray_tracing_cilk.cpp: This version uses Cilk (http://www.cilkplus.org/) library to exploit parallelism.
  - ray_tracing_omp.cpp: This version uses OpenMP (http://openmp.org/wp/) and Advanced Vector Extensions (AVX) (https://software.intel.com/en-us/articles/introduction-to-intel-advanced-vector-extensions) to exploit the parallelism.
  - smvm_base_main.cpp: A basic case of SMVM.
  - smvm_cilk.cpp: SMVM with Cilk.
  - smvm_omp.cpp: SMVM with OpenMP.
  - smvm_omp.cpp SMVM with OpenMP, tweaked version.
