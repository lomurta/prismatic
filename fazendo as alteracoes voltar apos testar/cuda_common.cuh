#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

__global__ void scan_efficient_1G(int * input, int* auxiliry_array, int input_size);
__global__ void scan_summation(int * input, int * auxiliry_array, int input_size);

#endif // !CUDA_COMMON_H

//void query_device();

__device__ int atomicSub2(int *address, int val)
{
	int *teste = (int *)address;
	int old = *address;
	int assumed;
	do
	{
		assumed = old;
		old = atomicCAS(teste, assumed, (-val) + assumed);
	} while (assumed != old);
	return old;
}

__device__ double atomicAdd2(double *address, double val)
{
	unsigned long long int *address_as_ull =
		(unsigned long long int *)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
						__double_as_longlong(val +
											 __longlong_as_double(assumed)));
		// Note: uses integer comparison to avoid hang in case of NaN (since NaN !=NaN)
	} while (assumed != old);
	return __longlong_as_double(old);
}

__device__ int atomicAdd2(int *address, int val)
{
	int *teste = (int *)address;
	int old = *address;
	int assumed;
	do
	{
		assumed = old;
		old = atomicCAS(teste, assumed, val + assumed);
	} while (assumed != old);
	return old;
}



