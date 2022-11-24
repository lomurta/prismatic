#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <stdint.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include "cuda_runtime.h"
#include <curand.h>
#include <curand_kernel.h>
#include "device_launch_parameters.h"
using namespace std;
#include "hd_structs.cuh"
#include "ISEEDS.h"
#include "cuda_common.cuh"
#include "cuda_functions.cuh"
#include "cpp_functions.h"


int main() {

	//Alocando memoria para atributos das structs
	memoryAllocCPU();

	
	int simGPU = 1;

	//Lendo os arquivos de entrada e inicializando os pacotes de simulação.
	pmrdr2_();

	if (simGPU){//Simulação na GPU

		int sizeTrack = 0;

		if (*CSOUR0_.JOBEND != 0)
			goto L103;

		//alocando memooria na GPU
		memoryAllocGPU();

		//Inicializa as semnentes para GPU e par CPU
		initializeISSEDS_();
		*RSEED_.ISEED1 = IS1[0];
		*RSEED_.ISEED2 = IS2[0];

		sizeTrack = pilhaPart;

		dim3 block(blockSize);
		dim3 grid(ceil(sizeTrack / block.x)+1);

		//inicializa gerando de numeros aleatorios cuRand
		initializeRand<<<grid, block>>>(sizeTrack);
		gpuErrchk(cudaDeviceSynchronize());

		//trasnferindo dados para GPU
		transfCPU_to_GPU();

		while ((*CNTRL_.TSEC < *CNTRL_.TSECA) && (*CNTRL_.SHN < *CNTRL_.DSHN)){
			
			cleans2GPU_();
			simPriTrack_G();
			transfnTRACKSGPU_to_CPU();

			// simulação de particulas secundarias
			while ((nTRACKS_.nSECTRACK_E > 0) || (nTRACKS_.nSECTRACK_G > 0) || (nTRACKS_.nSECTRACK_P > 0))
			{
				if (nTRACKS_.nSECTRACK_E > 0)
				{
					simSecTrack_E();
				}

				if (nTRACKS_.nSECTRACK_P > 0)
				{
					simSecTrack_P();
				}

				if (nTRACKS_.nSECTRACK_G > 0)
				{
					simSecTrack_G();
				}
			}
			gpuErrchk(cudaDeviceSynchronize());

			//contabilizacao da contribuição das particulas para a dose geral
			showers_cont<<<grid,block>>>(sizeTrack);
			gpuErrchk(cudaDeviceSynchronize());

			printf("Simulado: %f\n", *CNTRL_.SHN);
			
			if (*CSOUR0_.JOBEND != 0)
				goto L202;

			timer2_(*CNTRL_.TSEC);

			//verifica tempo do DUMP
			//if (*CDUMP_.LDUMP) {
				if ((*CNTRL_.TSEC - *CNTRL_.TSECAD > *CNTRL_.DUMPP) || (*CNTRL_.SHN == 1e4) || (*CNTRL_.SHN == 1e5) || (*CNTRL_.SHN == 1e6) || (*CNTRL_.SHN == 1e7) || (*CNTRL_.SHN == 1e8) || (*CNTRL_.SHN == 1e9)) {
				
					//retorna os dados da gpu para cpu
					transfGPU_to_CPU();
					gpuErrchk(cudaDeviceSynchronize());

					*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;

					//imprime os resultados parciais da simulação
					pmwrt2_(1);
					plotdose2_();

					*CNTRL_.TSECAD = *CNTRL_.TSEC;
					*CNTRL_.CPUT0 = cputim2_();

					printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
				}
			//}
		}


L202:;
	
		*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
		
		//retorna os dados da GPU para imprimir o DUMP
		transfGPU_to_CPU();
		
	}else{ //Simulação na CPU

		if (*CSOUR0_.JOBEND != 0)
			goto L103;

L101:;
		//Simulação de uma nova ducha e pontuação.
		shower2_();
		if (*CSOUR0_.JOBEND != 0)
			goto L102;

		timer2_(*CNTRL_.TSEC);

		//Terminar a simulação após o tempo previsto ou após completar Chuveiros DSHN.
		if ((*CNTRL_.TSEC < *CNTRL_.TSECA) && (*CNTRL_.SHN < *CNTRL_.DSHN)) {
			//Escreva os resultados parciais após cada período de despejo.
			if (*CDUMP_.LDUMP) {
				if (*CNTRL_.TSEC - *CNTRL_.TSECAD > *CNTRL_.DUMPP) {
					*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
					pmwrt2_(-1);
					printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
					*CNTRL_.TSECAD = *CNTRL_.TSEC;
					*CNTRL_.CPUT0 = cputim2_();
					goto L101;
				}
			}
			goto L101;
		}
L102:;
		*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
	}

L103:;//Imprimir resultados Finais
	printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
	if (simGPU == 0){
		pmwrt2_(1);
		plotdose2_();
	}else{
		memoryFreeGPU();
		gpuErrchk(cudaDeviceReset());
	}
	memoryFreeCPU();
	printf("  *** END ***\n");
	return 0;
}

