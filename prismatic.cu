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
#include "device_launch_parameters.h"
using namespace std;
#include "hd_structs.cuh"
#include "cuda_common.cuh"
#include "cpp_functions.h"
#include "cuda_functions.cuh"



int main() {

	int sizeTrack = 0;
	int block_size = 128;
	int simGPU = 1;

    //Alocando memoria para atributos das structs
	memoryAllocCPU();
    
	//Lendo os arquivos de entrada e inicializando os pacotes de simulação.
	pmrdr2_();

	if (simGPU){//Simulação na GPU

		if (*CSOUR0_.JOBEND != 0)
			goto L103;
		//Resete da GPUFF
		//gpuErrchk(cudaDeviceReset());

		//alocando memooria na GPU
		memoryAllocGPU();

		//aloca vetores das particulas primarias e secundarias
		bool btransfCPU_to_GPU = false;
		//bool btransfGPU_to_CPU = false;
		cleans2GPU_();
		
		while ((*CNTRL_.TSEC < *CNTRL_.TSECA) && (*CNTRL_.SHN < *CNTRL_.DSHN)){
			
			//criar vetor de particulas primarias inicial
			iniPRITRACK(); 

			//transferindo os structs da CPU para GPU
			if (!btransfCPU_to_GPU){
				transfCPU_to_GPU();
				btransfCPU_to_GPU = true;
			}
			//seta o tamanho da pilha a ser simulada
			sizeTrack = pilhaPart;

			//transfere a pilha de particulas secundarias para GPU
			transfSecTracksCPU_to_GPU();

			//transfere as particulas a serem simuladas para a GPU
			gpuErrchk(cudaMalloc(&d_TRACK_mod, sizeof(hd_TRACK_MOD)*pilhaPart));
    		//gpuErrchk(cudaMemcpy(d_TRACK_mod, PRITRACK, sizeof(hd_TRACK_MOD)*pilhaPart, cudaMemcpyHostToDevice));
    		//gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, d_TRACK_mod, sizeof(hd_TRACK_MOD*)*pilhaPart));
			gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, PRITRACK, sizeof(hd_TRACK_MOD)*pilhaPart,0));
			
			//Quantidade de blocos no grid e de threads nos blocos
			dim3 block(block_size);
			dim3 grid((sizeTrack / block.x));

			printf("[0]: %lf\n", PRITRACK[0].E);
			printf("[1]: %lf\n", PRITRACK[1].E);
			printf("[256]: %lf\n", PRITRACK[256].E);
			printf("[4095]: %lf\n", PRITRACK[4095].E);
			//printf("[256]: %lf\n", PRITRACK[4096].E);
			//chamada do kernel para simulacao das particulas primarias
			showers_pri<<<grid,block>>>(sizeTrack);
			//Aguarda o termino da simulação das particulas primarias enviadas
			gpuErrchk(cudaDeviceSynchronize());
			gpuErrchk(cudaFree(d_TRACK_mod));

			//resgata o pacote de particulas primarias da gpu
			transfSecTracksGPU_to_CPU();
		
			printf("\nquantidade de parricula secundaria photon %d\n\n", nTRACKS_.nSECTRACK_G);
			printf("\nquantidade de parricula secundaria eletron %d\n\n", nTRACKS_.nSECTRACK_E);
			printf("\nquantidade de parricula secundaria positron %d\n\n", nTRACKS_.nSECTRACK_P);
			
			if (nTRACKS_.nSECTRACK_E == pilhaSec){
				while (nTRACKS_.nSECTRACK_E > 1024){
					nTRACKS_.nSECTRACK_E = 0;
					sizeTrack = pilhaSec;

					//transfere a pilha de particulas secundarias para GPU
					transfSecTracksCPU_to_GPU();
					
					
					//transfere as particulas a serem simuladas para a GPU
					gpuErrchk(cudaMalloc(&d_TRACK_mod, sizeof(hd_TRACK_MOD)*pilhaSec));
					//PRITRACK = SECTRACK_E;
					//gpuErrchk(cudaMemcpy(d_TRACK_mod, SECTRACK_E, sizeof(hd_TRACK_MOD)*pilhaSec, cudaMemcpyHostToDevice));
					//gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, d_TRACK_mod, sizeof(hd_TRACK_MOD)*pilhaSec));

					//memcpy(vTrack_Simular, SECTRACK_E, sizeof(hd_TRACK_MOD)*pilhaPart);

					gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, SECTRACK_E, sizeof(hd_TRACK_MOD)*pilhaPart,0));

					
				

					dim3 blockSec(block_size);
					dim3 gridSec((sizeTrack / block.x));

					printf("[0]: %lf\n", SECTRACK_E[0].E);
					printf("[1]: %lf\n", SECTRACK_E[1].E);
					printf("[256]: %lf\n", SECTRACK_E[256].E);
					printf("[4095]: %lf\n", SECTRACK_E[4095].E);
					printf("[4096]: %lf\n", SECTRACK_E[4096].E);

					//chamada do kernel para simulacao das particulas primarias
					showers_sec<<<gridSec,blockSec>>>(sizeTrack);
						//Aguarda o termino da simulação das particulas primarias enviadas
					gpuErrchk(cudaDeviceSynchronize());
						//resgata o pacote de particulas primarias da gpu
					gpuErrchk(cudaFree(d_TRACK_mod));

					transfSecTracksGPU_to_CPU();
				}
			}

			if (*CSOUR0_.JOBEND != 0)
				goto L202;

			timer2_(*CNTRL_.TSEC);

			//verifica tempo do DUMP
			if (*CDUMP_.LDUMP) {
				if (*CNTRL_.TSEC - *CNTRL_.TSECAD > *CNTRL_.DUMPP) {
					//retorna os dados da GPU para imprimir o DUMP
					
					transfGPU_to_CPU();
					gpuErrchk(cudaDeviceSynchronize());
					*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
					pmwrt2_(-1);
					printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
					*CNTRL_.TSECAD = *CNTRL_.TSEC;
					*CNTRL_.CPUT0 = cputim2_();
				}
			}
		}

L202:;
		
		*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
		//retorna os dados da GPU para imprimir o DUMP
		transfGPU_to_CPU();
		gpuErrchk(cudaDeviceSynchronize());
		memoryFreeGPU();
		printf("aqui\n");
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
	pmwrt2_(1);
	printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
	plotdose2_();
	memoryFreeCPU();
	printf("  *** END ***\n");
	return 0;
	
}

