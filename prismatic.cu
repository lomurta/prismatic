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
#include "cuda_functions.cuh"
#include "cpp_functions.h"




int main() {

	int sizeTrack = 0;
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
			quickSort(PRITRACK, 0, pilhaPart - 1);

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
			//gpuErrchk(cudaMalloc(&d_TRACK_mod, sizeof(hd_TRACK_MOD)*pilhaPart));
    		//gpuErrchk(cudaMemcpy(d_TRACK_mod, PRITRACK, sizeof(hd_TRACK_MOD)*pilhaPart, cudaMemcpyHostToDevice));
    		//gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, d_TRACK_mod, sizeof(hd_TRACK_MOD*)*pilhaPart));
			gpuErrchk(cudaMemcpyToSymbol(dg_TRACK_mod_, PRITRACK, sizeof(hd_TRACK_MOD)*pilhaPart,0));
			
			//Quantidade de blocos no grid e de threads nos blocos
			dim3 block(blockSize);
			dim3 grid(ceil(sizeTrack / block.x)+1);

			/*printf("[0]: %lf\n", PRITRACK[0].E);
			printf("[1]: %lf\n", PRITRACK[1].E);
			printf("[256]: %lf\n", PRITRACK[256].E);
			printf("[4095]: %lf\n", PRITRACK[4095].E);*/
			//printf("[256]: %lf\n", PRITRACK[4096].E);
			//chamada do kernel para simulacao das particulas primarias
			showers_pri<<<grid,block>>>(sizeTrack);
			//Aguarda o termino da simulação das particulas primarias enviadas
			gpuErrchk(cudaDeviceSynchronize());
			//gpuErrchk(cudaFree(d_TRACK_mod));

			//resgata o pacote de particulas primarias da gpu
			transfSecTracksGPU_to_CPU();
			
		
			//printf("\nquantidade de parricula secundaria photon %d\n\n", nTRACKS_.nSECTRACK_G);
			//printf("\nquantidade de parricula secundaria eletron %d\n\n", nTRACKS_.nSECTRACK_E);
			//printf("\nquantidade de parricula secundaria positron %d\n\n", nTRACKS_.nSECTRACK_P);
			//Esvaziando parcialmente particulas secundarias.
			/*if (nTRACKS_.nSECTRACK_E >= pilhaPart/2){
				while (nTRACKS_.nSECTRACK_E > pilhaPart/4){
					simSecTrack_E();
				}
			}
			if (nTRACKS_.nSECTRACK_G >= pilhaPart/2){
				while (nTRACKS_.nSECTRACK_G > pilhaPart/4){
					simSecTrack_G();
				}
			}
			if (nTRACKS_.nSECTRACK_P >= pilhaPart/2){
				while (nTRACKS_.nSECTRACK_P > pilhaPart/4){
					simSecTrack_P();
				}
			}*/

			//Zerando particulas segundarias
			while ((nTRACKS_.nSECTRACK_E > 0) || (nTRACKS_.nSECTRACK_G > 0) || (nTRACKS_.nSECTRACK_P > 0)){
			if (nTRACKS_.nSECTRACK_E > 0)
				simSecTrack_E();
			if (nTRACKS_.nSECTRACK_G > 0)
				simSecTrack_G();
			if (nTRACKS_.nSECTRACK_P > 0)
				simSecTrack_P();
			}
			gpuErrchk(cudaDeviceSynchronize());

			sizeTrack = pilhaPart;
			//Quantidade de blocos no grid e de threads nos blocos
			dim3 blockCont(blockSize);
			dim3 gridCont(ceil(sizeTrack / block.x)+1);
			
			showers_cont<<<gridCont,blockCont>>>(sizeTrack);
			//Aguarda o termino da simulação das particulas primarias enviadas
			gpuErrchk(cudaDeviceSynchronize());

			printf("Simulado: %f\n", *CNTRL_.SHN);


			
		/*	if ((nTRACKS_.nSECTRACK_E >= pilhaSec/2) || (nTRACKS_.nSECTRACK_G >= pilhaSec/2) || (nTRACKS_.nSECTRACK_P >= pilhaSec/2)){
				while ((nTRACKS_.nSECTRACK_E > 0) || (nTRACKS_.nSECTRACK_G > 0) || (nTRACKS_.nSECTRACK_P > 0 )){
				/*	printf("Antes da simulacao\n");
					printf("Quantidade de parricula secundaria photon: %d\n", nTRACKS_.nSECTRACK_G);
					printf("Quantidade de parricula secundaria eletron: %d\n", nTRACKS_.nSECTRACK_E);
					printf("Quantidade de parricula secundaria positron: %d\n\n", nTRACKS_.nSECTRACK_P);*/
				////	if (nTRACKS_.nSECTRACK_E > 0){
				//		simSecTrack_E();
				/*	printf("Apos simular Eletrons\n");
					printf("Quantidade de parricula secundaria photon: %d\n", nTRACKS_.nSECTRACK_G);
					printf("Quantidade de parricula secundaria eletron: %d\n", nTRACKS_.nSECTRACK_E);
					printf("Quantidade de parricula secundaria positron: %d\n\n", nTRACKS_.nSECTRACK_P);*/
				//	}
				//	if (nTRACKS_.nSECTRACK_G > 0){
				//		simSecTrack_G();
				/*	printf("Apos simular Fotons\n");
					printf("Quantidade de parricula secundaria photon: %d\n", nTRACKS_.nSECTRACK_G);
					printf("Quantidade de parricula secundaria eletron: %d\n", nTRACKS_.nSECTRACK_E);
					printf("Quantidade de parricula secundaria positron: %d\n\n", nTRACKS_.nSECTRACK_P);*/
					//}
				//	if (nTRACKS_.nSECTRACK_P > 0){
				//		simSecTrack_P();
				/*	printf("Apos simular Positrons\n");
					printf("Quantidade de parricula secundaria photon: %d\n", nTRACKS_.nSECTRACK_G);
					printf("Quantidade de parricula secundaria eletron: %d\n", nTRACKS_.nSECTRACK_E);
					printf("Quantidade de parricula secundaria positron: %d\n\n", nTRACKS_.nSECTRACK_P);*/
				//	}
			//	}
		//	}*/

			

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

		//Zerando particulas segundarias
		/*while ((nTRACKS_.nSECTRACK_E > 0) || (nTRACKS_.nSECTRACK_G > 0) || (nTRACKS_.nSECTRACK_P > 0 )){
			if (nTRACKS_.nSECTRACK_E > 0)
				simSecTrack_E();
			if (nTRACKS_.nSECTRACK_G > 0)
				simSecTrack_G();
			if (nTRACKS_.nSECTRACK_P > 0)
				simSecTrack_P();
		}*/

L202:;
		
		*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
		//retorna os dados da GPU para imprimir o DUMP
		transfGPU_to_CPU();
		gpuErrchk(cudaDeviceSynchronize());
		memoryFreeGPU();
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
	pmwrt2_(1);
	
	plotdose2_();
	memoryFreeCPU();
	printf("  *** END ***\n");
	return 0;
	
}

