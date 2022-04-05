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

	int simGPU = 1;

    //Alocando memoria para atributos das structs
	inicializarStructs();
    
	//Lendo os arquivos de entrada e inicializando os pacotes de simulação.
	pmrdr2_();


	if (simGPU){//Simulação na GPU

		if (*CSOUR0_.JOBEND != 0)
			goto L103;
		//Resete da GPU
		cudaDeviceReset();

		//aloca vetores das particulas primarias e secundarias

	
		PRITRACK =  (hd_TRACK_MOD*)malloc(pilhaPart * sizeof(hd_TRACK_MOD)); //vetor de particulas primarias
		SECTRACK_G = (hd_TRACK_MOD*)malloc(pilhaPart * sizeof(hd_TRACK_MOD)); //vetor de particulas secundarias de fotons
		SECTRACK_E = (hd_TRACK_MOD*)malloc(pilhaPart * sizeof(hd_TRACK_MOD)); //vetor de particulas secudarias de eletrons
		SECTRACK_P = (hd_TRACK_MOD*)malloc(pilhaPart * sizeof(hd_TRACK_MOD)); //vetor de particulas secundarias de protons


		//transferindo os structs da CPU para GPU
		transfCPU_to_GPU();

		while ((*CNTRL_.TSEC < *CNTRL_.TSECA) && (*CNTRL_.SHN < *CNTRL_.DSHN)){

			//criar vetor de particulas primarias
			iniPRITRACK(); 
			//transfere a pilha de particulas secundarias para GPU
			transfSecTracksCPU_to_GPU();
			//chamada do kernel para simulacao das particulas primarias
		//	showers<<<1,1>>>();
			//Aguarda o termino da simulação das particulas primarias enviadas
			cudaDeviceSynchronize();
			//resgata o pacote de particulas primarias da gpu
			transfSecTracksGPU_to_CPU();
			//while para realizar simulação das particulas secundarias
			//verifica se o pacote de particulas secundarias esta grande o suficinete para uma simulacao
			//Se for grande o suficiente, ondena o vetor de particulas secundarias
			//chama a simulacao das particulas secundarias se for grande o suficiente	
			//	showers<<<1,1>>>();

			//verifica tempo do DUMP
			if (*CDUMP_.LDUMP) {
				if (*CNTRL_.TSEC - *CNTRL_.TSECAD > *CNTRL_.DUMPP) {
					//retorna os dados da GPU para imprimir o DUMP
					transfGPU_to_CPU();
					cudaDeviceSynchronize();
					*CNTRL_.TSIM = *CNTRL_.TSIM + cputim2_() - *CNTRL_.CPUT0;
					pmwrt2_(-1);
					printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
					*CNTRL_.TSECAD = *CNTRL_.TSEC;
					*CNTRL_.CPUT0 = cputim2_();
				}
			}
		}
	
		memoryFreeGPU();
		free(PRITRACK);
		free(SECTRACK_G);
		free(SECTRACK_E);
		free(SECTRACK_P);


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
		transfCPU_to_GPU();
		cudaDeviceSynchronize();
		transfGPU_to_CPU();
		memoryFreeGPU();
		cudaDeviceSynchronize();
		memoryFree();
		printf("  *** END ***\n");
		return 0;
	

}

