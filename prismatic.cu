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
#include "cpp_functions.h"
#include "cuda_common.cuh"
#include "cuda_functions.cuh"

int main() {

    //Alocando memoria para atributos das structs
	inicializarStructs();
    
	//Lendo os arquivos de entrada e inicializando os pacotes de simulação.
	pmrdr2_();

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
L103:;
	pmwrt2_(1);
	printf("  Number of simulated showers = %.6E\n", *CNTRL_.SHN);
	plotdose2_();
//Resete da GPU
    cudaDeviceReset();

    //transferindo os structs da CPU para GPU
    transfCPU_to_GPU();

    //teste<<<1,1>>>(); CHAMADA DA FUNCAO KERNEL

    cudaDeviceSynchronize();

    transfGPU_to_CPU();

    memoryFreeGPU();
    memoryFree();
	printf("  *** END ***\n");
	return 0;
}

