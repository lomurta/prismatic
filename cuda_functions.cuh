
__device__ double atomicAdd2(double *address, double val);

__device__ int atomicAdd2(int *address, int val);

__device__ void d_gcone02_(double THETA, double PHI, double ALPHA);

__device__ void d_gcone2_(double &UF, double &VF, double &WF);

__device__ void d_step2_(double &DS, double &DSEF, int &NCROSS);

__global__ void g_step2_(int size);

__device__ void d_steplb2_(int &KB, int &IERR);

//__device__ void d_stepsi2_(int &KB, double *S, int *IS, int &NSC);

__device__ void d_stepsi2_(int &KB, int &NSC);

__device__ void d_fsurf2_(int &KS, double &A, double &B, double &C);

__device__ void d_locate2_();

__device__ void d_sendet2_(double &ED, int &ID);

__device__ void d_secpar2_(int &LEFT);

__device__ void d_pana2_(double &E, double &E1, double &CDT1, double &E2, double &CDT2, int &M);

__device__ void d_paux2_();

__device__ void d_tenang2_(int &IEXIT, int &N);

__device__ double d_rand2_(double DUMMY);

__device__ void d_schiff2_(double &B, double &G1, double &G2);

__device__ void d_gaux2_();

__device__ void d_peld2_(double &RNDC, double &RMU);

__device__ void d_pina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC);

__device__ void d_psia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

__device__ void d_gpha2_(double &ES, int &IZZ, int &ISH);

__device__ void d_sauter2_(double &ES, double &CDTS);

__device__ void d_gppa2_(double &EE, double &CDTE, double &EP, double &CDTP, int &IZZ, int &ISH);

__device__ void d_eaux2_();

__device__ void d_graa2_(double &E, double &CDT, int &IEFF, int &M);

__device__ void d_dirpol2_(double &CDT, double &DF, double &CONS, double &SP1, double &SP2, double &SP3, double &U, double &V, double &W);

__device__ void d_gcoa2_(double &E, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

__device__ void d_ebraa2_(double &E, double &DE, double &CDT, int &M);

__device__ void d_esia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH);

__device__ void d_relax2_(int &IZ, int &IS);

__device__ void d_stores2_(double &EI, double &XI, double &YI, double &ZI, double &UI, double &VI, double &WI, double &WGHTI, int &KPARI, int *ILBI, int &IPOLI);

__device__ void d_eeld2_(double &RNDC, double &RMU);

__device__ void d_eina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC);

__device__ void d_ebra2_(double &E, double &W, int &M);

__device__ void d_fimdet2_(int &N, int &ID, double &DSEF);

__device__ void d_fimdes2_(int &N, int &ID, double &EI, double &DECSD, double &DSEF);

__device__ void d_knock2_(double &DE, int &ICOL);

__device__ void d_eela2_(double &A, double &B, double &RNDC, double &RMU);

__device__ void d_direct2_(double &CDT, double &DF, double &U, double &V, double &W);

__device__ void d_panar2_(double &ECUT);

__device__ void d_pimfp2_(int &IEND);

__device__ void d_eimfp2_(int &IEND);

__device__ void d_gimfp2_();

__device__ void d_start2_();

__device__ void d_jump2_(double &DSMAX, double &DS);

__device__ double d_rndg32_();

__device__ void d_gcone2_(double &UF, double &VF, double &WF);

__device__ void d_simdet2_(int &N, int &ID);

__device__ void d_sdose2_(double &DEP, double &XD, double &YD, double &ZD, int &MATC, int &N);

__device__ void d_cleans2_();

__global__ void showers_sec_G(int size);

__global__ void showers_sec_E(int size);

__global__ void showers_sec_P(int size);

__global__ void showers_cont(int size);

__device__ void d_knock2_G(double &DE, int &ICOL);

__global__ void g_knock2_G(int size);

__device__ void d_knock2_E(double &DE, int &ICOL);

__global__ void g_knock2_E(int size);

__device__ void d_knock2_P(double &DE, int &ICOL);

__global__ void g_knock2_P(int size);

__device__ void showers_step1_G(int size); // L302

__device__ void d_jump2_G(double &DSMAX, double &DS);

__global__ void g_jump2_G(int size);

__global__ void dg_jump2_G(double &DSMAX, double &DS);

__device__ void d_jump2_E(double &DSMAX, double &DS);

__global__ void g_jump2_E(int size);

__device__ void d_jump2_P(double &DSMAX, double &DS);

__global__ void g_jump2_P(int size);

// tentativa de sequencia

__global__ void g_showers_step1(int size); // primeiro passo

__device__ void d_showers_step2_G(int size); // segundo passo

__device__ void d_showers_step2_E(int size); // segundo passo

__device__ void d_showers_step2_P(int size); // segundo passo

__global__ void g_showers_step3_G(int size); // terceiro passo

__global__ void g_showers_step3_E(int size); // terceiro passo

__global__ void g_showers_step3_P(int size); // terceiro passo

__global__ void g_showers_step4_G(int size); // contadores de particulas

__global__ void g_showers_step4_E(int size); // contadores de particulas

__global__ void g_showers_step4_P(int size); // contadores de particulas

__global__ void g_showers_step5_G(int size);

__global__ void g_showers_step5_E(int size);

__global__ void g_showers_step5_P(int size);

__global__ void g_showers_step6(int size); // contador de particulas

__global__ void g_showers_sec_G(int size);

__global__ void g_showers_sec_E(int size);

__global__ void g_showers_sec_P(int size);

//	__global__ void g_showers_step6_E(int size); // particulas secundarias eletrons

//	__global__ void g_showers_step6_P(int size); // particulas secundarias positrons

__global__ void g_showers_step8(int size); // contabilizadores de deposito de dose


__global__ void g_cpyTrack_simular(int tipo, int size);

__global__ void g_cpyTrack_complementar(int tipo, int size);

__global__ void g_imprimir(int size);






__global__ void g_showers_step1(int size)
{ // Inicio da simulacao das particulas primarias
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if ((index < size) && (dg_TRACK_mod_[index].STEP == 1))
	{
		dg_TRACK_mod_[index].INDEX = index;
		dg_TRACK_mod_[index].STEP = 2; //proximo passo 2: d_showers_step2_G
		d_showers_step2_G(size);
	}
}

__device__ void d_showers_step2_G(int size)
{ // segundo passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// Verifique se a trajetória cruza o sistema de materiais.

	// L302:;

	if (dg_TRACK_mod_[index].STEP == 2)
	{
		dg_TRACK_mod_[index].STEP = 3; //g_showers_stepd_G

		dg_wSHOWERS_[trackIndex].IEXIT = 0;

		d_locate2_();
		if (dg_TRACK_mod_[index].MAT == 0)
		{
			dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
			dg_wSHOWERS_[trackIndex].DS = 1.0e30;
			d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS);
			/*if (dg_TRACK_mod_[index].LAGE)
				DPAGE(dg_wSHOWERS_[trackIndex].DSEF,DSTOT)*/
			// Funcao para contabilizar o tempo de vida de uma particula, não sera implementada

			if (dg_TRACK_mod_[index].MAT == 0)
			{ // A particula não entrou no sistema
				if (dg_TRACK_mod_[index].W > 0)
				{
					dg_wSHOWERS_[trackIndex].IEXIT = 1; // Rotula partículas ascendentes emergentes.
				}
				else
				{
					dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				dg_TRACK_mod_[index].STEP = 9; //g_showers_step6

				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}

				return;
			}
			// Detetores de Impacto

			dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
			if (dg_wSHOWERS_[trackIndex].IDET != 0)
			{
				if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
				{
					// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
					/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
						NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
						wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
						dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
						dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

					d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

					if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
					{
						atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
						//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
						dg_wSHOWERS_[trackIndex].IEXIT = 3;
						// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
						dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
						dg_TRACK_mod_[index].STEP = 9; //g_showers_step6
						if (d_cudaCapability >= 6)
						{
							atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
						}
						else
						{
							dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
						}
						return;
					}
				}
			}
		}
		// Aniquiliação de positron quando a energia da particula é muito pequena

		if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{ // energia é muito baixa
			dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
			{ // aniquilação de positrion
				d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
				dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
			}
			dg_TRACK_mod_[index].E = 0.0e0;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
			if (dg_CNT6_.LDOSEM)
			{
				d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
			dg_wSHOWERS_[trackIndex].IEXIT = 3;
			// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
			dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
			dg_TRACK_mod_[index].STEP = 9; //g_showers_step6
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
			}
			else
			{
				dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
			}
			return;
		}
		
	}
}

__device__ void d_showers_step2_E(int size)
{ // segundo passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// Verifique se a trajetória cruza o sistema de materiais.

	// L302:;

	dg_wSHOWERS_[trackIndex].IEXIT = 0;

	d_locate2_();
	if (dg_TRACK_mod_[index].MAT == 0)
	{
		dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
		dg_wSHOWERS_[trackIndex].DS = 1.0e30;
		d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS);
		/*if (dg_TRACK_mod_[index].LAGE)
			DPAGE(dg_wSHOWERS_[trackIndex].DSEF,DSTOT)*/
		// Funcao para contabilizar o tempo de vida de uma particula, não sera implementada

		if (dg_TRACK_mod_[index].MAT == 0)
		{ // A particula não entrou no sistema
			if (dg_TRACK_mod_[index].W > 0)
			{
				dg_wSHOWERS_[trackIndex].IEXIT = 1; // Rotula partículas ascendentes emergentes.
			}
			else
			{
				dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
			}
			// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
			dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
			}
			else
			{
				dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
			}

			return;
		}
		// Detetores de Impacto

		dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
		if (dg_wSHOWERS_[trackIndex].IDET != 0)
		{
			if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
			{
				// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
				/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
					NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
					wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
					dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
					dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

				d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

				if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
				{
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
					//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
					dg_wSHOWERS_[trackIndex].IEXIT = 3;
					// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
					dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
					}
					else
					{
						dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
					}
					return;
				}
			}
		}
	}
	// Aniquiliação de positron quando a energia da particula é muito pequena

	if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
	{ // energia é muito baixa
		dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
		if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
		{ // aniquilação de positrion
			d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
			dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
		}
		dg_TRACK_mod_[index].E = 0.0e0;
		atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
		//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
		if (dg_CNT6_.LDOSEM)
		{
			d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
		}
		dg_wSHOWERS_[trackIndex].IEXIT = 3;
		// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
		dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
		if (d_cudaCapability >= 6)
		{
			atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
		}
		else
		{
			dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
		}
		return;
	}
}

__device__ void d_showers_step2_P(int size)
{ // segundo passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// Verifique se a trajetória cruza o sistema de materiais.

	// L302:;

	dg_wSHOWERS_[trackIndex].IEXIT = 0;

	d_locate2_();
	if (dg_TRACK_mod_[index].MAT == 0)
	{
		dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
		dg_wSHOWERS_[trackIndex].DS = 1.0e30;
		d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS);
		/*if (dg_TRACK_mod_[index].LAGE)
			DPAGE(dg_wSHOWERS_[trackIndex].DSEF,DSTOT)*/
		// Funcao para contabilizar o tempo de vida de uma particula, não sera implementada

		if (dg_TRACK_mod_[index].MAT == 0)
		{ // A particula não entrou no sistema
			if (dg_TRACK_mod_[index].W > 0)
			{
				dg_wSHOWERS_[trackIndex].IEXIT = 1; // Rotula partículas ascendentes emergentes.
			}
			else
			{
				dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
			}
			// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
			dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
			}
			else
			{
				dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
			}

			return;
		}
		// Detetores de Impacto

		dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
		if (dg_wSHOWERS_[trackIndex].IDET != 0)
		{
			if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
			{
				// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
				/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
					NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
					wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
					dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
					dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

				d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

				if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
				{
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
					//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
					dg_wSHOWERS_[trackIndex].IEXIT = 3;
					// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
					dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
					}
					else
					{
						dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
					}
					return;
				}
			}
		}
	}
	// Aniquiliação de positron quando a energia da particula é muito pequena

	if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
	{ // energia é muito baixa
		dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
		if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
		{ // aniquilação de positrion
			d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
			dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
		}
		dg_TRACK_mod_[index].E = 0.0e0;
		atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
		//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
		if (dg_CNT6_.LDOSEM)
		{
			d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
		}
		dg_wSHOWERS_[trackIndex].IEXIT = 3;
		// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
		dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
		if (d_cudaCapability >= 6)
		{
			atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
		}
		else
		{
			dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
		}
		return;
	}
}

__global__ void g_showers_step3_G(int size)
{
	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// if ((index == 1) || (index == 2))
	// printf("valor do trackIndex %d\n", trackIndex);
	if ((index < size) && (dg_TRACK_mod_[index].STEP == 3))
	{
	//	if (dg_TRACK_mod_[index].IEXIT == 0)
	//	{
			//	L102:;
			dg_TRACK_mod_[index].STEP = 4; //g_jump2_G
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				dg_TRACK_mod_[index].STEP = 9; //g_showers_step6
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}
			d_start2_(); // Inicia a simulação no meio atual.
			
	//	}
	}
}

__device__ void d_showers_step3_G(int size)
{ // terceiro passo

	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	if (index < size)
	{
		//if (dg_TRACK_mod_[index].IEXIT == 0)
		//{
			//	L102:;
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			d_start2_(); // Inicia a simulação no meio atual.
	//	}
	}
}

__global__ void g_showers_step3_E(int size)
{
	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 3))
	{
		int trackIndex = dg_TRACK_mod_[index].INDEX;
		//if (dg_TRACK_mod_[index].IEXIT == 0)
		//{
			//	L102:;
			dg_TRACK_mod_[index].STEP = 4; //g_jump2_E
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				dg_TRACK_mod_[index].STEP = 9; //g_showers_ste6
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			d_start2_(); // Inicia a simulação no meio atual.
		//}
	}
}

__device__ void d_showers_step3_E(int size)
{ // terceiro passo

	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	// L102:;

	if (index < size)
	{
		int trackIndex = dg_TRACK_mod_[index].INDEX;
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{
			//	L102:;
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			d_start2_(); // Inicia a simulação no meio atual.
		}
	}
}

__global__ void g_showers_step3_P(int size)
{
	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	if (index < size)
	{
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{
			//	L102:;
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			d_start2_(); // Inicia a simulação no meio atual.
		}
	}
}

__device__ void d_showers_step3_P(int size)
{ // terceiro passo

	// Simulação da historia da particula inicia aqui.

	// Divisão de particulas e roleta russa . Apenas para particulas lidas de um arquivoco Phase-Space.
	// Essa versão não irá implementar reduções de variancia e leitura de arquivo Phase-Space.
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	// L102:;
	if (index < size)
	{
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{
			//	L102:;
			dg_wSHOWERS_[trackIndex].LINTF = false;

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // energia é muito baixa
				dg_wSHOWERS_[trackIndex].DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{ // aniquilação de positrion
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]);
					dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DEP + d_TREV * dg_TRACK_mod_[index].WGHT;
				}
				dg_TRACK_mod_[index].E = 0.0e0;

				// A energia é depositada localmente no material.
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
				//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;
				if (dg_CNT6_.LDOSEM)
				{ // Partícula dentro da caixa de dose.
					d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
				}
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			d_start2_(); // Inicia a simulação no meio atual.
		}
	}
}

__global__ void g_showers_step4_G(int size)
{ // terceiro passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 6))
	{
		//if (dg_TRACK_mod_[index].IEXIT == 0)
		//{
			 dg_TRACK_mod_[index].STEP = 7; //g_knock2

			// Comprimento do caminho livre para o próximo evento de interação.
			//	L103:;
			dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
			dg_wSHOWERS_[trackIndex].MATL = dg_TRACK_mod_[index].MAT;
			dg_wSHOWERS_[trackIndex].XL = dg_TRACK_mod_[index].X;
			dg_wSHOWERS_[trackIndex].YL = dg_TRACK_mod_[index].Y;
			dg_wSHOWERS_[trackIndex].ZL = dg_TRACK_mod_[index].Z;

			/*if ((CFORCI_.LFORCE[KPAR-1][IBODY-1]) && ((WGHT > WLOW[KPAR-1][IBODY-1]) && (WGHT < WHIG[KPAR-1][IBODY-1]))){
				//JUMPF(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  //Força de Interação
				//dg_wSHOWERS_[trackIndex].LINTF=true;
				//Esse trecho nao sera implementado pois a versão nao irá tratar redução de variancia.
			}else`{
				jumtp2(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  ! Analogue simulation.
				dg_wSHOWERS_[trackIndex].LINTF=false
			}*/

			// d_jump2_G(dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DS);

			// d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS); // Determina a posição final do passo.

			// volta aqui
			dg_wSHOWERS_[trackIndex].LINTF = false;

			// Distribuição de energia de fluência.
			dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1];
			if (dg_wSHOWERS_[trackIndex].IDET != 0)
			{
				if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 2)
				{
					if (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1)
					{
						if (dg_TRACK_mod_[index].KPAR == 2)
						{
							d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
						}
						else
						{
							dg_wSHOWERS_[trackIndex].DECSD = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
							if (dg_wSHOWERS_[trackIndex].DECSD > 1.0e-12)
							{																																										   // A distribuição pode se estender
								d_fimdes2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_PENELOPE_mod_.E0STEP[trackIndex], dg_wSHOWERS_[trackIndex].DECSD, dg_wSHOWERS_[trackIndex].DSEF); // abaixo do EABS.
							}
							else
							{
								d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
							}
						}
					}
				}
			}
			

			 dg_wSHOWERS_[trackIndex].CROSS = false;
			

			// A partícula cruzou uma interface.

			if (dg_wSHOWERS_[trackIndex].NCROSS > 0)
			{
				dg_wSHOWERS_[trackIndex].CROSS = true;
				dg_TRACK_mod_[index].STEP = 3; //g_showers_step3_G/E/P
				// Correção da perda de energia suave (CSDA).
				if (dg_TRACK_mod_[index].KPAR != 2)
				{
					dg_TRACK_mod_[index].E = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
					if (dg_CJUMP1_[trackIndex].MHINGE == 0)
					{
						dg_wSHOWERS_[trackIndex].DEP = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							dg_wSHOWERS_[trackIndex].DSEFR = d_rand2_(8.0e0) * dg_wSHOWERS_[trackIndex].DSEF;
							dg_wSHOWERS_[trackIndex].XD = dg_wSHOWERS_[trackIndex].XL + dg_TRACK_mod_[index].U * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].YD = dg_wSHOWERS_[trackIndex].YL + dg_TRACK_mod_[index].V * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].ZD = dg_wSHOWERS_[trackIndex].ZL + dg_TRACK_mod_[index].W * dg_wSHOWERS_[trackIndex].DSEFR;
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XD, dg_wSHOWERS_[trackIndex].YD, dg_wSHOWERS_[trackIndex].ZD, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					else
					{
						dg_wSHOWERS_[trackIndex].DEP = -dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_wSHOWERS_[trackIndex].DS - dg_wSHOWERS_[trackIndex].DSEF) * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XL, dg_wSHOWERS_[trackIndex].YL, dg_wSHOWERS_[trackIndex].ZL, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
					//dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] = dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] + dg_wSHOWERS_[trackIndex].DEP;
				}
				// Verifique se a partícula está fora do invólucro.
				if (dg_TRACK_mod_[index].MAT == 0)
				{ // A partícula está fora do recinto.
					if (dg_TRACK_mod_[index].W > 0.0e0)
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 1; // Marca partículas ascendentes emergentes.
					}
					else
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
					}
					dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
					dg_TRACK_mod_[index].STEP = 9;

					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
					}
					else
					{
						dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
					}
					// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //Saida
					return;
				}

				// Detectores de Impacto

				dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
				if (dg_wSHOWERS_[trackIndex].IDET != 0)
				{
					if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
					{
						// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
						/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
							NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
							wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
							dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
							dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

						d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

						if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
						{
							atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
							//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
							dg_wSHOWERS_[trackIndex].IEXIT = 3;
							dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
							dg_TRACK_mod_[index].STEP = 9;
							if (d_cudaCapability >= 6)
							{
								atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
							}
							else
							{
								dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
							}
							// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
							return;
						}
					}
				}
				
				//d_showers_step3_G(size);

				// if (dg_TRACK_mod_[index].IEXIT == 0)
				// goto L103;
			}
		//}
	}
}

__global__ void g_showers_step4_E(int size)
{ // terceiro passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 6))
	{
		int trackIndex = dg_TRACK_mod_[index].INDEX;
		//if (dg_TRACK_mod_[index].IEXIT == 0)
		//{
			dg_TRACK_mod_[index].STEP = 7; //g_knock2_E

			// Comprimento do caminho livre para o próximo evento de interação.
			//		L103:;
			dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
			dg_wSHOWERS_[trackIndex].MATL = dg_TRACK_mod_[index].MAT;
			dg_wSHOWERS_[trackIndex].XL = dg_TRACK_mod_[index].X;
			dg_wSHOWERS_[trackIndex].YL = dg_TRACK_mod_[index].Y;
			dg_wSHOWERS_[trackIndex].ZL = dg_TRACK_mod_[index].Z;

			/*if ((CFORCI_.LFORCE[KPAR-1][IBODY-1]) && ((WGHT > WLOW[KPAR-1][IBODY-1]) && (WGHT < WHIG[KPAR-1][IBODY-1]))){
				//JUMPF(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  //Força de Interação
				//dg_wSHOWERS_[trackIndex].LINTF=true;
				//Esse trecho nao sera implementado pois a versão nao irá tratar redução de variancia.
			}else`{
				jumtp2(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  ! Analogue simulation.
				dg_wSHOWERS_[trackIndex].LINTF=false
			}*/

			// d_jump2_G(dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DS);

			// d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS); // Determina a posição final do passo.

			// volta aqui
			dg_wSHOWERS_[trackIndex].LINTF = false;

			// Distribuição de energia de fluência.
			dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1];
			if (dg_wSHOWERS_[trackIndex].IDET != 0)
			{
				if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 2)
				{
					if (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1)
					{
						if (dg_TRACK_mod_[index].KPAR == 2)
						{
							d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
						}
						else
						{
							dg_wSHOWERS_[trackIndex].DECSD = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
							if (dg_wSHOWERS_[trackIndex].DECSD > 1.0e-12)
							{																																										   // A distribuição pode se estender
								d_fimdes2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_PENELOPE_mod_.E0STEP[trackIndex], dg_wSHOWERS_[trackIndex].DECSD, dg_wSHOWERS_[trackIndex].DSEF); // abaixo do EABS.
							}
							else
							{
								d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
							}
						}
					}
				}
			}
			 dg_wSHOWERS_[trackIndex].CROSS = false;
			//  A partícula cruzou uma interface.

			if (dg_wSHOWERS_[trackIndex].NCROSS > 0)
			{
				dg_wSHOWERS_[trackIndex].CROSS = true;
				dg_TRACK_mod_[index].STEP = 3; //g_showers_step3_G/E/P
				// Correção da perda de energia suave (CSDA).
				if (dg_TRACK_mod_[index].KPAR != 2)
				{
					dg_TRACK_mod_[index].E = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
					if (dg_CJUMP1_[trackIndex].MHINGE == 0)
					{
						dg_wSHOWERS_[trackIndex].DEP = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							dg_wSHOWERS_[trackIndex].DSEFR = d_rand2_(8.0e0) * dg_wSHOWERS_[trackIndex].DSEF;
							dg_wSHOWERS_[trackIndex].XD = dg_wSHOWERS_[trackIndex].XL + dg_TRACK_mod_[index].U * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].YD = dg_wSHOWERS_[trackIndex].YL + dg_TRACK_mod_[index].V * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].ZD = dg_wSHOWERS_[trackIndex].ZL + dg_TRACK_mod_[index].W * dg_wSHOWERS_[trackIndex].DSEFR;
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XD, dg_wSHOWERS_[trackIndex].YD, dg_wSHOWERS_[trackIndex].ZD, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					else
					{
						dg_wSHOWERS_[trackIndex].DEP = -dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_wSHOWERS_[trackIndex].DS - dg_wSHOWERS_[trackIndex].DSEF) * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XL, dg_wSHOWERS_[trackIndex].YL, dg_wSHOWERS_[trackIndex].ZL, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
					//dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] = dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] + dg_wSHOWERS_[trackIndex].DEP;
				}
				// Verifique se a partícula está fora do invólucro.
				if (dg_TRACK_mod_[index].MAT == 0)
				{ // A partícula está fora do recinto.
					if (dg_TRACK_mod_[index].W > 0.0e0)
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 1; // Marca partículas ascendentes emergentes.
					}
					else
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
					}
					dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
					dg_TRACK_mod_[index].STEP = 9; //g_showers_step3_G/E/P
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
					}
					else
					{
						dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
					}
					// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //Saida
					return;
				}

				// Detectores de Impacto

				dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
				if (dg_wSHOWERS_[trackIndex].IDET != 0)
				{
					if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
					{
						// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
						/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
							NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
							wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
							dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
							dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

						d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

						if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
						{
							atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
							//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
							dg_wSHOWERS_[trackIndex].IEXIT = 3;
							dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
							dg_TRACK_mod_[index].STEP = 9; //g_showers_step6
							if (d_cudaCapability >= 6)
							{
								atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
							}
							else
							{
								dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
							}
							// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
							return;
						}
					}
				}
				//dg_wSHOWERS_[index].CROSS = true;
				//d_showers_step3_E(size);
				// if (dg_TRACK_mod_[index].IEXIT == 0)
				// goto L103;
			}
	//	}
	}
}

__global__ void g_showers_step4_P(int size)
{ // terceiro passo

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{

			// Comprimento do caminho livre para o próximo evento de interação.
			// L103:;
			dg_wSHOWERS_[trackIndex].IBODYL = dg_TRACK_mod_[index].IBODY;
			dg_wSHOWERS_[trackIndex].MATL = dg_TRACK_mod_[index].MAT;
			dg_wSHOWERS_[trackIndex].XL = dg_TRACK_mod_[index].X;
			dg_wSHOWERS_[trackIndex].YL = dg_TRACK_mod_[index].Y;
			dg_wSHOWERS_[trackIndex].ZL = dg_TRACK_mod_[index].Z;

			/*if ((CFORCI_.LFORCE[KPAR-1][IBODY-1]) && ((WGHT > WLOW[KPAR-1][IBODY-1]) && (WGHT < WHIG[KPAR-1][IBODY-1]))){
				//JUMPF(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  //Força de Interação
				//dg_wSHOWERS_[trackIndex].LINTF=true;
				//Esse trecho nao sera implementado pois a versão nao irá tratar redução de variancia.
			}else`{
				jumtp2(DSMAX(IBODY),dg_wSHOWERS_[trackIndex].DS)  ! Analogue simulation.
				dg_wSHOWERS_[trackIndex].LINTF=false
			}*/

			// d_jump2_G(dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DS);

			// d_step2_(dg_wSHOWERS_[trackIndex].DS, dg_wSHOWERS_[trackIndex].DSEF, dg_wSHOWERS_[trackIndex].NCROSS); // Determina a posição final do passo.

			// volta aqui
			dg_wSHOWERS_[trackIndex].LINTF = false;

			// Distribuição de energia de fluência.
			dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1];
			if (dg_wSHOWERS_[trackIndex].IDET != 0)
			{
				if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 2)
				{
					if (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1)
					{
						if (dg_TRACK_mod_[index].KPAR == 2)
						{
							d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
						}
						else
						{
							dg_wSHOWERS_[trackIndex].DECSD = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
							if (dg_wSHOWERS_[trackIndex].DECSD > 1.0e-12)
							{																																										   // A distribuição pode se estender
								d_fimdes2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_PENELOPE_mod_.E0STEP[trackIndex], dg_wSHOWERS_[trackIndex].DECSD, dg_wSHOWERS_[trackIndex].DSEF); // abaixo do EABS.
							}
							else
							{
								d_fimdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET, dg_wSHOWERS_[trackIndex].DSEF);
							}
						}
					}
				}
			}
			// dg_wSHOWERS_[trackIndex].CROSS = false;
			//  A partícula cruzou uma interface.

			if (dg_wSHOWERS_[trackIndex].NCROSS > 0)
			{
				// Correção da perda de energia suave (CSDA).
				if (dg_TRACK_mod_[index].KPAR != 2)
				{
					dg_TRACK_mod_[index].E = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF;
					if (dg_CJUMP1_[trackIndex].MHINGE == 0)
					{
						dg_wSHOWERS_[trackIndex].DEP = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_wSHOWERS_[trackIndex].DSEF * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							dg_wSHOWERS_[trackIndex].DSEFR = d_rand2_(8.0e0) * dg_wSHOWERS_[trackIndex].DSEF;
							dg_wSHOWERS_[trackIndex].XD = dg_wSHOWERS_[trackIndex].XL + dg_TRACK_mod_[index].U * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].YD = dg_wSHOWERS_[trackIndex].YL + dg_TRACK_mod_[index].V * dg_wSHOWERS_[trackIndex].DSEFR;
							dg_wSHOWERS_[trackIndex].ZD = dg_wSHOWERS_[trackIndex].ZL + dg_TRACK_mod_[index].W * dg_wSHOWERS_[trackIndex].DSEFR;
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XD, dg_wSHOWERS_[trackIndex].YD, dg_wSHOWERS_[trackIndex].ZD, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					else
					{
						dg_wSHOWERS_[trackIndex].DEP = -dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_wSHOWERS_[trackIndex].DS - dg_wSHOWERS_[trackIndex].DSEF) * dg_TRACK_mod_[index].WGHT;
						if (dg_CNT6_.LDOSEM)
						{
							d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_wSHOWERS_[trackIndex].XL, dg_wSHOWERS_[trackIndex].YL, dg_wSHOWERS_[trackIndex].ZL, dg_wSHOWERS_[trackIndex].MATL, dg_TRACK_mod_[index].N);
						}
					}
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
					//dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] = dg_CNT1_.DEBO[dg_wSHOWERS_[trackIndex].IBODYL - 1] + dg_wSHOWERS_[trackIndex].DEP;
				}
				// Verifique se a partícula está fora do invólucro.
				if (dg_TRACK_mod_[index].MAT == 0)
				{ // A partícula está fora do recinto.
					if (dg_TRACK_mod_[index].W > 0.0e0)
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 1; // Marca partículas ascendentes emergentes.
					}
					else
					{
						dg_wSHOWERS_[trackIndex].IEXIT = 2; // Rotula partículas descendentes emergentes.
					}
					dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
					}
					else
					{
						dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
					}
					// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //Saida
					return;
				}

				// Detectores de Impacto

				dg_wSHOWERS_[trackIndex].IDET = dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1];
				if (dg_wSHOWERS_[trackIndex].IDET != 0)
				{
					if ((dg_PENGEOM_mod_.KDET[dg_wSHOWERS_[trackIndex].IBODYL - 1] != dg_wSHOWERS_[trackIndex].IDET) && (dg_CNT4_.KKDI[dg_TRACK_mod_[index].KPAR - 1][dg_wSHOWERS_[trackIndex].IDET - 1] == 1))
					{
						// Esse trecho faz gravação no arquivo Phase-Space que não será implementado nesta versao
						/*if (dg_CNT4_.IPSF[dg_wSHOWERS_[trackIndex].IDET-1] == 1){
							NSHJ=dg_CNTRL_.SHN - dg_CNT4_.RLAST;
							wrpsf2_(dg_CNT4_.IPSFO,NSHJ,0);
							dg_CNT4_.RWRITE=dg_CNT4_.RWRITE+1.0e0;
							dg_CNT4_.RLAST=dg_CNTRL_.SHN;*/

						d_simdet2_(dg_TRACK_mod_[index].N, dg_wSHOWERS_[trackIndex].IDET);

						if (dg_CNT4_.IDCUT[dg_wSHOWERS_[trackIndex].IDET - 1] == 0)
						{
							atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT);
							//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
							dg_wSHOWERS_[trackIndex].IEXIT = 3;
							dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
							if (d_cudaCapability >= 6)
							{
								atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
							}
							else
							{
								dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
							}
							// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT);
							return;
						}
					}
				}
				dg_wSHOWERS_[trackIndex].CROSS = true;
				d_showers_step3_P(size);
				// if (dg_TRACK_mod_[index].IEXIT == 0)
				// goto L103;
			}
		}
	}
}

__global__ void g_showers_step5_G(int size)
{

	// Simulação dos eventos de interação da particula com a materia.
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// if ((index >= 20) && (index <= 21))
	// printf("IEXIT %d\n", dg_TRACK_mod_[index].IEXIT);

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 8))
	{
		//if ((dg_wSHOWERS_[trackIndex].CROSS == false) && (dg_TRACK_mod_[index].IEXIT == 0))
		//{
			dg_TRACK_mod_[index].STEP = 4;
			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + dg_TRACK_mod_[index].E;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{																		// Aniquilação de positron
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]); // Quando Absorvida
					dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + d_TREV;
				}
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DE * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;

			if (dg_CNT6_.LDOSEM)
				d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{										////A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				dg_TRACK_mod_[index].STEP = 9; //g_shower_step6
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}
			

			// goto L103;

			// A simulação da particula termina aqui.
		//}
		//dg_wSHOWERS_[trackIndex].CROSS = false;
	}
}

__global__ void g_showers_step5_E(int size)
{

	// Simulação dos eventos de interação da particula com a materia.
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// if ((index >= 20) && (index <= 21))
	//	printf("IEXIT %d\n", dg_TRACK_mod_[index].IEXIT);

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 8))
	{
		// if ((dg_TRACK_mod_[index].IEXIT == 0))
		// printf("Step_5_E NCROSS %d, IEXET %d  Energia %f \n", dg_wSHOWERS_[trackIndex].CROSS, dg_TRACK_mod_[index].IEXIT, dg_TRACK_mod_[index].E);
		//if ((dg_wSHOWERS_[trackIndex].CROSS == false) && (dg_TRACK_mod_[index].IEXIT == 0))
		//{
			dg_TRACK_mod_[index].STEP = 4;
			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + dg_TRACK_mod_[index].E;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{																		// Aniquilação de positron
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]); // Quando Absorvida
					dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + d_TREV;
				}
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DE * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;

			if (dg_CNT6_.LDOSEM)
				d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{										////A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				dg_TRACK_mod_[index].STEP = 9;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			// goto L103;

			// A simulação da particula termina aqui.
		//}
		//dg_wSHOWERS_[trackIndex].CROSS = false;
	}
}

__global__ void g_showers_step5_P(int size)
{

	// Simulação dos eventos de interação da particula com a materia.
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// if ((index >= 20) && (index <= 21))
	// printf("IEXIT %d\n", dg_TRACK_mod_[index].IEXIT);

	if (index < size)
	{
		if ((dg_wSHOWERS_[trackIndex].CROSS == false) && (dg_TRACK_mod_[index].IEXIT == 0))
		{

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{ // A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + dg_TRACK_mod_[index].E;
				if ((dg_TRACK_mod_[index].KPAR == 3) && (dg_TRACK_mod_[index].E > 1.0e-6))
				{																		// Aniquilação de positron
					d_panar2_(dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][2 - 1]); // Quando Absorvida
					dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE + d_TREV;
				}
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			dg_wSHOWERS_[trackIndex].DEP = dg_wSHOWERS_[trackIndex].DE * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], dg_wSHOWERS_[trackIndex].DEP);
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] + dg_wSHOWERS_[trackIndex].DEP;

			if (dg_CNT6_.LDOSEM)
				d_sdose2_(dg_wSHOWERS_[trackIndex].DEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);

			if (dg_TRACK_mod_[index].E < dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
			{										////A partícula foi absorvida.
				dg_wSHOWERS_[trackIndex].IEXIT = 3; // Marca partículas absorvidas.
				dg_TRACK_mod_[index].IEXIT = dg_wSHOWERS_[trackIndex].IEXIT;
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
				}
				else
				{
					dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
				}
				// showers_step3(size, dg_wSHOWERS_[trackIndex].IEXIT); //saida
				return;
			}

			// goto L103;

			// A simulação da particula termina aqui.
		}
	}
}

__global__ void g_showers_step6(int size)
{ // contadores de particulas

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	if ((index < size) && (dg_TRACK_mod_[index].STEP == 9))
	{ 
		// Incrementar contadores de partículas.
		if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
		{
			dg_CNT0_.DPRIM[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] = dg_CNT0_.DPRIM[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] + dg_TRACK_mod_[index].WGHT;
			/*if (*CSOUR0_.LPSF){  não sera implementado redução de variancia como divisão de particulas ou roleta russa
				if (NSPL1 > 1)
					dg_CNT0_.DPRIM[indexx][IEXIT-1]=dg_CNT0_.DPRIM[indexx][IEXIT-1]+dg_TRACK_mod_[indexx].WGHT*(NSPL1-1);
			}
			*/

			if (dg_TRACK_mod_[index].IEXIT < 3)
			{
				dg_CNT0_.DAVW[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] = dg_CNT0_.DAVW[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] + dg_TRACK_mod_[index].W * dg_TRACK_mod_[index].WGHT;
				dg_CNT0_.DAVA[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] = dg_CNT0_.DAVA[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] + acos(dg_TRACK_mod_[index].W) * dg_TRACK_mod_[index].WGHT;
				dg_CNT0_.DAVE[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] = dg_CNT0_.DAVE[trackIndex][dg_TRACK_mod_[index].IEXIT - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			}
		}
		else
		{
			dg_CNT0_.DSEC[trackIndex][dg_TRACK_mod_[index].IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT0_.DSEC[trackIndex][dg_TRACK_mod_[index].IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
		}

		if (dg_TRACK_mod_[index].IEXIT < 3)
		{
			d_tenang2_(dg_TRACK_mod_[index].IEXIT, dg_TRACK_mod_[index].N); // Energia e distribuições angulares.
		}

		dg_TRACK_mod_[index].STEP = 10; //g_showers_step7
	}
	
}

__global__ void g_showers_sec_G(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;
//	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		// int KEn;
		// double DEP, WS, US, VS, SDTS, DF;
		double DEP;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*	if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
			{ // Fonte da particula primaria
				KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
				dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
				dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
				// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
				//	page02_();
				// goto L302; //A energia não é removida do local.
				// printf("passou aqui e chamou step1\n");
				d_showers_step2_G(size); // 302
				return;
			}*/
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			/*if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}*/
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			}
			else
			{
				dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			}
			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		dg_TRACK_mod_[index].INDEX = index;
		dg_TRACK_mod_[index].IEXIT = 0;

		/*else
		{
			return; // L202;
		}
		// goto L102;
		// d_showers_step31_G(size);*/
	}
}

__global__ void g_showers_sec_E(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 11))
	{
//		int trackIndex = dg_TRACK_mod_[index].INDEX;

		// int KEn;
		// double DEP, WS, US, VS, SDTS, DF;
		double DEP;

		//if (index == 0)
		//printf("Energia %lf Step %d\n", dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].STEP);
		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*	if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
			{ // Fonte da particula primaria
				KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
				dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
				dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
				// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
				//	page02_();
				// goto L302; //A energia não é removida do local.
				// printf("passou aqui e chamou step1\n");
				d_showers_step2_G(size); // 302
				return;
			}*/
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			/*if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}*/
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			}
			else
			{
				dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			}
			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
			dg_TRACK_mod_[index].INDEX = index;
			dg_TRACK_mod_[index].IEXIT = 0;
			dg_wSHOWERS_[index].CROSS = false;
			dg_TRACK_mod_[index].STEP = 3;
		}
		else
		{
			dg_TRACK_mod_[index].INDEX = index;

			dg_wSHOWERS_[index].CROSS = false;
			dg_TRACK_mod_[index].IEXIT = 3;
			dg_TRACK_mod_[index].STEP = 20;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nFINISH, -1);
			}
			else
			{
				dg_nTRACKS_.nFINISH = dg_nTRACKS_.nFINISH - 1;
			}
		}

		 //g_showers_step3_E

		/*else
		{
			return; // L202;
		}
		// goto L102;
		// d_showers_step31_G(size);*/
	}
}

__global__ void g_showers_sec_P(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	//int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		// int KEn;
		// double DEP, WS, US, VS, SDTS, DF;
		double DEP;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*	if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
			{ // Fonte da particula primaria
				KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
				dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
				dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
				// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
				//	page02_();
				// goto L302; //A energia não é removida do local.
				// printf("passou aqui e chamou step1\n");
				d_showers_step2_G(size); // 302
				return;
			}*/
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			/*if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}*/
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			if (d_cudaCapability >= 6)
			{
					atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], - DEP);
			}
			else
			{
				dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			}
			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		dg_TRACK_mod_[index].INDEX = index;
		dg_TRACK_mod_[index].IEXIT = 0;

		/*else
		{
			return; // L202;
		}
		// goto L102;
		// d_showers_step31_G(size);*/
	}
}

__global__ void g_showers_step7_G(int size)
{ // particula secundaria fotons

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		int KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
		{ // Fonte da particula primaria
			KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
			dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
			dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
			// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			// goto L302; //A energia não é removida do local.
			// printf("passou aqui e chamou step1\n");
			d_showers_step2_G(size); // 302
			return;
		}
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}
		// goto L102;
		// d_showers_step31_G(size);
	}
}

__global__ void g_showers_step7_E(int size)
{ // particula secundaria fotons

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		int KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
		{ // Fonte da particula primaria
			KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
			dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
			dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
			// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			// goto L302; //A energia não é removida do local.
			// printf("passou aqui e chamou step1\n");
			// showers_step2_E(size);//302
			return;
		}
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}
		// goto L102;
		//	d_showers_step31_E(size);
	}
}

__global__ void g_showers_step7_P(int size)
{ // particula secundaria fotons

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		int KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		if (dg_TRACK_mod_[index].ILB[1 - 1] == 1)
		{ // Fonte da particula primaria
			KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
			dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
			dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
			// if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			// goto L302; //A energia não é removida do local.
			// printf("passou aqui e chamou step1\n");
			d_showers_step2_G(size); // 302
			return;
		}
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}
		// goto L102;
		// d_showers_step31_P(size);
	}
}

__global__ void g_showers_step8(int size)
{ // contadores de deposito de dose

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	int IDET;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 10))
	{

		/*Energias depositadas em diferentes corpos e detectores.
			Calculando os espectros dos detectores de deposição de energia.*/

		if (dg_CNT5_.NED > 0)
		{
			for (int KD = 1; KD <= dg_CNT5_.NED; KD++)
			{
				dg_CNT5_.DEDE[trackIndex][KD - 1] = 0.0e0;
			}
			for (int KB = 1; KB <= dg_PENGEOM_mod_.NBODY; KB++)
			{
				IDET = dg_CNT5_.KBDE[KB - 1];
				if (IDET != 0)
					dg_CNT5_.DEDE[trackIndex][IDET - 1] = dg_CNT5_.DEDE[trackIndex][IDET - 1] + dg_CNT1_.DEBO[KB - 1];
			}

			for (int IDET = 1; IDET <= dg_CNT5_.NED; IDET++)
			{
				d_sendet2_(dg_CNT5_.DEDE[trackIndex][IDET - 1], IDET);
			}
		}
		//__syncthreads();

		 if (index == 0){
		for (int KB = 1; KB <= dg_PENGEOM_mod_.NBODY; KB++)
		{
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_CNT1_.TDEBO[KB - 1], dg_CNT1_.DEBO[KB - 1]);
				atomicAdd2(&dg_CNT1_.TDEBO2[KB - 1], pow(dg_CNT1_.DEBO[KB - 1], 2));
			}
			else
			{
				dg_CNT1_.TDEBO[KB - 1] = dg_CNT1_.TDEBO[KB - 1] + dg_CNT1_.DEBO[KB - 1];
				dg_CNT1_.TDEBO2[KB - 1] = dg_CNT1_.TDEBO2[KB - 1] + (dg_CNT1_.DEBO[KB - 1] * dg_CNT1_.DEBO[KB - 1]);
			}
			//	printf("TDBEBO: %f\n", dg_CNT1_.TDEBO[KB-1]);
			//	printf("TDBEBO2: %f\n\n", dg_CNT1_.TDEBO2[KB-1]);

			//}
		}
		 }

		// Contadores de estado final
		for (int I = 1; I <= 3; I++)
		{
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_CNT0_.PRIM[I - 1], dg_CNT0_.DPRIM[trackIndex][I - 1]);
				atomicAdd2(&dg_CNT0_.PRIM2[I - 1], pow(dg_CNT0_.DPRIM[trackIndex][I - 1], 2));
			}
			else
			{
				dg_CNT0_.PRIM[I - 1] = dg_CNT0_.PRIM[I - 1] + dg_CNT0_.DPRIM[trackIndex][I - 1];
				dg_CNT0_.PRIM2[I - 1] = dg_CNT0_.PRIM2[I - 1] + (dg_CNT0_.DPRIM[trackIndex][I - 1] * dg_CNT0_.DPRIM[trackIndex][I - 1]);
			}
			for (int K = 1; K <= 3; K++)
			{
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_CNT0_.SEC[I - 1][K - 1], dg_CNT0_.DSEC[trackIndex][I - 1][K - 1]);
					atomicAdd2(&dg_CNT0_.SEC2[I - 1][K - 1], pow(dg_CNT0_.DSEC[trackIndex][I - 1][K - 1], 2));
				}
				else
				{
					dg_CNT0_.SEC[I - 1][K - 1] = dg_CNT0_.SEC[I - 1][K - 1] + dg_CNT0_.DSEC[trackIndex][I - 1][K - 1];
					dg_CNT0_.SEC2[I - 1][K - 1] = dg_CNT0_.SEC2[I - 1][K - 1] + (dg_CNT0_.DSEC[trackIndex][I - 1][K - 1] * dg_CNT0_.DSEC[trackIndex][I - 1][K - 1]);
				}
			}
		}

		for (int I = 1; I <= 2; I++)
		{
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_CNT0_.AVW[I - 1], dg_CNT0_.DAVW[trackIndex][I - 1]);
				atomicAdd2(&dg_CNT0_.AVW2[I - 1], pow(dg_CNT0_.DAVW[trackIndex][I - 1], 2));
				atomicAdd2(&dg_CNT0_.AVA[I - 1], dg_CNT0_.DAVA[trackIndex][I - 1]);
				atomicAdd2(&dg_CNT0_.AVA2[I - 1], pow(dg_CNT0_.DAVA[trackIndex][I - 1], 2));
				atomicAdd2(&dg_CNT0_.AVE[I - 1], dg_CNT0_.DAVE[trackIndex][I - 1]);
				atomicAdd2(&dg_CNT0_.AVE2[I - 1], pow(dg_CNT0_.DAVE[trackIndex][I - 1], 2));
			}
			else
			{
				dg_CNT0_.AVW[I - 1] = dg_CNT0_.AVW[I - 1] + dg_CNT0_.DAVW[trackIndex][I - 1];
				dg_CNT0_.AVW2[I - 1] = dg_CNT0_.AVW2[I - 1] + (dg_CNT0_.DAVW[trackIndex][I - 1] * dg_CNT0_.DAVW[trackIndex][I - 1]);
				dg_CNT0_.AVA[I - 1] = dg_CNT0_.AVA[I - 1] + dg_CNT0_.DAVA[trackIndex][I - 1];
				dg_CNT0_.AVA2[I - 1] = dg_CNT0_.AVA2[I - 1] + (dg_CNT0_.DAVA[trackIndex][I - 1] * dg_CNT0_.DAVA[trackIndex][I - 1]);
				dg_CNT0_.AVE[I - 1] = dg_CNT0_.AVE[I - 1] + dg_CNT0_.DAVE[trackIndex][I - 1];
				dg_CNT0_.AVE2[I - 1] = dg_CNT0_.AVE2[I - 1] + (dg_CNT0_.DAVE[trackIndex][I - 1] * dg_CNT0_.DAVE[trackIndex][I - 1]);
			}
		}
	}
}

__global__ void d1_jump2_G(double &DSMAX, double &DS)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
}

__global__ void d2_jump2_G(double &DSMAX, double &DS)
{
	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// Fotons

	if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
	{
		dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
		d_gimfp2_();
		dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[1 - 1] + dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	}

	DS = -log(d_rand2_(1.0e0)) / dg_CJUMP0_[trackIndex].ST;
}

__global__ void showers_sec_G(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		/*		 if ((index == 0) || (index == 1) || (index == 4095) || (index == 4096) || (index == 256)){
				printf("index %d, 	energia %lf\n", index, dg_TRACK_mod_[index].E);
				}*/

		/*	if ((index == 0) || (index == 1)){
				printf("SHOWER SEC DEBO[%d]: %f\n", index, dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1]);
			}*/

		//	int  KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*if (dg_TRACK_mod_[index].ILB[1 - 1] == 1) { //Fonte da particula primaria
			KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
			dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
			dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
			//if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			//goto L302; //A energia não é removida do local.
			//printf("passou aqui e chamou step1\n");
			showers_step1(size);//302
			return;
		} */
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}

		// goto L102;
		// printf("passou aqui e chamou step2\n");
		d_showers_step2_G(size);
		// chmar showers_step 4 para contabilizar os depositos de dose apos todas as particulas secundasrias terem sido simuladas.
		// showers_step4(size);
	}
}

__global__ void showers_sec_E(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		/*		 if ((index == 0) || (index == 1) || (index == 4095) || (index == 4096) || (index == 256)){
				printf("index %d, 	energia %lf\n", index, dg_TRACK_mod_[index].E);
				}*/

		/*	if ((index == 0) || (index == 1)){
				printf("SHOWER SEC DEBO[%d]: %f\n", index, dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1]);
			}*/

		// int  KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*	if (dg_TRACK_mod_[index].ILB[1 - 1] == 1) { //Fonte da particula primaria
				KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
				dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
				dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
				//if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
				//	page02_();
				//goto L302; //A energia não é removida do local.
				//printf("passou aqui e chamou step1\n");
				showers_step1(size);//302
				return;
			} */
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}

		// goto L102;
		// printf("passou aqui e chamou step2\n");
		//	showers_step2_E(size);
		// chmar showers_step 4 para contabilizar os depositos de dose apos todas as particulas secundasrias terem sido simuladas.
		// showers_step4(size);
	}
}

__global__ void showers_sec_P(int size)
{

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if (index < size)
	{

		/*		 if ((index == 0) || (index == 1) || (index == 4095) || (index == 4096) || (index == 256)){
				printf("index %d, 	energia %lf\n", index, dg_TRACK_mod_[index].E);
				}*/

		/*	if ((index == 0) || (index == 1)){
				printf("SHOWER SEC DEBO[%d]: %f\n", index, dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1]);
			}*/

		// int  KEn;
		double DEP, WS, US, VS, SDTS, DF;

		// Particulas Secundarias

		// L202:;
		//	d_secpar2_(LEFT);
		// if (LEFT > 0) {

		/*if (dg_TRACK_mod_[index].ILB[1 - 1] == 1) { //Fonte da particula primaria
			KEn = int(dg_TRACK_mod_[index].E * dg_CNT3_.RDSDE + 1.0e0);
			dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
			dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CNT3_.SEDS2[trackIndex][KEn - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_TRACK_mod_[index].WGHT, 2);
			//if (dg_TRACK_mod_[index].LAGE) //nao irá contabilizar a idade da particula
			//	page02_();
			//goto L302; //A energia não é removida do local.
			//printf("passou aqui e chamou step1\n");
			showers_step1(size);//302
			return;
		} */
		if (dg_TRACK_mod_[index].E > dg_CSPGEO_.EABSB[dg_TRACK_mod_[index].IBODY - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			// Divisoes de raios X
			if (dg_TRACK_mod_[index].KPAR == 2)
			{
				if (dg_TRACK_mod_[index].ILB[4 - 1] > 0)
				{ // caracteristica de raio X
					if ((dg_TRACK_mod_[index].ILB[1 - 1] == 2) && (dg_TRACK_mod_[index].ILB[3 - 1] < 9))
					{
						if (dg_CXRSPL_.LXRSPL[dg_TRACK_mod_[index].IBODY - 1])
						{
							// printf("processando secundaria\n");
							dg_TRACK_mod_[index].WGHT = dg_TRACK_mod_[index].WGHT / (dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]);
							dg_CXRSPL_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1];
							dg_CXRSPL_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].ILB[2 - 1];
							dg_CXRSPL_.ILBA[trackIndex][3 - 1] = 9;
							dg_CXRSPL_.ILBA[trackIndex][4 - 1] = dg_TRACK_mod_[index].ILB[4 - 1];
							dg_CXRSPL_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
							for (int I = 2; I <= dg_CXRSPL_.IXRSPL[dg_TRACK_mod_[index].IBODY - 1]; I++)
							{
								WS = -1.0e0 + 2.0e0 * d_rand2_(9.0e0);
								SDTS = sqrt(1.0e0 - WS * WS);
								DF = d_TWOPI * d_rand2_(10.0e0);
								US = cos(DF) * SDTS;
								VS = sin(DF) * SDTS;
								d_stores2_(dg_TRACK_mod_[index].E, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, dg_TRACK_mod_[index].KPAR, dg_CXRSPL_.ILBA[trackIndex], d_wIPOLI);
							}
						}
					}
				}
			}
			// printf("\nprocessando secundaria finaliznado.\n");
			DEP = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
			atomicAdd2(&dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1], -DEP);
			// dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;
			//dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] = dg_CNT1_.DEBO[dg_TRACK_mod_[index].IBODY - 1] - DEP;

			if (dg_CNT6_.LDOSEM)
			{
				double wDEP = -DEP;
				d_sdose2_(wDEP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, dg_TRACK_mod_[index].MAT, dg_TRACK_mod_[index].N);
			}
		}
		else
		{
			return; // L202;
		}

		// goto L102;
		// printf("passou aqui e chamou step2\n");
		//	showers_step2_P(size);
		// chmar showers_step 4 para contabilizar os depositos de dose apos todas as particulas secundasrias terem sido simuladas.
		// showers_step4(size);
	}
}

__global__ void g_step2_(int size)
{
	/*
	Esta sub-rotina lida com a parte geom�trica da simulacao de pista
   ��o A particula come�a no ponto (X, Y, Z) e percorre um comprimento
   DS na dire��o (U, V, W) dentro do material para onde se move. Quando
   a trilha deixa o material inicial, a part�cula � parada apenas
   depois de entrar no pr�ximo corpo material (regi�es vazias com MAT = 0 s�o
   cruzado automaticamente). Al�m disso, quando a part�cula chega de
   uma regi�o vazia, ele � interrompido logo ap�s entrar no primeiro material
   corpo.

   Valores de entrada (m�dulo TRACK_mod):
	  X, Y, Z ... coordenadas do ponto inicial,
	  U, V, W ... cossenos de dire��o do deslocamento,
	  IBODY ..... corpo onde o ponto inicial est� localizado,
	  MAT ....... material no corpo IBODY.
   NB: Quando uma trilha de part�culas � iniciada, as vari�veis ??IBODY e MAT
   deve ser definido chamando a sub-rotina LOCATE.

   Argumento de entrada:
	  DS ........ comprimento do caminho a percorrer.

   Argumentos de sa�da:
	  DSEF ....... percorreu o comprimento do caminho antes de deixar o inicial
				  material ou completar o salto (menos do que DS se o
				  a trilha atravessa uma interface),
	  NCROSS .... = 0 se toda a etapa estiver contida no inicial
					material,
				  .gt.0 se a part�cula cruzou uma interface, ou seja,
					se entrou em um novo material.

   Valores de sa�da (m�dulo TRACK_mod):
	  X, Y, Z ... coordenadas da posi��o final,
	  IBODY ..... corpo onde o ponto final est� localizado,
	  MAT ....... material em IBODY. O valor MAT = 0 indica que o
				  part�cula escapou do sistema.

   Valores de sa�da (m�dulo dg_PENGEOM_mod):
	  DSTOT ..... comprimento do caminho percorrido, incluindo segmentos de caminho no vazio
				  volumes.
	  KSLAST .... quando NCROSS.ne.0, o valor de sa�da de KSLAST � o
		 r�tulo da �ltima superf�cie cruzada pela part�cula antes
		 entrar em um corpo material. KSLAST � usado para renderiza��o em 3D.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int trackIndex = dg_TRACK_mod_[index].INDEX;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 5))
	{

		if (dg_TRACK_mod_[index].IEXIT == 0)
		{

			dg_TRACK_mod_[index].STEP = 6; //g_showers_step4_G/E/P
			dg_wSHOWERS_[trackIndex].DSEF = 0.0e0;
			dg_PENGEOM_mod_.DSTOT[trackIndex] = 0.0e0;
			dg_wSHOWERS_[trackIndex].NCROSS = 0;
			dg_PENGEOM_mod_.KSLAST[trackIndex] = 0;
			double DSRES;
			int KB1;
			// double d_S[trackIndex][NS2M];
			// int d_ISNS2M];

			// double *d_S[trackIndex];
			// int *d_IS[trackIndex];

			//__host__​__device__​cudaError_t cudaMalloc ( void** devPtr, size_t size )
			// gpuErrchk(cudaMalloc(&d_S[trackIndex], sizeof(double)*NS2M));
			// gpuErrchk(cudaMalloc(&d_IS[trackIndex], sizeof(int)*NS2M));

			//double *d_S[trackIndex][trackIndex] = (double *)malloc(dg_QSURF_.NSURF * 2 * sizeof(double));
			//int *d_IS[trackIndex] = (int *)malloc(dg_QSURF_.NSURF * 2 * sizeof(int));

			/*	if (index == 0){

					printf("\n step2 aqui %lf\n", dg_TRACK_mod_[0].E);

				}*/

			int NST;
			double DSP;
			int KS1;
			int KF;
			int IERR;
			int IBODYL;
			int NERR;
			int KS;
			int KFLO;
			double SW;
			int MATL;
			double A, B, C;

			int NSC = 0; // Número de cruzamentos da superfície à frente da partícula.
			int NSCT;
			int MAT0;

			//	printf("1\n");

			//	printf("LVERB %d\n", dg_PENGEOM_mod_.LVERB);
			//	printf(dg_PENGEOM_mod_.LVERB ? "true\n" : "false\n")
			for (int I = 1; I <= dg_QSURF_.NSURF; I++)
			{
				dg_QTREE_.KSP[trackIndex][I - 1] = 0; // Ponteiros laterais das superfícies avaliadas.
			}

			//	printf("2\n");

			MAT0 = dg_TRACK_mod_[index].MAT; // Material Inicial

			if (dg_TRACK_mod_[index].MAT == 0)
			{
				DSRES = 1.0e35; // No vácuo, as partículas voam livremente.
			}
			else
			{
				DSRES = dg_wSHOWERS_[trackIndex].DS; // comprimento do camimho residual
			}

			// A partícula entra de fora do recinto.

			if (dg_TRACK_mod_[index].IBODY > dg_QTREE_.NBODYS)
			{
				KB1 = dg_QTREE_.NBODYS;
				d_stepsi2_(KB1, NSC);
				if (NSC == 0)
					goto L300;
				NSCT = NSC;
				NST = dg_QTREE_.KSURF[NXG - 1][KB1 - 1];

				for (int KI = NSCT; KI >= 1; KI--)
				{
					// a particula atravessa uma superficie
					dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][KI - 1];
					if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
						dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
					else
						dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;

					DSP = d_S[trackIndex][KI - 1];
					dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSP;
					dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
					dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
					dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
					dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
					NSC = NSC - 1;

					if (NSC > 0)
					{
						for (int I = 1; I <= NSC; I++)
						{
							d_S[trackIndex][I - 1] = d_S[trackIndex][I - 1] - DSP;
						}
					}

					for (int KSS = 1; KSS <= NST; KSS++)
					{
						KS1 = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
						KF = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
						if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS1 - 1] != KF))
							goto L101;
					}
					// A partícula entra no invólucro.
				L100:;
					d_steplb2_(KB1, IERR);
					// A partícula entra em um submódulo.
					if (IERR == -1)
					{
						KB1 = dg_TRACK_mod_[index].IBODY;
						d_stepsi2_(KB1, NSC);
						goto L100;
					}
					else
					{
						// A particula entrou em um corpo material
						if (dg_TRACK_mod_[index].MAT != 0)
						{
							dg_wSHOWERS_[trackIndex].NCROSS = 1;
							//free(d_S[trackIndex]);
							//free(d_IS[trackIndex]);
							return;
						}
						else
						{
							KB1 = dg_TRACK_mod_[index].IBODY;
							d_stepsi2_(KB1,  NSC);
							goto L200;
						}
					}
					// Neste ponto, o programa saiu do ciclo DO.
				L101:;
				}
				goto L300;
			}

			// Cruzamentos de superfície.
			IBODYL = dg_TRACK_mod_[index].IBODY;
			if (dg_PENGEOM_mod_.LVERB)
				NERR = 0;
		L102:;
			KB1 = dg_TRACK_mod_[index].IBODY;
			d_stepsi2_(KB1,  NSC);
			d_steplb2_(KB1, IERR);

			// Evidência de erros de arredondamento.
			if (IERR != 0)
			{
				if (NSC > 0)
				{
					// Quando uma superfície está muito próxima, movemos a partícula além dela.
					if (d_S[trackIndex][NSC - 1] < 1.0e-10)
					{
						dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][NSC - 1];
						if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
						{
							dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
						}
						else
						{
							dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;
						}

						DSP = d_S[trackIndex][NSC - 1];
						dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
						dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
						dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
						if (dg_TRACK_mod_[index].MAT == MAT0)
						{
							dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSP;
							DSRES = DSRES - DSP;
						}

						dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
						NSC = NSC - 1;

						if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
							goto L102;
					}
				}
				if (dg_PENGEOM_mod_.LVERB)
				{

					NERR = NERR + 1;
					if ((dg_QTREE_.NWARN[trackIndex] < 100) && (dg_TRACK_mod_[index].MAT != 0))
					{

						printf("WARNING, STEP: Accidental undershot or r");

						for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB1 - 1]; KSS++)
						{
							KS = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
							KFLO = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
							if (KFLO < 3)
							{
								for (int KI = NSC; KI >= 1; KI--)
								{
									if (KS == d_IS[trackIndex][KI - 1])
									{
										SW = d_S[trackIndex][KI - 1];
										goto L103;
									}
								}

								SW = 0.0e0;
							L103:;
								d_fsurf2_(KS, A, B, C);
								if (KFLO == dg_QTREE_.KSP[trackIndex][KS - 1])
								{
									printf("KS, KFLO, KSP, SW %lf", SW);
								}
								else
								{
									printf("KS, KFLO, KSP, SW %lf", SW);
									dg_PENGEOM_mod_.KSLAST[trackIndex] = KS;
								}
							}
						}
						dg_QTREE_.NWARN[trackIndex] = dg_QTREE_.NWARN[trackIndex] + 1;
					}
				}
				if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
					goto L102;
			}

			if ((dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1] != dg_PENGEOM_mod_.KDET[IBODYL - 1]) || (dg_TRACK_mod_[index].MAT != MAT0))
			{
				dg_wSHOWERS_[trackIndex].NCROSS = 1;
				dg_wSHOWERS_[trackIndex].DSEF = 0.0e0;
				//free(d_S[trackIndex]);
				//free(d_IS[trackIndex]);
				return;
			}

			// A particula permanece no mesmo material

			if ((dg_TRACK_mod_[index].MAT != 0) && (DSRES < d_S[trackIndex][NSC - 1]))
			{
				if (dg_TRACK_mod_[index].MAT == MAT0)
					dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSRES;
				dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
				dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
				dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
				dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
				//free(d_S[trackIndex]);
				//free(d_IS[trackIndex]);
				return;
			}

			// Nova posição

		L200:;
			if (NSC == 0)
			{
				if (dg_TRACK_mod_[index].MAT == MAT0)
					dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSRES;
				dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
				dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
				dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
				dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
				//free(d_S[trackIndex]);
				//free(d_IS[trackIndex]);
				return;
			}
			NSCT = NSC;
			MATL = dg_TRACK_mod_[index].MAT;
			IBODYL = dg_TRACK_mod_[index].IBODY;
			for (int KI = NSCT; KI >= 1; KI--)
			{
				// O passo termina dentro do corpo
				if (DSRES < d_S[trackIndex][KI - 1])
				{
					if (dg_TRACK_mod_[index].MAT == MAT0)
						dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSRES;
					dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
					dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
					dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
					dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
					//free(d_S[trackIndex]);
					//free(d_IS[trackIndex]);
					return;
				}

				// A particula atravessa uma superfice
				dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][KI - 1];
				if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
					dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
				else
					dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;

				DSP = d_S[trackIndex][KI - 1];
				dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
				dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
				dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
				if (dg_TRACK_mod_[index].MAT == MAT0)
				{
					dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSP;
					DSRES = DSRES - DSP;
				}
				dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
				NSC = NSC - 1;
				if (NSC > 0)
				{
					for (int I = 1; I <= NSC; I++)
					{
						d_S[trackIndex][I - 1] = d_S[trackIndex][I - 1] - DSP;
					}
				}
				d_steplb2_(KB1, IERR);

			L201:;
				KB1 = dg_TRACK_mod_[index].IBODY;
				if (IERR == -1)
				{
					// A particula entrou em um submodulo
					d_stepsi2_(KB1, NSC);
					d_steplb2_(KB1, IERR);
					goto L201;
				}
				else if (IERR == 1)
				{
					// A partícula deixa o corpo ou módulo.
					if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
					{
						d_stepsi2_(KB1, NSC);
						d_steplb2_(KB1, IERR);
						goto L201;
					}
					else
					{
						// A partícula sai do recinto.
						if (dg_TRACK_mod_[index].MAT != MATL)
							dg_wSHOWERS_[trackIndex].NCROSS = dg_wSHOWERS_[trackIndex].NCROSS + 1;
						goto L300;
					}
				}

				// A partícula continua voando quando entra em uma região vazia
				if (dg_TRACK_mod_[index].MAT == 0)
				{
					if (MATL == MAT0)
						dg_wSHOWERS_[trackIndex].NCROSS = dg_wSHOWERS_[trackIndex].NCROSS + 1;
					MATL = 0;
					DSRES = 1.0e35;
					goto L202;
					
				}
				// A partícula continua voando quando entra em um novo corpo do
					// mesmo material que não faz parte de um detector diferente ...
				else if (dg_TRACK_mod_[index].MAT == MATL)
				{
					if (dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1] == dg_PENGEOM_mod_.KDET[IBODYL - 1])
					{
						goto L202;
					}
					else
					{
						dg_wSHOWERS_[trackIndex].NCROSS = dg_wSHOWERS_[trackIndex].NCROSS + 1;
						//free(d_S[trackIndex]);
						//free(d_IS[trackIndex]);
						return;
					}
					//.. e para quando penetra um novo corpo material ou um Detector
				}
				else
				{
					dg_wSHOWERS_[trackIndex].NCROSS = dg_wSHOWERS_[trackIndex].NCROSS + 1;
					//free(d_S[trackIndex]);
					//free(d_IS[trackIndex]);
					return;
				}
			L202:;
				d_stepsi2_(KB1, NSC);
				goto L200;
				// Neste ponto, o programa saiu do ciclo DO.
				// L203:; indica o final do loop
			}
			// A particula sai do recinto.
		L300:;

			DSP = 1.0e36;
			dg_TRACK_mod_[index].IBODY = dg_QTREE_.NBODYS + 1;
			dg_TRACK_mod_[index].MAT = 0;
			if (dg_TRACK_mod_[index].MAT == MAT0)
				dg_wSHOWERS_[trackIndex].DSEF = dg_wSHOWERS_[trackIndex].DSEF + DSP;
			dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
			dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
			dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
			dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
			//free(d_S[trackIndex]);
			//free(d_IS[trackIndex]);
		}
	}
}

__device__ void d_step2_(double &DS, double &DSEF, int &NCROSS)
{
	/*
	Esta sub-rotina lida com a parte geom�trica da simulacao de pista
   ��o A particula come�a no ponto (X, Y, Z) e percorre um comprimento
   DS na dire��o (U, V, W) dentro do material para onde se move. Quando
   a trilha deixa o material inicial, a part�cula � parada apenas
   depois de entrar no pr�ximo corpo material (regi�es vazias com MAT = 0 s�o
   cruzado automaticamente). Al�m disso, quando a part�cula chega de
   uma regi�o vazia, ele � interrompido logo ap�s entrar no primeiro material
   corpo.

   Valores de entrada (m�dulo TRACK_mod):
	  X, Y, Z ... coordenadas do ponto inicial,
	  U, V, W ... cossenos de dire��o do deslocamento,
	  IBODY ..... corpo onde o ponto inicial est� localizado,
	  MAT ....... material no corpo IBODY.
   NB: Quando uma trilha de part�culas � iniciada, as vari�veis ??IBODY e MAT
   deve ser definido chamando a sub-rotina LOCATE.

   Argumento de entrada:
	  DS ........ comprimento do caminho a percorrer.

   Argumentos de sa�da:
	  DSEF ....... percorreu o comprimento do caminho antes de deixar o inicial
				  material ou completar o salto (menos do que DS se o
				  a trilha atravessa uma interface),
	  NCROSS .... = 0 se toda a etapa estiver contida no inicial
					material,
				  .gt.0 se a part�cula cruzou uma interface, ou seja,
					se entrou em um novo material.

   Valores de sa�da (m�dulo TRACK_mod):
	  X, Y, Z ... coordenadas da posi��o final,
	  IBODY ..... corpo onde o ponto final est� localizado,
	  MAT ....... material em IBODY. O valor MAT = 0 indica que o
				  part�cula escapou do sistema.

   Valores de sa�da (m�dulo dg_PENGEOM_mod):
	  DSTOT ..... comprimento do caminho percorrido, incluindo segmentos de caminho no vazio
				  volumes.
	  KSLAST .... quando NCROSS.ne.0, o valor de sa�da de KSLAST � o
		 r�tulo da �ltima superf�cie cruzada pela part�cula antes
		 entrar em um corpo material. KSLAST � usado para renderiza��o em 3D.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	DSEF = 0.0e0;
	dg_PENGEOM_mod_.DSTOT[trackIndex] = 0.0e0;
	NCROSS = 0;
	dg_PENGEOM_mod_.KSLAST[trackIndex] = 0;
	double DSRES;
	int KB1;
	// double d_S[trackIndex][NS2M];
	// int d_ISNS2M];

	// double *d_S[trackIndex];
	// int *d_IS[trackIndex];

	//__host__​__device__​cudaError_t cudaMalloc ( void** devPtr, size_t size )
	// gpuErrchk(cudaMalloc(&d_S[trackIndex], sizeof(double)*NS2M));
	// gpuErrchk(cudaMalloc(&d_IS[trackIndex], sizeof(int)*NS2M));

	//double *d_S[trackIndex] = (double *)malloc(dg_QSURF_.NSURF * 2 * sizeof(double));
	//int *d_IS[trackIndex] = (int *)malloc(dg_QSURF_.NSURF * 2 * sizeof(int));

	/*	if (index == 0){

			printf("\n step2 aqui %lf\n", dg_TRACK_mod_[0].E);

		}*/

	int NST;
	double DSP;
	int KS1;
	int KF;
	int IERR;
	int IBODYL;
	int NERR;
	int KS;
	int KFLO;
	double SW;
	int MATL;
	double A, B, C;

	int NSC = 0; // Número de cruzamentos da superfície à frente da partícula.
	int NSCT;
	int MAT0;

	//	printf("1\n");

	//	printf("LVERB %d\n", dg_PENGEOM_mod_.LVERB);
	//	printf(dg_PENGEOM_mod_.LVERB ? "true\n" : "false\n")
	for (int I = 1; I <= dg_QSURF_.NSURF; I++)
	{
		dg_QTREE_.KSP[trackIndex][I - 1] = 0; // Ponteiros laterais das superfícies avaliadas.
	}

	//	printf("2\n");

	MAT0 = dg_TRACK_mod_[index].MAT; // Material Inicial

	if (dg_TRACK_mod_[index].MAT == 0)
	{
		DSRES = 1.0e35; // No vácuo, as partículas voam livremente.
	}
	else
	{
		DSRES = DS; // comprimento do camimho residual
	}

	//		printf("3\n");
	// A partícula entra de fora do recinto.

	if (dg_TRACK_mod_[index].IBODY > dg_QTREE_.NBODYS)
	{
		KB1 = dg_QTREE_.NBODYS;
		d_stepsi2_(KB1,  NSC);
		// d_stepsi2_(KB1, NSC);
		if (NSC == 0)
			goto L300;
		NSCT = NSC;
		NST = dg_QTREE_.KSURF[NXG - 1][KB1 - 1];

		for (int KI = NSCT; KI >= 1; KI--)
		{
			// a particula atravessa uma superficie
			dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][KI - 1];
			if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
				dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
			else
				dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;

			DSP = d_S[trackIndex][KI - 1];
			DSEF = DSEF + DSP;
			dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
			dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
			dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
			dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
			NSC = NSC - 1;

			if (NSC > 0)
			{
				for (int I = 1; I <= NSC; I++)
				{
					d_S[trackIndex][I - 1] = d_S[trackIndex][I - 1] - DSP;
				}
			}

			for (int KSS = 1; KSS <= NST; KSS++)
			{
				KS1 = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
				KF = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
				if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS1 - 1] != KF))
					goto L101;
			}
			// A partícula entra no invólucro.
		L100:;
			d_steplb2_(KB1, IERR);
			// A partícula entra em um submódulo.
			if (IERR == -1)
			{
				KB1 = dg_TRACK_mod_[index].IBODY;
				d_stepsi2_(KB1,  NSC);
				// d_stepsi2_(KB1, NSC);
				goto L100;
			}
			else
			{
				// A particula entrou em um corpo material
				if (dg_TRACK_mod_[index].MAT != 0)
				{
					NCROSS = 1;
					//free(d_S[trackIndex]);
					//free(d_IS[trackIndex]);
					return;
				}
				else
				{
					KB1 = dg_TRACK_mod_[index].IBODY;
					d_stepsi2_(KB1,  NSC);
					// d_stepsi2_(KB1, NSC);
					goto L200;
				}
			}

			// Neste ponto, o programa saiu do ciclo DO.
		L101:;
		}
		//	printf("6\n");
		goto L300;
	}

	// Cruzamentos de superfície.

	IBODYL = dg_TRACK_mod_[index].IBODY;
	if (dg_PENGEOM_mod_.LVERB)
		NERR = 0;
L102:;
	KB1 = dg_TRACK_mod_[index].IBODY;
	d_stepsi2_(KB1,  NSC);
	// d_stepsi2_(KB1, NSC);
	d_steplb2_(KB1, IERR);

	// Evidência de erros de arredondamento.
	if (IERR != 0)
	{
		if (NSC > 0)
		{
			// Quando uma superfície está muito próxima, movemos a partícula além dela.
			if (d_S[trackIndex][NSC - 1] < 1e-10)
			{
				dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][NSC - 1];
				if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
				{
					dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
				}
				else
				{
					dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;
				}

				DSP = d_S[trackIndex][NSC - 1];
				dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
				dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
				dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
				if (dg_TRACK_mod_[index].MAT == MAT0)
				{
					DSEF = DSEF + DSP;
					DSRES = DSRES - DSP;
				}

				dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
				NSC = NSC - 1;

				if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
					goto L102;
			}
		}
		if (dg_PENGEOM_mod_.LVERB)
		{

			NERR = NERR + 1;
			if ((dg_QTREE_.NWARN[trackIndex] < 100) && (dg_TRACK_mod_[index].MAT != 0))
			{

				printf("WARNING, STEP: Accidental undershot or r");

				for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB1 - 1]; KSS++)
				{
					KS = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
					KFLO = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
					if (KFLO < 3)
					{
						for (int KI = NSC; KI >= 1; KI--)
						{
							if (KS == d_IS[trackIndex][KI - 1])
							{
								SW = d_S[trackIndex][KI - 1];
								goto L103;
							}
						}

						SW = 0.0e0;
					L103:;
						d_fsurf2_(KS, A, B, C);
						if (KFLO == dg_QTREE_.KSP[trackIndex][KS - 1])
						{
							printf("KS, KFLO, KSP, SW %lf", SW);
						}
						else
						{
							printf("KS, KFLO, KSP, SW %lf", SW);
							dg_PENGEOM_mod_.KSLAST[trackIndex] = KS;
						}
					}
				}
				dg_QTREE_.NWARN[trackIndex] = dg_QTREE_.NWARN[trackIndex] + 1;
			}
		}
		if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
			goto L102;
	}

	if ((dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1] != dg_PENGEOM_mod_.KDET[IBODYL - 1]) || (dg_TRACK_mod_[index].MAT != MAT0))
	{
		NCROSS = 1;
		DSEF = 0.0e0;
		//free(d_S[trackIndex]);
		//free(d_IS[trackIndex]);
		return;
	}

	// A particula permanece no mesmo material

	if ((dg_TRACK_mod_[index].MAT != 0) && (DSRES < d_S[trackIndex][NSC - 1]))
	{
		if (dg_TRACK_mod_[index].MAT == MAT0)
			DSEF = DSEF + DSRES;
		dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
		dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
		dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
		dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
		//free(d_S[trackIndex]);
		//free(d_IS[trackIndex]);
		return;
	}

	// Nova posição

L200:;
	if (NSC == 0)
	{
		if (dg_TRACK_mod_[index].MAT == MAT0)
			DSEF = DSEF + DSRES;
		dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
		dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
		dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
		dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
		return;
	}
	NSCT = NSC;
	MATL = dg_TRACK_mod_[index].MAT;
	IBODYL = dg_TRACK_mod_[index].IBODY;
	for (int KI = NSCT; KI >= 1; KI--)
	{
		// A etapa termina dentro do corpo
		if (DSRES < d_S[trackIndex][KI - 1])
		{
			if (dg_TRACK_mod_[index].MAT == MAT0)
				DSEF = DSEF + DSRES;
			dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSRES;
			dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSRES * dg_TRACK_mod_[index].U;
			dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSRES * dg_TRACK_mod_[index].V;
			dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSRES * dg_TRACK_mod_[index].W;
			//free(d_S[trackIndex]);
			//free(d_IS[trackIndex]);
			return;
		}

		// A particula atravessa uma superfice
		dg_PENGEOM_mod_.KSLAST[trackIndex] = d_IS[trackIndex][KI - 1];
		if (dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] == 1)
			dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 2;
		else
			dg_QTREE_.KSP[trackIndex][dg_PENGEOM_mod_.KSLAST[trackIndex] - 1] = 1;

		DSP = d_S[trackIndex][KI - 1];
		dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
		dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
		dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
		if (dg_TRACK_mod_[index].MAT == MAT0)
		{
			DSEF = DSEF + DSP;
			DSRES = DSRES - DSP;
		}
		dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
		NSC = NSC - 1;
		if (NSC > 0)
		{
			for (int I = 1; I <= NSC; I++)
			{
				d_S[trackIndex][I - 1] = d_S[trackIndex][I - 1] - DSP;
			}
		}
		d_steplb2_(KB1, IERR);

	L201:;
		KB1 = dg_TRACK_mod_[index].IBODY;
		if (IERR == -1)
		{
			// A particula entrou em um submodulo
			d_stepsi2_(KB1, NSC);
			// d_stepsi2_(KB1, NSC);
			d_steplb2_(KB1, IERR);
			goto L201;
		}
		else if (IERR == 1)
		{
			// A partícula deixa o corpo ou módulo.
			if (dg_TRACK_mod_[index].IBODY <= dg_QTREE_.NBODYS)
			{
				d_stepsi2_(KB1, NSC);
				// d_stepsi2_(KB1, NSC);
				d_steplb2_(KB1, IERR);
				goto L201;
			}
			else
			{
				// A partícula sai do recinto.
				if (dg_TRACK_mod_[index].MAT != MATL)
					NCROSS = NCROSS + 1;
				goto L300;
			}
		}

		// A partícula continua voando quando entra em uma região vazia
		if (dg_TRACK_mod_[index].MAT == 0)
		{
			if (MATL == MAT0)
				NCROSS = NCROSS + 1;
			MATL = 0;
			DSRES = 1.0e35;
			goto L202;
			// A partícula continua voando quando entra em um novo corpo do
			// mesmo material que não faz parte de um detector diferente ...
		}
		else if (dg_TRACK_mod_[index].MAT == MATL)
		{
			if (dg_PENGEOM_mod_.KDET[dg_TRACK_mod_[index].IBODY - 1] == dg_PENGEOM_mod_.KDET[IBODYL - 1])
			{
				goto L202;
			}
			else
			{
				NCROSS = NCROSS + 1;
				//free(d_S[trackIndex]);
				//free(d_IS[trackIndex]);
				return;
			}
			//.. e para quando penetra um novo corpo material ou umDetector
		}
		else
		{
			NCROSS = NCROSS + 1;
			//free(d_S[trackIndex]);
			//free(d_IS[trackIndex]);
			return;
		}
	L202:;
		d_stepsi2_(KB1,  NSC);
		// d_stepsi2_(KB1, NSC);
		goto L200;
		// Neste ponto, o programa saiu do ciclo DO.
		// L203:; indica o final do loop
	}
	// A particula sai do recinto.
L300:;

	DSP = 1.0e36;
	dg_TRACK_mod_[index].IBODY = dg_QTREE_.NBODYS + 1;
	dg_TRACK_mod_[index].MAT = 0;
	if (dg_TRACK_mod_[index].MAT == MAT0)
		DSEF = DSEF + DSP;
	dg_PENGEOM_mod_.DSTOT[trackIndex] = dg_PENGEOM_mod_.DSTOT[trackIndex] + DSP;
	dg_TRACK_mod_[index].X = dg_TRACK_mod_[index].X + DSP * dg_TRACK_mod_[index].U;
	dg_TRACK_mod_[index].Y = dg_TRACK_mod_[index].Y + DSP * dg_TRACK_mod_[index].V;
	dg_TRACK_mod_[index].Z = dg_TRACK_mod_[index].Z + DSP * dg_TRACK_mod_[index].W;
	//free(d_S[trackIndex]);
	//free(d_IS[trackIndex]);
}

__device__ void d_stepsi2_(int &KB, int &NSC)
{
	/*Calcula as interseções da trajetória com o limite
		Superfícies do corpo KB. Os cruzamentos são adicionados à lista e
		classificados em ordem decrescente.
		Esta sub-rotina funciona apenas quando chamada de dentro da sub-rotina STEP.*/

	// stepsi_(KB, S, IS, NSC);
	// return;

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	int KFL;
	int KS;
	double ABSA;
	double ABSB;
	double A, B, C;
	double FUZZL = 1.0e-12;
	double T1;
	double DISCR;
	double FUZZ;
	int IAMBIG;
	double R2A;
	double DELTA;
	double SH;
	double T2;
	int KMAX;
	double SMAX;
	int KKMAX;

	// Determine cruzamentos de superfície.

	for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++)
	{
		// printf("\nKSURF: %d\n",  dg_QTREE_.KSURF[NXG - 1][KB -1]);

		/*As interseções com uma determinada superfície são calculadas apenas uma vez.
		O ponteiro lateral de uma superfície deve ser alterado cada vez que o a superfície está cruzada.*/

		KFL = dg_QTREE_.KFLAG[KSS - 1][KB - 1];
		if (KFL > 4)
			goto L100;
		KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
		if (dg_QTREE_.KSP[trackIndex][KS - 1] != 0)
			goto L100;
		d_fsurf2_(KS, A, B, C);
		ABSA = fabs(A);
		ABSB = fabs(B);

		// Plano, unica raiz
		if (ABSA < 1.0e-36)
		{
			if (ABSB > 0.0e0)
			{
				if (C < -FUZZL)
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
				}
				else if (C > FUZZL)
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
				}
				else
				{
					if (B < 0.0e0)
					{
						dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
					}
					else
					{
						dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
					}
					goto L100;
				}
				T1 = -C / B;
				if (T1 > 0.0e0)
				{
					NSC = NSC + 1;
					d_IS[trackIndex][NSC - 1] = KS;
					d_S[trackIndex][NSC - 1] = T1;
				}
			}
			else
			{
				if (C < 0.0e0)
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
				}
				else
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
				}
			}

			// Superficie não plana, duas raizes
		}
		else
		{
			DISCR = B * B - 4.0e0 * A * C;
			FUZZ = FUZZL * DISCR / ABSA;
			if (C < -FUZZ)
			{
				IAMBIG = 0;
				dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
			}
			else if (C > FUZZ)
			{
				IAMBIG = 0;
				dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
			}
			else
			{
				IAMBIG = 1;
				if (B < 0.0e0)
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
				}
				else
				{
					dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
				}
			}

			if (DISCR < 1.0e-36)
				goto L100;

			if (IAMBIG == 0)
			{
				R2A = 0.5e0 / A;
				DELTA = sqrt(DISCR) * fabs(R2A);
				SH = -B * R2A;
				T1 = SH - DELTA;
				if (T1 > 0.0e0)
				{
					NSC = NSC + 1;
					d_IS[trackIndex][NSC - 1] = KS;
					d_S[trackIndex][NSC - 1] = T1;
				}
				T2 = SH + DELTA;
				if (T2 > 0.0e0)
				{
					NSC = NSC + 1;
					d_IS[trackIndex][NSC - 1] = KS;
					d_S[trackIndex][NSC - 1] = T2;
				}
			}
			else
			{
				if (B * A < 0.0e0)
				{
					R2A = 0.5e0 / A;
					DELTA = sqrt(DISCR) * fabs(R2A);
					SH = -B * R2A;
					T2 = SH + DELTA;
					NSC = NSC + 1;
					d_IS[trackIndex][NSC - 1] = KS;
					d_S[trackIndex][NSC - 1] = fmax(T2, 0.0e0);
				}
			}
		}
	L100:;
	}

	// Classifique as distâncias da superfície em ordem decrescente.
	// printf("\nNSC %d\n", NSC);
	if (NSC > 1)
	{
		for (int KI = 1; KI <= (NSC - 1); KI++)
		{
			SMAX = d_S[trackIndex][KI - 1];
			KMAX = KI;
			for (int KJ = (KI + 1); KJ <= NSC; KJ++)
			{
				if (d_S[trackIndex][KJ - 1] > SMAX)
				{
					SMAX = d_S[trackIndex][KJ - 1];
					KMAX = KJ;
				}
			}
			if (KMAX != KI)
			{
				SMAX = d_S[trackIndex][KI - 1];
				d_S[trackIndex][KI - 1] = d_S[trackIndex][KMAX - 1];
				d_S[trackIndex][KMAX - 1] = SMAX;
				KKMAX = d_IS[trackIndex][KI - 1];
				d_IS[trackIndex][KI - 1] = d_IS[trackIndex][KMAX - 1];
				d_IS[trackIndex][KMAX - 1] = KKMAX;
			}
		}
	}
}

__device__ void d_steplb2_(int &KB, int &IERR)
{

	/*Ajuda a encontrar o corpo ou módulo que contém os ponteiros laterais fornecidos para as superfícies analisadas.
	 A subrotina STEPLB funciona apenas quando invocada de dentro da sub-rotina STEP.
	 Ele se move através da árvore de módulos em um passo unico*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	int NLBOD;
	int KBS;
	int KS;
	int KF;
	int KBD;

	// Analisa o corpo ou o modulo atual.

	if (dg_QBODY_.KBOMO[KB - 1] == 0)
	{
		// Corpo
		NLBOD = dg_QBODY_.KBODY[NXG - 1][KB - 1];
		if (NLBOD > 0)
		{
			for (int KBB = 1; KBB <= NLBOD; KBB++)
			{
				KBS = dg_QBODY_.KBODY[KBB - 1][KB - 1];
				for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KBS - 1]; KSS++)
				{
					KS = dg_QTREE_.KSURF[KSS - 1][KBS - 1];
					KF = dg_QTREE_.KFLAG[KSS - 1][KBS - 1];
					if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS - 1] != KF))
						goto L100;
				}
				dg_TRACK_mod_[index].IBODY = KBS;
				if (dg_QTREE_.KDGHT[NXG - 1][dg_TRACK_mod_[index].IBODY - 1] > 1)
				{
					IERR = -1; // A particula está dentro de um módulo irmã
				}
				else
				{
					IERR = 0; // A particula está dentro de um corpo irma
					dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_[index].IBODY - 1];
				}
				return;
			L100:;
			}
		}

		for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++)
		{
			KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
			KF = dg_QTREE_.KFLAG[KSS - 1][KB - 1];
			if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS - 1] != KF))
				goto L300;
		}
		dg_TRACK_mod_[index].IBODY = KB;
		IERR = 0; //A partícula permanece no mesmo corpo.
		dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_[index].IBODY - 1];
		return;
	}
	else
	{
		// Modulo
		for (int KBB = 1; KBB <= dg_QTREE_.KDGHT[NXG - 1][KB - 1]; KBB++)
		{
			KBD = dg_QTREE_.KDGHT[KBB - 1][KB - 1];
			for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KBD - 1]; KSS++)
			{
				KS = dg_QTREE_.KSURF[KSS - 1][KBD - 1];
				KF = dg_QTREE_.KFLAG[KSS - 1][KBD - 1];
				if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS - 1] != KF))
					goto L200;
			}
			dg_TRACK_mod_[index].IBODY = KBD;
			if (KBD == KB)
			{
				IERR = 0; // A partícula permanece dentro do módulo atual.
				dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_[index].IBODY - 1];
			}
			else
			{
				if (dg_QTREE_.KDGHT[NXG - 1][KBD - 1] > 1)
				{
					IERR = -1; // A particula esta dentro de um submodulo
				}
				else
				{
					IERR = 0; // A particula esta dentro de um corpo simples;
					dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_[index].IBODY - 1];
				}
			}
			return;
		L200:;
		}
	}

	// A partícula está fora do corpo ou módulo atual.
L300:;
	IERR = 1;
	dg_TRACK_mod_[index].IBODY = dg_QTREE_.KMOTH[KB - 1];
	if (dg_TRACK_mod_[index].IBODY == 0)
	{
		dg_TRACK_mod_[index].IBODY = dg_QTREE_.NBODYS + 1;
		dg_TRACK_mod_[index].MAT = 0;
	}
}

__device__ void d_fsurf2_(int &KS, double &A, double &B, double &C)
{

	/*
	Calcula os par�metros da fun��o mestre da superf�cie KS e o raio (X, Y, Z) + S * (U, V, W).
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double XXX, YYY, ZZZ;

	if (dg_QSURF_.KPLANE[KS - 1] == 0)
	{
		A = dg_TRACK_mod_[index].U * (dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_[index].U + dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_[index].V + dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_[index].W) +
			dg_TRACK_mod_[index].V * (dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_[index].V + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_[index].W) + dg_TRACK_mod_[index].W * dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_[index].W;
		XXX = dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_[index].X + dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_[index].Y + dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_[index].Z + dg_QSURF_.AX[KS - 1];
		YYY = dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_[index].Y + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_[index].Z + dg_QSURF_.AY[KS - 1];
		ZZZ = dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_[index].Z + dg_QSURF_.AZ[KS - 1];

		B = dg_TRACK_mod_[index].U * (dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_[index].X + XXX) + dg_TRACK_mod_[index].V * (dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_[index].X + dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_[index].Y + YYY) +
			dg_TRACK_mod_[index].W * (dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_[index].X + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_[index].Y + dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_[index].Z + ZZZ);

		C = dg_TRACK_mod_[index].X * XXX + dg_TRACK_mod_[index].Y * YYY + dg_TRACK_mod_[index].Z * ZZZ + dg_QSURF_.A0[KS - 1];
	}
	else
	{
		A = 0.0e0;
		B = dg_TRACK_mod_[index].U * dg_QSURF_.AX[KS - 1] + dg_TRACK_mod_[index].V * dg_QSURF_.AY[KS - 1] + dg_TRACK_mod_[index].W * dg_QSURF_.AZ[KS - 1];
		C = dg_TRACK_mod_[index].X * dg_QSURF_.AX[KS - 1] + dg_TRACK_mod_[index].Y * dg_QSURF_.AY[KS - 1] + dg_TRACK_mod_[index].Z * dg_QSURF_.AZ[KS - 1] + dg_QSURF_.A0[KS - 1];
	}
}

__device__ void d_locate2_()
{
	/*

	Esta sub-rotina determina o corpo que cont�m o ponto com
   coordenadas (X, Y, Z). Os efeitos dos erros de arredondamento num�ricos s�o
   evitado considerando superf�cies difusas, que aumentam ou encolhem ligeiramente
   quando a part�cula os cruza.

   Valores de entrada (m�dulo TRACK_mod):
	  X, Y, Z ... coordenadas da part�cula,
	  U, V, W ... dire��o do movimento.

   Valores de sa�da (m�dulo TRACK_mod):
	  IBODY ..... corpo onde a part�cula se move,
	  MAT ...... material em IBODY,
					 MAT = 0, indica uma regi�o vazia


	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double A = 0.0;
	double B = 0.0;
	double C = 0.0;
	double FUZZ = 0.0;
	const double FUZZL = 1.0e-12;
	;
	int KS;
	int KB;
	int KF;
	double ABSA = 0.0;

	for (int I = 1; I <= dg_QSURF_.NSURF; I++)
	{
		dg_QTREE_.KSP[trackIndex][I - 1] = 0;
	}
	int KB0 = dg_QTREE_.NBODYS;

d100:
	for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB0 - 1]; KSS++)
	{
		KS = dg_QTREE_.KSURF[KSS - 1][KB0 - 1];
		if ((dg_QTREE_.KSP[trackIndex][KS - 1] != 0) || (dg_QTREE_.KFLAG[KSS - 1][KB0 - 1] > 4))
		{
			goto d101;
		}

		d_fsurf2_(KS, A, B, C);

		ABSA = fabs(A);

		if (ABSA > 1.0e-36)
			FUZZ = FUZZL * (B * B - 4.0e0 * A * C) / ABSA;
		else
			FUZZ = FUZZL * fabs(B);

		if (C < (-FUZZ))
		{
			dg_QTREE_.KSP[trackIndex][KS - 1] = 1;
		}
		else if (C > FUZZ)
		{
			dg_QTREE_.KSP[trackIndex][KS - 1] = 2;
		}
		else
		{
			if (B < 0.0e0)
				dg_QTREE_.KSP[trackIndex][KS - 1] = 1; // particula movendo-se para dentro
			else
				dg_QTREE_.KSP[trackIndex][KS - 1] = 2; // particula movendo-se para fora
		}
	d101:;
	}
	for (int KBB = 1; KBB <= dg_QTREE_.KDGHT[NXG - 1][KB0 - 1]; KBB++)
	{
		KB = dg_QTREE_.KDGHT[KBB - 1][KB0 - 1];
		for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++)
		{
			KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
			KF = dg_QTREE_.KFLAG[KSS - 1][KB - 1];

			if ((KF < 3) && (dg_QTREE_.KSP[trackIndex][KS - 1] != KF))
			{
				goto d102;
			}
		}
		if (KB == KB0)
		{
			dg_TRACK_mod_[index].IBODY = KB;						  // a particula está dentro do corpo ou modulo KB
			dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[KB - 1]; // a particula está dentro do MATERial KB
			return;
		}
		else if (dg_QTREE_.KDGHT[NXG - 1][KB - 1] > 1)
		{
			KB0 = KB; // o ponto está dentro de um submodulo
			goto d100;
		}
		else
		{
			dg_TRACK_mod_[index].IBODY = KB; // a particula esta dentro de um corpo ou modulo irmão
			dg_TRACK_mod_[index].MAT = dg_PENGEOM_mod_.MATER[KB - 1];
			return;
		}
	d102:;
	}
	dg_TRACK_mod_[index].IBODY = dg_QTREE_.NBODYS + 1;
	dg_TRACK_mod_[index].MAT = 0;
}

__device__ void d_sendet2_(double &ED, int &ID)
{
	/*
	Espectro de energia depositada.
	ED inclui o peso da partícula.
	*/
	int IE;

	dg_CENDET_.EDEP[ID - 1] = dg_CENDET_.EDEP[ID - 1] + ED;
	// dg_CENDET_.EDEP2[ID - 1] = dg_CENDET_.EDEP2[ID - 1] + pow(ED, 2);
	dg_CENDET_.EDEP2[ID - 1] = dg_CENDET_.EDEP2[ID - 1] + (ED * ED);
	if (ED > 1.0e-5)
	{
		if (dg_CENDET_.LLE[ID - 1] == 1)
		{
			IE = (int)(1.0e0 + (log(ED) - dg_CENDET_.EL[ID - 1]) * dg_CENDET_.RBSE[ID - 1]);
		}
		else
		{
			IE = (int)(1.0e0 + (ED - dg_CENDET_.EL[ID - 1]) * dg_CENDET_.RBSE[ID - 1]);
		}
		if ((IE > 0) && (IE <= dg_CENDET_.NE[ID - 1]))
			dg_CENDET_.DET[IE - 1][ID - 1] = dg_CENDET_.DET[IE - 1][ID - 1] + 1.0e0;
	}
}

__device__ void d_secpar2_(int &LEFT)
{
	/*
	Esta sub-rotina entrega o estado inicial de uma partícula secundária
	produzido durante a simulação anterior do chuveiro. Esta partícula
	é removido da pilha secundária, de modo que será perdido se um novo
	A chamada  para SECPAR é realizada antes de simular sua trajetória até
	o fim.

	LEFT é o número de partículas na pilha secundária na chamada
	hora . Quando LEFT=0, a simulação do chuveiro foi concluída.
	*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;
	/*	if (dg_SECST_.NSEC > 0) {
			LEFT = dg_SECST_.NSEC;
			dg_TRACK_mod_[index].E = dg_SECST_.ES[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].X = dg_SECST_.XS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].Y = dg_SECST_.YS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].Z = dg_SECST_.ZS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].U = dg_SECST_.US[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].V = dg_SECST_.VS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].W = dg_SECST_.WS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].WGHT = dg_SECST_.WGHTS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].KPAR = dg_SECST_.KS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].IBODY = dg_SECST_.IBODYS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].MAT = dg_SECST_.MS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].IPOL = dg_SECST_.IPOLS[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].SP1 = dg_SECST_.SP1S[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].SP2 = dg_SECST_.SP2S[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].SP3 = dg_SECST_.SP3S[dg_SECST_.NSEC - 1];
			dg_TRACK_mod_[index].PAGE = dg_SECST_.PAGES[dg_SECST_.NSEC - 1];
			for (int I = 1; I <= 5; I++) {
				dg_TRACK_mod_[index].ILB[I - 1] = dg_SECST_.ILBS[dg_SECST_.NSEC - 1][I - 1];
			}
			dg_SECST_.NSEC = dg_SECST_.NSEC - 1;
		}
		else {
			LEFT = 0;
		}*/
}

__device__ void d_pana2_(double &E, double &E1, double &CDT1, double &E2, double &CDT2, int &M)
{
	/*
	Simulação de aniquilação de pósitrons (em repouso ou em voo) em
	material M. Ei e CDTi são as energias e os cossenos da direção polar
	dos dois fótons de aniquilação.
	*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double GAM, GAM21, ANI, CHIMIN, RCHI, GT0, CHI, GREJ, DET, CHIP;

	// Pósitrons lentos (assumidos em repouso).

	if (E < dg_PENELOPE_mod_.EABS[M - 1][3 - 1])
	{
		E1 = 0.5e0 * (E + d_TREV);
		E2 = E1;
		CDT1 = -1.0e0 + 2.0e0 * d_rand2_(1.0e0);
		CDT2 = -CDT1;
	}
	else
	{
		/*
		Aniquilação em vôo (dois fótons com energia e direções
		determinado a partir da dcs e conservação energia-momento).
		*/
		GAM = 1.0e0 + fmax(E, 1.0e0) / d_REV;
		GAM21 = sqrt(GAM * GAM - 1.0e0);
		ANI = 1.0e0 + GAM;
		CHIMIN = 1.0e0 / (ANI + GAM21);
		RCHI = (1.0e0 - CHIMIN) / CHIMIN;
		GT0 = ANI * ANI - 2.0e0;
	L1:;
		CHI = CHIMIN * pow(RCHI, d_rand2_(2.0e0));
		GREJ = ANI * ANI * (1.0e0 - CHI) + GAM + GAM - 1.0e0 / CHI;
		if (d_rand2_(3.0e0) * GT0 > GREJ)
			goto L1;

		DET = E + d_TREV;
		E1 = CHI * DET;
		CDT1 = (ANI - 1.0e0 / CHI) / GAM21;
		CHIP = 1.0e0 - CHI;
		E2 = DET - E1;
		CDT2 = (ANI - 1.0e0 / CHIP) / GAM21;
	}
}

__device__ void d_paux2_()
{
	/*
	Mecanismo de interação auxiliar para fotons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	// printf("Warning: Subroutine PAUX has been entered.\n");
}

__device__ void d_tenang2_(int &IEXIT, int &N)
{
	/*
	Calcula energia e distribuições angulares de partículas emergentes,
	grava e carrega arquivos de despejo, acumula arquivos de despejo de diferentes
	é executado e grava os resultados.
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double RA2DE = 180.0e0 / d_PI;

	// Pontue as contribuições de uma nova partícula.

	// Distribuição de energia de partículas emergentes.

	int KEn, KTH, KPH;
	double THETA, PHI;

	if (dg_CENANG_.LLE == 1)
	{
		KEn = (int)(1.0e0 + (log(dg_TRACK_mod_[index].E) - dg_CENANG_.EL) * dg_CENANG_.RBSE);
	}
	else
	{
		KEn = (int)(1.0e0 + (dg_TRACK_mod_[index].E - dg_CENANG_.EL) * dg_CENANG_.RBSE);
	}
	if ((KEn > 0) && (KEn <= dg_CENANG_.NE))
	{
		if (N != dg_CENANG_.LPDE[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1])
		{
			dg_CENANG_.PDE[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDE[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_CENANG_.PDEP[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1];
			dg_CENANG_.PDE2[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDE2[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_CENANG_.PDEP[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1], 2);
			dg_CENANG_.PDEP[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_TRACK_mod_[index].WGHT;
			dg_CENANG_.LPDE[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = N;
		}
		else
		{
			dg_CENANG_.PDEP[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDEP[KEn - 1][IEXIT - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
		}
	}

	// Distribuição angular de partículas emergentes.
	THETA = acos(dg_TRACK_mod_[index].W);
	if (dg_CENANG_.LLTH == 1)
	{
		KTH = (int)(1.0e0 + (log(fmax(THETA, 1.0e-12) * RA2DE) - dg_CENANG_.THL) * dg_CENANG_.RBSTH);
		if (KTH < 1)
			KTH = 1;
	}
	else
	{
		KTH = (int)(1.0e0 + THETA * RA2DE * dg_CENANG_.RBSTH);
	}
	if (fabs(dg_TRACK_mod_[index].U) > 1.0e-16)
	{
		PHI = atan2(dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].U);
	}
	else if (fabs(dg_TRACK_mod_[index].V) > 1.0e-16)
	{
		PHI = atan2(dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].U);
	}
	else
	{
		PHI = 0.0e0;
	}
	if (PHI < 0.0e0)
		PHI = d_TWOPI + PHI;
	KPH = (int)(1.0e0 + PHI * RA2DE * dg_CENANG_.RBSPH);
	if (N != dg_CENANG_.LPDA[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1])
	{
		dg_CENANG_.PDA[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDA[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_CENANG_.PDAP[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1];
		dg_CENANG_.PDA2[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDA2[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] + pow(dg_CENANG_.PDAP[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1], 2);
		dg_CENANG_.PDAP[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_TRACK_mod_[index].WGHT;
		dg_CENANG_.LPDA[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] = N;
	}
	else
	{
		dg_CENANG_.PDAP[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] = dg_CENANG_.PDAP[KPH - 1][KTH - 1][dg_TRACK_mod_[index].KPAR - 1] + dg_TRACK_mod_[index].WGHT;
	}
}

__device__ double d_rand2_(double DUMMY)
{ // gerador de numeros aleatorios
	/*
	Esta é uma versão adaptada da sub-rotina RANECU escrita por F. James
	(Comput. Phys. Commun. 60 (1990) 329-344), que foi modificado para
	dá um único número aleatório em cada chamada.

	As 'sementes' ISEED1 e ISEED2 devem ser inicializadas no programa principal
	e transferido através do bloco comum nomeado /RSEED/.

	Alguns compiladores incorporam um gerador de números aleatórios intrínseco com
	o mesmo nome (mas com diferentes listas de argumentos). Para evitar conflitos,
	é aconselhável declarar RAND como uma função externa em todas as sub-
	programas em que o chamam.
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double USCALE = 1.0e0 / 2.147483563e9;
	int I1, I2, IZ;

	I1 = dg_RSEED2_[index % 1000].ISEED1 / 53668;
	dg_RSEED2_[index % 1000].ISEED1 = 40014 * (dg_RSEED2_[index % 1000].ISEED1 - I1 * 53668) - I1 * 12211;
	if (dg_RSEED2_[index % 1000].ISEED1 < 0)
		dg_RSEED2_[index % 1000].ISEED1 = dg_RSEED2_[index % 1000].ISEED1 + 2147483563;

	I2 = dg_RSEED2_[index % 1000].ISEED2 / 52774;
	dg_RSEED2_[index % 1000].ISEED2 = 40692 * (dg_RSEED2_[index % 1000].ISEED2 - I2 * 52774) - I2 * 3791;
	if (dg_RSEED2_[index % 1000].ISEED2 < 0)
		dg_RSEED2_[index % 1000].ISEED2 = dg_RSEED2_[index % 1000].ISEED2 + 2147483399;

	IZ = dg_RSEED2_[index % 1000].ISEED1 - dg_RSEED2_[index % 1000].ISEED2;
	if (IZ < 1)
		IZ = IZ + 2147483562;

	return IZ * USCALE;
}

/*__device__ double d_rand2_(double DUMMY) { //gerador de numeros aleatorios
	/*
	Esta é uma versão adaptada da sub-rotina RANECU escrita por F. James
	(Comput. Phys. Commun. 60 (1990) 329-344), que foi modificado para
	dá um único número aleatório em cada chamada.

	As 'sementes' ISEED1 e ISEED2 devem ser inicializadas no programa principal
	e transferido através do bloco comum nomeado /RSEED/.

	Alguns compiladores incorporam um gerador de números aleatórios intrínseco com
	o mesmo nome (mas com diferentes listas de argumentos). Para evitar conflitos,
	é aconselhável declarar RAND como uma função externa em todas as sub-
	programas em que o chamam.
	*/
//	int index = blockDim.x * blockIdx.x + threadIdx.x;
/*	double USCALE = 1.0e0 / 2.147483563e9;
	int I1, I2, IZ;

	I1 = dg_RSEED_.ISEED1 / 53668;
	dg_RSEED_.ISEED1 = 40014 * (dg_RSEED_.ISEED1 - I1 * 53668) - I1 * 12211;
	if (dg_RSEED_.ISEED1 < 0)
		dg_RSEED_.ISEED1 = dg_RSEED_.ISEED1 + 2147483563;

	I2 = dg_RSEED_.ISEED2 / 52774;
	dg_RSEED_.ISEED2 = 40692 * (dg_RSEED_.ISEED2 - I2 * 52774) - I2 * 3791;
	if (dg_RSEED_.ISEED2 < 0)
		dg_RSEED_.ISEED2 = dg_RSEED_.ISEED2 + 2147483399;

	IZ = dg_RSEED_.ISEED1 - dg_RSEED_.ISEED2;
	if (IZ < 1)
		IZ = IZ + 2147483562;

	return IZ * USCALE;

}*/

__device__ void d_schiff2_(double &B, double &G1, double &G2)
{

	/*
	Funções de triagem F1(B) e F2(B) no diferencial Bethe-Heitler
	seção transversal para produção de pares.
	*/

	double F1, F2, A0, B2;

	B2 = B * B;
	F1 = 2.0e0 - 2.0e0 * log(1.0e0 + B2);
	F2 = F1 - 6.666666666666666e-1;
	if (B < 1.0e-10)
	{
		F1 = F1 - d_TWOPI * B;
	}
	else
	{
		A0 = 4.0e0 * B * atan2(1.0e0, B);
		F1 = F1 - A0;
		F2 = F2 + 2.0e0 * B2 * (4.0e0 - A0 - 3.0e0 * log((1.0e0 + B2) / B2));
	}
	G1 = 0.5e0 * (3.0e0 * F1 - F2);
	G2 = 0.25e0 * (3.0e0 * F1 + F2);
}

__device__ void d_gaux2_()
{
	/*
	Mecanismo de interação auxiliar para fotons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	// printf("Warning: Subroutine GAUX has been entered.\n");
}

__device__ void d_peld2_(double &RNDC, double &RMU)
{

	/*
	Simulação de eventos elásticos rígidos de positrons. Seções transversais de
	a base de dados numérica ELSEPA.

	Valor do argumento :
	RNDC ... valor de corte do número aleatório uniforme
	(apenas eventos difíceis são simulados).
	RMU .... deflexão angular amostrada, =(1-CDT)/2.
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	static const int NP = 128;
	static const int NPM1 = NP - 1;

	double PK, RU, RR, PP, XX, AA, BB, D;
	int ITN, JE, I, J, K;

	// Ponto da rede de energia

	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;

	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Ponto
	RU = RNDC + d_rand2_(2.0e0) * (1.0e0 - RNDC);

	// Selection of the interval (binary search in a restricted interval).
	ITN = (int)(RU * NPM1 + 1);
	I = dg_CPELDB_.ITLP[dg_TRACK_mod_[index].MAT - 1][JE - 1][ITN - 1];
	J = dg_CPELDB_.ITUP[dg_TRACK_mod_[index].MAT - 1][JE - 1][ITN - 1];
	if ((J - I) < 2)
		goto L2;
L1:;
	K = (I + J) / 2;
	if (RU > dg_CPELDB_.PSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][K - 1])
	{
		I = K;
	}
	else
	{
		J = K;
	}

	if ((J - I) > 1)
		goto L1;

	// Amostragem da distribuição cumulativa inversa racional.
L2:;
	PP = dg_CPELDB_.PSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
	RR = RU - PP;
	if (RR > 1.0e-16)
	{
		XX = dg_CPELDB_.XSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		AA = dg_CPELDB_.ASP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		BB = dg_CPELDB_.BSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		D = dg_CPELDB_.PSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I + 1 - 1] - PP;
		RMU = XX + ((1.0e0 + AA + BB) * D * RR / (D * D + (AA * D + BB * RR) * RR)) * (dg_CPELDB_.XSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I + 1 - 1] - XX);
	}
	else
	{
		RMU = dg_CPELDB_.XSP[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
	}
}

__device__ void d_pina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC)
{

	/*
	Amostragem aleatória de colisões inelásticas duras de positrons.

	Modelo C Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IOSC .... índice do oscilador que foi 'ionizado'.

	*/

	/*if (imprimiu==0){
		printf("\n\nPINA2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	bool LDIST;

	static const int NO = 512;

	double WCCM, PK, TST, UK, WK, WTHR, WM, WKP, QKP, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, XHC, XHTOT, TS1, RK, PHI;
	double QS, Q, QTREV, G12, BHA1, BHA2, BHA3, BHA4;
	int IO, JO, IT, JE;

	WCCM = dg_PENELOPE_mod_.WCC[M - 1];

	if (WCCM > E)
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}
	// Ponto da rede de energia
	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;
	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Seleção do oscilador ativo.
	TST = d_rand2_(2.0e0);
	// Busca binaria
	IO = 1;
	JO = dg_CPINAC_.NPIN[M - 1] + 1;
L1:;
	IT = (IO + JO) / 2;
	if (TST > dg_CPINAC_.PINAC[IT - 1][JE - 1][M - 1])
	{
		IO = IT;
	}
	else
	{
		JO = IT;
	}
	if (JO - IO > 1)
		goto L1;
	IOSC = dg_CPINAC_.IPIN[IO - 1][M - 1];
	UK = dg_CEIN_.UI[IOSC - 1][M - 1];
	WK = dg_CEIN_.WRI[IOSC - 1][M - 1];

	if (UK > 1.0e-3)
	{
		WTHR = fmax(WCCM, UK);
	}
	else
	{
		WTHR = fmax(WCCM, WK);
	}

	if (E < (WTHR + 1.0e-6))
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	conchas internas são variadas para produzir um limiar suave.
	*/

	LDIST = true;
	if (UK > 1.0e-3)
	{
		WM = 3.0e0 * WK - 2.0e0 * UK;
		if (E > WM)
		{
			WKP = WK;
			QKP = UK;
		}
		else
		{
			WKP = (E + 2.0e0 * UK) / 3.0e0;
			QKP = UK * (E / WM);
			WM = E;
		}
		if (WCCM > WM)
			LDIST = false;
		WCMAX = E;
		WDMAX = fmin(WM, WCMAX);
		if (WTHR > WDMAX)
			LDIST = false;
	}
	else
	{
		if (WCCM > WK)
			LDIST = false;
		WKP = WK;
		QKP = WK;
		WM = E;
		WCMAX = E;
		WDMAX = WKP + 1.0e0;
	}

	// Constantes

	RB = E + d_TREV;
	GAM = 1.0e0 + E * d_RREV;
	GAM2 = GAM * GAM;
	BETA2 = (GAM2 - 1.0e0) / GAM2;
	G12 = pow((GAM + 1.0e0), 2);
	AMOL = pow(((GAM - 1.0e0) / GAM), 2);
	BHA1 = AMOL * (2.0e0 * G12 - 1.0e0) / (GAM2 - 1.0e0);
	BHA2 = AMOL * (3.0e0 + 1.0e0 / G12);
	BHA3 = AMOL * 2.0e0 * GAM * (GAM - 1.0e0) / G12;
	BHA4 = AMOL * pow((GAM - 1.0e0), 2) / G12;
	CPS = E * RB;
	CP = sqrt(CPS);

	// Seções transversais parciais do oscilador ativo.
	// Excitalçoes Distantes
	if (LDIST)
	{
		CPPS = (E - WKP) * (E - WKP + d_TREV);
		CPP = sqrt(CPPS);
		if (WKP > 1.0e-6 * E)
		{
			QM = sqrt(pow((CP - CPP), 2) + d_REV * d_REV) - d_REV;
		}
		else
		{
			QM = pow(WKP, 2) / (BETA2 * d_TREV);
			QM = QM * (1.0e0 - QM * d_RTREV);
		}
		if (QM < QKP)
		{
			RWKP = 1.0e0 / WKP;
			XHDL = log(QKP * (QM + d_TREV) / (QM * (QKP + d_TREV))) * RWKP;
			XHDT = fmax(log(GAM2) - BETA2 - DELTA, 0.0e0) * RWKP;
			if (UK > 1.0e-3)
			{
				F0 = (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR) / pow((WM - UK), 2);
				XHDL = F0 * XHDL;
				XHDT = F0 * XHDT;
			}
		}
		else
		{
			XHDL = 0.0e0;
			XHDT = 0.0e0;
		}
	}
	else
	{
		QM = 0.0e0; // Definido para evitar avisos de compilação.
		CPP = 0.0e0;
		CPPS = 0.0e0;
		XHDL = 0.0e0;
		XHDT = 0.0e0;
	}

	// Colisoes Fechadas
	RCL = WTHR / E;
	RL1 = 1.0e0 - RCL;
	XHC = ((1.0e0 / RCL - 1.0e0) + BHA1 * log(RCL) + BHA2 * RL1 + (BHA3 / 2.0e0) * (pow(RCL, 2) - 1.0e0) + (BHA4 / 3.0e0) * (1.0e0 - pow(RCL, 3))) / E;

	XHTOT = XHC + XHDL + XHDT;
	if (XHTOT < 1.0e-35)
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}

	// Amostragem de variáveis ​​de estado final.

	TST = d_rand2_(3.0e0) * XHTOT;

	// Colisão fechada dura

	TS1 = XHC;
	if (TST < TS1)
	{
	L2:;
		RK = RCL / (1.0e0 - d_rand2_(4.0e0) * RL1);
		PHI = 1.0e0 - RK * (BHA1 - RK * (BHA2 - RK * (BHA3 - BHA4 * RK)));
		if (d_rand2_(5.0e0) > PHI)
			goto L2;
		// Energia e ângulo de espalhamento (positron primário).
		DE = RK * E;
		EP = E - DE;
		CDT = sqrt(EP * RB / (E * (RB - DE)));
		// Energia e ângulo de emissão do raio delta.
		if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
		{
			if (UK > dg_CECUTR_.ECUTR[M - 1])
			{
				ES = DE - UK;
			}
			else
			{
				ES = DE;
			}
		}
		else
		{
			ES = DE;
		}
		CDTS = sqrt(DE * RB / (E * (DE + d_TREV)));
		return;
	}

	// Interação longitudinal dura distante.
	TS1 = TS1 + XHDL;
	if (UK > 1.0e-3)
	{
		DE = WM - sqrt(pow((WM - WTHR), 2) - d_rand2_(7.0e0) * (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR));
	}
	else
	{
		DE = WKP;
	}
	EP = E - DE;
	if (TST < TS1)
	{
		QS = QM / (1.0e0 + QM * d_RTREV);
		Q = QS / (pow(((QS / QKP) * (1.0e0 + QKP * d_RTREV)), d_rand2_(6.0e0)) - (QS * d_RTREV));
		QTREV = Q * (Q + d_TREV);
		CDT = (CPPS + CPS - QTREV) / (2.0e0 * CP * CPP);
		if (CDT > 1.0e0)
		{
			CDT = 1.0e0;
		}
		// Energia e ângulo de emissão do raio delta.
		if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
		{
			ES = DE - UK; // Apenas conchas internas.
		}
		else
		{
			ES = DE;
		}
		CDTS = 0.5e0 * (WKP * (E + RB - WKP) + QTREV) / sqrt(CPS * QTREV);
		if (CDTS > 1.0e0)
			CDTS = 1.0e0;
		return;
	}

	// Interação transversal distante difícil.
	CDT = 1.0e0;
	// Energia e ângulo de emissão do raio delta.
	if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
	{
		if (UK > dg_CECUTR_.ECUTR[M - 1])
		{
			ES = DE - UK; // Apenas conchas internas.
		}
		else
		{
			ES = DE;
		}
	}
	else
	{
		ES = DE;
	}
	CDTS = 1.0e0;
}

__device__ void d_psia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH)
{

	/*
		Amostragem aleatória de ionização de camada interna por impacto de positrons.

		Modelo Sternheimer-Liljequist GOS.

		Argumentos de entrada:
		E ....... energia do elétron (eV).
		M ....... material onde os elétrons se propagam.
		DELTA ... Correção do efeito de densidade de Fermi.
		Argumentos de saída:
		DE ...... perda de energia (eV).
		EP ...... energia do elétron espalhado (eV).
		CDT ..... cosseno do ângulo de dispersão polar.
		ES ...... energia do elétron secundário emitido (eV).
		CDTS .... cosseno polar de direção do elétron secundário.
		IZZ ..... número atômico do elemento onde a ionização ocorreu.
		ISH ..... camada de elétrons atômicos que foi ionizada.
	*/

	/*if (imprimiu==0){
		printf("\n\nPSIA2\n\n");
		imprimiu++;
	}*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double PK, TST, UK, WK, WTHR, WM, WKP, QKP, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, XHC, XHTOT, TS1, RK, PHI;
	double QS, Q, QTREV, G12, BHA1, BHA2, BHA3, BHA4;
	int IO, JO, IT, IOSC, JE;

	// Ponto da rede de energia

	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;
	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Seleção do oscilador ativo.

	TST = d_rand2_(2.0e0);

	// Busca binaria
	IO = 1;
	JO = dg_CPSIAC_.NPSI[M - 1] + 1;
L1:;
	IT = (IO + JO) / 2;
	if (TST > dg_CPSIAC_.PSIAC[IT - 1][JE - 1][M - 1])
	{
		IO = IT;
	}
	else
	{
		JO = IT;
	}
	if ((JO - IO) > 1)
		goto L1;

	IOSC = dg_CPSIAC_.IPSI[IO - 1][M - 1];
	IZZ = dg_CEIN_.KZ[IOSC - 1][M - 1];
	ISH = dg_CEIN_.KS[IOSC - 1][M - 1];
	UK = dg_CEIN_.UI[IOSC - 1][M - 1];
	WK = dg_CEIN_.WRI[IOSC - 1][M - 1];

	if (UK > 1.0e-3)
	{
		WTHR = UK;
	}
	else
	{
		WTHR = WK;
	}

	if (E < WTHR + 1.0e-6)
	{
		DE = UK;
		EP = E - DE;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	 conchas internas são variadas para produzir um limiar suave.
	*/

	WM = 3.0e0 * WK - 2.0e0 * UK;
	if (E > WM)
	{
		WKP = WK;
		QKP = UK;
	}
	else
	{
		WKP = (E + 2.0e0 * UK) / 3.0e0;
		QKP = UK * (E / WM);
		WM = E;
	}

	WCMAX = E;
	WDMAX = fmin(WM, WCMAX);

	// Constantes

	RB = E + d_TREV;
	GAM = 1.0e0 + E * d_RREV;
	GAM2 = GAM * GAM;
	BETA2 = (GAM2 - 1.0e0) / GAM2;
	G12 = pow((GAM + 1.0e0), 2);
	AMOL = pow(((GAM - 1.0e0) / GAM), 2);
	BHA1 = AMOL * (2.0e0 * G12 - 1.0e0) / (GAM2 - 1.0e0);
	BHA2 = AMOL * (3.0e0 + 1.0e0 / G12);
	BHA3 = AMOL * 2.0e0 * GAM * (GAM - 1.0e0) / G12;
	BHA4 = AMOL * pow((GAM - 1.0e0), 2) / G12;
	CPS = E * RB;
	CP = sqrt(CPS);

	// Seções transversais parciais do oscilador ativo.

	// Excitações distantes.

	CPPS = (E - WKP) * (E - WKP + d_TREV);
	CPP = sqrt(CPPS);

	if (WKP > 1.0e-6 * E)
	{
		QM = sqrt(pow((CP - CPP), 2) + d_REV * d_REV) - d_REV;
	}
	else
	{
		QM = pow(WKP, 2) / (BETA2 * d_TREV);
		QM = QM * (1.0e0 - QM * d_RTREV);
	}

	if (QM < QKP)
	{
		RWKP = 1.0e0 / WKP;
		XHDL = log(QKP * (QM + d_TREV) / (QM * (QKP + d_TREV))) * RWKP;
		XHDT = fmax(log(GAM2) - BETA2 - DELTA, 0.0e0) * RWKP;
		F0 = (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR) / pow((WM - UK), 2);
		XHDL = F0 * XHDL;
		XHDT = F0 * XHDT;
	}
	else
	{
		XHDL = 0.0e0;
		XHDT = 0.0e0;
	}

	// Colisoes Proximas
	RCL = WTHR / E;
	RL1 = 1.0e0 - RCL;

	XHC = ((1.0e0 / RCL - 1.0e0) + BHA1 * log(RCL) + BHA2 * RL1 + (BHA3 / 2.0e0) * (pow(RCL, 2) - 1.0e0) + (BHA4 / 3.0e0) * (1.0e0 - pow(RCL, 3))) / E;

	XHTOT = XHC + XHDL + XHDT;

	if (XHTOT < 1.0e-35)
	{
		DE = UK;
		EP = E - DE;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		return;
	}

	// Amostragem de variáveis ​​de estado final.

	TST = d_rand2_(3.0e0) * XHTOT;

	// Colisão aproximada dura.

	TS1 = XHC;
	if (TST < TS1)
	{
	L2:;
		RK = RCL / (1.0e0 - d_rand2_(4.0e0) * RL1);
		PHI = 1.0e0 - RK * (BHA1 - RK * (BHA2 - RK * (BHA3 - BHA4 * RK)));
		if (d_rand2_(5.0e0) > PHI)
			goto L2;
		// Energia e ângulo de espalhamento (elétron primário).
		DE = RK * E;
		EP = E - DE;
		CDT = sqrt(EP * RB / (E * (RB - DE)));
		// Energia e ângulo de emissão do raio delta.
		ES = DE - UK;
		CDTS = sqrt(DE * RB / (E * (DE + d_TREV)));
		return;
	}

	// Interação longitudinal dura distante.
	TS1 = TS1 + XHDL;
	DE = WM - sqrt(pow((WM - WTHR), 2) - d_rand2_(7.0e0) * (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR));
	EP = E - DE;

	if (TST < TS1)
	{
		QS = QM / (1.0e0 + QM * d_RTREV);
		Q = QS / (pow(((QS / QKP) * (1.0e0 + QKP * d_RTREV)), d_rand2_(6.0e0)) - (QS * d_RTREV));
		QTREV = Q * (Q + d_TREV);
		CDT = (CPPS + CPS - QTREV) / (2.0e0 * CP * CPP);
		if (CDT > 1.0e0)
			CDT = 1.0e0;

		// Energia e ângulo de emissão do raio delta.
		ES = DE - UK;
		CDTS = 0.5e0 * (WKP * (E + RB - WKP) + QTREV) / sqrt(CPS * QTREV);
		if (CDTS > 1.0e0)
			CDTS = 1.0e0;
		return;
	}

	// Interação transversal distante difícil.

	CDT = 1.0e0;

	// Energia e ângulo de emissão do raio delta.

	ES = DE - UK;
	CDTS = 1.0e0;
}

__device__ void d_gpha2_(double &ES, int &IZZ, int &ISH)
{

	/*
	   Simulação de absorção fotoelétrica no material M.

	   Argumentos de saída:
	   ES .... energia cinética do fotoelétron.
	   IZZ ... número atômico do átomo onde ocorreu a absorção.
	   ISH ... camada de elétrons atômica que foi ionizada.

	   NOTA: O JUMP usa uma seção transversal fotoelétrica ligeiramente maior
	   do que seu valor 'verdadeiro'. Para corrigir isso, o fóton pode
	   'sobrevive' a um evento fotoelétrico. A sobrevivência do fóton é sinalizada por
	   configuração IZZ=0, ISH=0, ES=0.0D0 (a energia E do fóton é mantida
	   inalterado.
	*/

	/*if (imprimiu==0){
		printf("\n\nGPHA2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double ACP[35];
	int IP[35];

	int I, IU, IT, IELAC, J;
	double PTOT, DEE, PCSL, TST, PIS, EBB;

	// Coeficientes de atenuação parcial.

	PTOT = 0.0e0;
	for (int IEL = 1; IEL <= dg_COMPOS_.NELEM[dg_TRACK_mod_[index].MAT - 1]; IEL++)
	{
		IZZ = dg_COMPOS_.IZ[IEL - 1][dg_TRACK_mod_[index].MAT - 1];
		// Busca Binaria
		I = dg_CGPH00_.IPHF[IZZ - 1];
		IU = dg_CGPH00_.IPHL[IZZ - 1];
	L1:;
		IT = (I + IU) / 2;
		if (dg_CEGRID_.XEL[trackIndex] > dg_CGPH00_.EPH[IT - 1])
			I = IT;
		else
			IU = IT;
		if (IU - I > 1)
			goto L1;

		IP[IEL - 1] = I;
		DEE = dg_CGPH00_.EPH[I + 1 - 1] - dg_CGPH00_.EPH[I - 1];
		if (DEE > 1.0e-15)
		{
			PCSL = dg_CGPH00_.XPH[1 - 1][I - 1] + (dg_CGPH00_.XPH[1 - 1][I + 1 - 1] - dg_CGPH00_.XPH[1 - 1][I - 1]) * (dg_CEGRID_.XEL[trackIndex] - dg_CGPH00_.EPH[I - 1]) / DEE;
		}
		else
		{
			PCSL = dg_CGPH00_.XPH[1 - 1][I - 1];
		}
		PTOT = PTOT + dg_COMPOS_.STF[IEL - 1][dg_TRACK_mod_[index].MAT - 1] * exp(PCSL);
		ACP[IEL - 1] = PTOT;
	}
	if (PTOT * dg_COMPOS_.VMOL[dg_TRACK_mod_[index].MAT - 1] > dg_CGIMFP_.SGPH[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1])
	{
		printf("WARNING: SGPH is less than the actual mac.\n");
	}

	// Faça uma amostra do elemento ativo.
	TST = d_rand2_(1.0e0) * dg_CGIMFP_.SGPH[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] / dg_COMPOS_.VMOL[dg_TRACK_mod_[index].MAT - 1];
	for (int IEL = 1; IEL <= dg_COMPOS_.NELEM[dg_TRACK_mod_[index].MAT - 1]; IEL++)
	{
		if (ACP[IEL - 1] > TST)
		{
			IELAC = IEL;
			IZZ = dg_COMPOS_.IZ[IEL - 1][dg_TRACK_mod_[index].MAT - 1];
			goto L2;
		}
	}

	// nteração delta. Introduzido para corrigir o uso de um limite superior do coeficiente de atenuação fotoelétrica.
	IZZ = 0;
	ISH = 0;
	ES = 0.0e0;
	return;

L2:;
	// Seleção do shell ativo.
	I = IP[IELAC - 1];
	DEE = dg_CGPH00_.EPH[I + 1 - 1] - dg_CGPH00_.EPH[I - 1];
	PIS = 0.0e0;

	if (DEE > 1.0e-15)
	{
		PTOT = exp(dg_CGPH00_.XPH[1 - 1][I - 1] + (dg_CGPH00_.XPH[1 - 1][I + 1 - 1] - dg_CGPH00_.XPH[1 - 1][I - 1]) * (dg_CEGRID_.XEL[trackIndex] - dg_CGPH00_.EPH[I - 1]) / DEE);
		TST = d_rand2_(2.0e0) * PTOT;

		for (int IS = 1; IS <= dg_CGPH00_.NPHS[IZZ - 1]; IS++)
		{
			J = IS + 1;
			PCSL = dg_CGPH00_.XPH[J - 1][I - 1] + (dg_CGPH00_.XPH[J - 1][I + 1 - 1] - dg_CGPH00_.XPH[J - 1][I - 1]) * (dg_CEGRID_.XEL[trackIndex] - dg_CGPH00_.EPH[I - 1]) / DEE;
			PIS = PIS + exp(PCSL);
			if (PIS > TST)
			{
				ISH = IS;
				goto L3;
			}
		}
	}
	else
	{
		PTOT = exp(dg_CGPH00_.XPH[1 - 1][I - 1]);
		TST = d_rand2_(2.0e0) * PTOT;
		for (int IS = 1; IS <= dg_CGPH00_.NPHS[IZZ - 1]; IS++)
		{
			PIS = PIS + exp(dg_CGPH00_.XPH[IS + 1 - 1][I - 1]);
			if (PIS > TST)
			{
				ISH = IS;
				goto L3;
			}
		}
	}
	ISH = 17;

	// Emissão do FotoEletron
L3:;
	if (ISH < 17)
	{
		EBB = dg_CADATA_.EB[ISH - 1][IZZ - 1];
		if (EBB > dg_CECUTR_.ECUTR[dg_TRACK_mod_[index].MAT - 1])
		{
			ES = dg_TRACK_mod_[index].E - EBB;
		}
		else
		{
			ES = dg_TRACK_mod_[index].E;
			ISH = 17;
		}
	}
	else
	{
		ES = dg_TRACK_mod_[index].E;
	}
}

__device__ void d_sauter2_(double &ES, double &CDTS)
{

	/*
	Amostragem aleatória da direção inicial dos fotoelétrons do
	Distribuição Sauter.
	*/

	/*if (imprimiu==0){
		printf("\n\nSAUTER2\n\n");
		imprimiu++;
	}*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double GAM, GAM2, BETA, AC, A1, A2, GTMAX, RU, TSAM, GTR;

	if (ES > 1.0e9)
	{
		CDTS = 1.0e0;
		return;
	}
	GAM = 1.0e0 + ES / d_REV;
	GAM2 = GAM * GAM;
	BETA = sqrt((GAM2 - 1.0e0) / GAM2);
	AC = 1.0e0 / BETA - 1.0e0;
	A1 = 0.5e0 * BETA * GAM * (GAM - 1.0e0) * (GAM - 2.0e0);
	A2 = AC + 2.0e0;
	GTMAX = 2.0e0 * (A1 + 1.0e0 / AC);
L1:;
	RU = d_rand2_(1.0e0);
	TSAM = 2.0e0 * AC * (2.0e0 * RU + A2 * sqrt(RU)) / (A2 * A2 - 4.0e0 * RU);
	GTR = (2.0e0 - TSAM) * (A1 + 1.0e0 / (AC + TSAM));
	if (d_rand2_(2.0e0) * GTMAX > GTR)
		goto L1;
	CDTS = 1.0e0 - TSAM;
}

__device__ void d_gppa2_(double &EE, double &CDTE, double &EP, double &CDTP, int &IZZ, int &ISH)
{

	/*
	   Amostragem aleatória do par elétron-pósitron e produção de tripletos por
	   fótons. Seção transversal diferencial Bethe-Heitler.

	   Valores de saída :
	   EE ..... energia cinética do elétron.
	   CDTE ... direção polar cosseno do elétron.
	   EP ..... energia cinética do pósitron.
	   CDTP ... cosseno de direção polar do pósitron.
	   IZZ ... número atômico do átomo onde ocorreu a absorção.
	   ISH .... camada de elétrons atômicos que foi ionizada.
	*/

	/*if (imprimiu==0){
		printf("\n\nGPPA2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double SL = 137.035999074e0;

	double EKI, EPS, ALZ, T, F00, G0, BMIN, G1MIN, G2MIN, XR, A1, A2, P1, RU2M1, B, G1, G2, TRIPL, TST;
	int I, ISHELL, JO;

	EKI = d_REV / dg_TRACK_mod_[index].E;

	if (dg_TRACK_mod_[index].E < 1.1e6)
	{
		EPS = EKI + (1.0e0 - 2.0e0 * EKI) * d_rand2_(1.0e0);
		goto L3;
	}

	// Low-energy and Coulomb corrections.
	ALZ = dg_CGPP00_.ZEQPP[dg_TRACK_mod_[index].MAT - 1] / SL;
	T = sqrt(2.0e0 * EKI);
	F00 = (-1.774e0 - 1.210e1 * ALZ + 1.118e1 * ALZ * ALZ) * T + (8.523e0 + 7.326e1 * ALZ - 4.441e1 * ALZ * ALZ) * pow(T, 2) - (1.352e1 + 1.211e2 * ALZ - 9.641e1 * ALZ * ALZ) * pow(T, 3) + (8.946e0 + 6.205e1 * ALZ - 6.341e1 * ALZ * ALZ) * pow(T, 4);
	G0 = dg_CGPP00_.F0[2 - 1][dg_TRACK_mod_[index].MAT - 1] + F00;
	BMIN = 4.0e0 * EKI / dg_CGPP00_.BCB[dg_TRACK_mod_[index].MAT - 1];
	d_schiff2_(BMIN, G1, G2);
	G1MIN = G1 + G0;
	G2MIN = G2 + G0;
	XR = 0.5e0 - EKI;
	A1 = 6.666666666666666e-1 * G1MIN * pow(XR, 2);
	P1 = A1 / (A1 + G2MIN);

	// Amostragem aleatória de EPS.

L1:;
	if (d_rand2_(2.020) > P1)
		goto L2;
	RU2M1 = 2.0e0 * d_rand2_(3.0e0) - 1.0e0;
	if (RU2M1 < 0.0e0)
		EPS = 0.5e0 - XR * pow(fabs(RU2M1), 3.333333333333333e-1);
	else
		EPS = 0.5e0 + XR * pow(RU2M1, 3.333333333333333e-1);
	B = EKI / (dg_CGPP00_.BCB[dg_TRACK_mod_[index].MAT - 1] * EPS * (1.0e0 - EPS));
	d_schiff2_(B, G1, G2);
	G1 = fmax(G1 + G0, 0.0e0);
	if (d_rand2_(4.0e0) * G1MIN > G1)
		goto L1;
	goto L3;
L2:;
	EPS = EKI + 2.0e0 * XR * d_rand2_(5.0e0);
	B = EKI / (dg_CGPP00_.BCB[dg_TRACK_mod_[index].MAT - 1] * EPS * (1.0e0 - EPS));
	d_schiff2_(B, G1, G2);
	G2 = fmax(G2 + G0, 0.0e0);
	if (d_rand2_(6.0e0) * G2MIN > G2)
		goto L1;
L3:;
	// Eletron
	EE = EPS * dg_TRACK_mod_[index].E - d_REV;
	CDTE = 2.0e0 * d_rand2_(7.0e0) - 1.0e0;
	A1 = EE + d_REV;
	A2 = sqrt(EE * (EE + d_TREV));
	CDTE = (CDTE * A1 + A2) / (A1 + CDTE * A2);
	// Positron
	EP = (1.0e0 - EPS) * dg_TRACK_mod_[index].E - d_REV;
	CDTP = 2.0e0 * d_rand2_(8.0e0) - 1.0e0;
	A1 = EP + d_REV;
	A2 = sqrt(EP * (EP + d_TREV));
	CDTP = (CDTP * A1 + A2) / (A1 + CDTP * A2);

	// Produção de Tripletos
	TRIPL = dg_CGPP01_.TRIP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CGPP01_.TRIP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CGPP01_.TRIP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	IZZ = 0;
	ISH = 30;
	if (TRIPL < 1.0e-5)
		return;
	if (d_rand2_(9.0e0) > TRIPL)
		return;
	TST = d_rand2_(10.0e0);
	// Busca binaria
	if (TST < dg_CGCO_.PTRSH[1 - 1][dg_TRACK_mod_[index].MAT - 1])
	{
		ISHELL = 1;
	}
	else
	{
		ISHELL = 1;
		JO = dg_CGCO_.NOSCCO[dg_TRACK_mod_[index].MAT - 1] + 1;
	L4:;
		I = (ISHELL + JO) / 2;
		if (TST > dg_CGCO_.PTRSH[I - 1][dg_TRACK_mod_[index].MAT - 1])
		{
			ISHELL = I;
		}
		else
		{
			JO = I;
		}
		if (JO - ISHELL > 1)
			goto L4;
		ISHELL = ISHELL + 1;
	}
	IZZ = dg_CGCO_.KZCO[ISHELL - 1][dg_TRACK_mod_[index].MAT - 1];
	ISH = dg_CGCO_.KSCO[ISHELL - 1][dg_TRACK_mod_[index].MAT - 1];
}

__device__ void d_eaux2_()
{
	/*
	Mecanismo de interação auxiliar para elétrons, definível pelo usuário.
	Geralmente não está ativo.
	*/
	// printf("Warning: Subroutine EAUX has been entered.\n");
}

__device__ void d_graa2_(double &E, double &CDT, int &IEFF, int &M)
{

	/*if (imprimiu==0){
		printf("\n\nGRAA2\n\n");
		imprimiu++;
	}*/

	// Amostragem aleatória de espalhamento coerente (Rayleigh)
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	static const int NP = 150;
	static const int NPM1 = NP - 1;

	int II, IU, IT, ITN, I, J, K;
	double XSE, QMAX, Q2MAX, G, RU, RR, D, XX;

	// Busca Binaria

	II = dg_CGRA01_.IED[dg_CEGRID_.KE[trackIndex] - 1];
	IU = dg_CGRA01_.IEU[dg_CEGRID_.KE[trackIndex] - 1];
L1:;
	IT = (II + IU) / 2;

	if (dg_CEGRID_.XEL[trackIndex] > dg_CGRA01_.ERA[IT - 1])
		II = IT;
	else
		IU = IT;
	if (IU - II > 1)
		goto L1;

	XSE = exp(dg_CGRA01_.XSRA[II - 1][M - 1] + (dg_CGRA01_.XSRA[II + 1 - 1][M - 1] - dg_CGRA01_.XSRA[II - 1][M - 1]) * (dg_CEGRID_.XEL[trackIndex] - dg_CGRA01_.ERA[II - 1]) / (dg_CGRA01_.ERA[II + 1 - 1] - dg_CGRA01_.ERA[II - 1]));

	if (d_rand2_(1.0e0) * dg_CGIMFP_.SGRA[dg_CEGRID_.KE[trackIndex] - 1][M - 1] > XSE)
	{
		IEFF = 0;
		CDT = 1.0e0;
		return;
	}

	IEFF = 1;
	QMAX = 2.0e0 * E * d_RREV;

	if (QMAX < 1.0e-10)
	{
	L2:;
		CDT = 1.0e0 - 2.0e0 * d_rand2_(1.0e0);
		G = 0.5e0 * (1.0e0 + CDT * CDT);
		if (d_rand2_(2.0e0) > G)
			goto L2;
		return;
	}
	Q2MAX = fmin(QMAX * QMAX, dg_CGRA03_.QRA[M - 1][NP - 1]);

L3:;
	RU = d_rand2_(3.0e0) * dg_CGRA03_.PMAX[M - 1][dg_CEGRID_.KE[trackIndex] + 1 - 1];

	// Seleção do intervalo (busca binária dentro de limites pré-calculados).

	ITN = (int)(RU * NPM1 + 1);
	I = dg_CGRA03_.ITLRA[M - 1][ITN - 1];
	J = dg_CGRA03_.ITURA[M - 1][ITN - 1];

	if ((J - I) < 2)
		goto L5;
L4:;
	K = (I + J) / 2;
	if (RU > dg_CGRA03_.PRA[M - 1][K - 1])
		I = K;
	else
		J = K;
	if ((J - I) > 1)
		goto L4;

	// Amostragem da distribuição cumulativa inversa racional.

L5:;
	RR = RU - dg_CGRA03_.PRA[M - 1][I - 1];
	if (RR > 1.0e-16)
	{
		D = dg_CGRA03_.DPRA[M - 1][I - 1];
		XX = dg_CGRA03_.QRA[M - 1][I - 1] + ((1.0e0 + dg_CGRA03_.ARA[M - 1][I - 1] + dg_CGRA03_.BRA[M - 1][I - 1]) * D * RR / (D * D + (dg_CGRA03_.ARA[M - 1][I - 1] * D + dg_CGRA03_.BRA[M - 1][I - 1] * RR) * RR)) * (dg_CGRA03_.QRA[M - 1][I + 1 - 1] - dg_CGRA03_.QRA[M - 1][I - 1]);
	}
	else
	{
		XX = dg_CGRA03_.QRA[M - 1][I - 1];
	}
	if (XX > Q2MAX)
		goto L3;
	CDT = 1.0e0 - 2.0e0 * XX / Q2MAX;
	// Rejeição
	G = 0.5e0 * (1.0e0 + CDT * CDT);
	if (d_rand2_(4.0e0) > G)
		goto L3;
}

__device__ void d_dirpol2_(double &CDT, double &DF, double &CONS, double &SP1, double &SP2, double &SP3, double &U, double &V, double &W)
{

	/*
	Esta sub-rotina calcula os cossenos de direção _e_ os Stokes
	Parâmetros de um fóton polarizado após espalhamento com um dado polar
	ângulo.

	Entrada: U,V,W ... cossenos da direção inicial.
	SP1,SP2,SP3 ... parâmetros Stokes iniciais.
	CDT ..... cosseno do ângulo de dispersão polar.
	CONS .... constante no PDF do ângulo azimutal.
	Saída: U,V,W ... novos cossenos de direção.
	SP1,SP2,SP3 ... novos parâmetros Stokes.
	DF ...... ângulo de dispersão azimutal.
	CDT e CONS permanecem inalterados.
	*/

	/*if (imprimiu==0){
		printf("\n\nDIRPOL2\n\n");
		imprimiu++;
	}*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double SP1P, RSP0, CDT2, CDT21, PHA, PHB, SP0MAX, SDF, CDF, S2DF, C2DF, SP3P, SP0P, UV, UVW, FNORM, SDT, SDTSDF, SDTCDF, SUV, UN, VN;

	// Amostragem do ângulo de espalhamento azimutal.

	CDT2 = CDT * CDT;
	CDT21 = CDT2 + 1.0e0;
	PHA = CDT21 + CONS;
	PHB = 1.0e0 - CDT2;
	SP0MAX = PHA + PHB * sqrt(SP1 * SP1 + SP3 * SP3 + 1.0e-35);
L1:;
	DF = d_rand2_(1.0e0) * d_TWOPI;
	SDF = sin(DF);
	CDF = cos(DF);
	S2DF = 2.0e0 * SDF * CDF;
	C2DF = CDF * CDF - SDF * SDF;
	SP3P = S2DF * SP1 + C2DF * SP3; // Parâmetro Stokes com novo zero azimute.
	SP0P = PHA - PHB * SP3P;

	if (d_rand2_(2.0e0) * SP0MAX > SP0P)
		goto L1;

	// Calcular novos parâmetros Stokes
	SP1P = C2DF * SP1 - S2DF * SP3; // Parâmetro Stokes com novo zero azimute.
	RSP0 = 1.0e0 / SP0P;
	SP1 = 2.0e0 * CDT * SP1P * RSP0;
	SP2 = (2.0e0 + CONS) * CDT * SP2 * RSP0;
	SP3 = (CDT21 * SP3P - PHB) * RSP0;

	// Garanta a normalizacao
	UV = U * U + V * V;
	UVW = UV + W * W;
	if (fabs(UVW - 1.0e0) > 1.0e-13)
	{
		FNORM = 1.0e0 / sqrt(UVW);
		U = FNORM * U;
		V = FNORM * V;
		W = FNORM * W;
		UV = U * U + V * V;
	}

	// Calcula a nova direção
	if (1.0e0 - fabs(CDT) > 1.0e-8)
		SDT = sqrt(PHB);
	else
		SDT = sqrt(2.0e0 * (1.0e0 - fabs(CDT)));

	if (SDT < 1.0e-13)
	{
		if (CDT < 0.0e0)
		{
			U = -U;
			V = -V;
			W = -W;
		}
	}
	else
	{
		SDTSDF = SDT * SDF;
		SDTCDF = SDT * CDF;
		if (UV > 1.0e-26)
		{
			SUV = sqrt(UV);
			UN = U / SUV;
			VN = V / SUV;
			U = U * CDT + (UN * W * SDTCDF - VN * SDTSDF);
			V = V * CDT + (VN * W * SDTCDF + UN * SDTSDF);
			W = W * CDT - SUV * SDTCDF;
		}
		else
		{
			if (W > 0.0e0)
			{
				U = SDTCDF;
				V = SDTSDF;
				W = CDT;
			}
			else
			{
				U = -SDTCDF;
				V = -SDTSDF;
				W = -CDT;
			}
		}
	}
}

__device__ void d_gcoa2_(double &E, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH)
{

	/*
	   Amostragem aleatória de espalhamento incoerente (Compton) de fótons.

	   Argumentos de entrada:
	   E ..... energia do fóton incidente (eV).
	   M ..... material onde os fótons se propagam.
	   Argumento de saída:
	   DE .... perda de energia (eV).
	   EP .... energia do fóton espalhado (eV).
	   CDT ... cosseno do ângulo de espalhamento polar.
	   ES .... energia do elétron emitido (eV).
	   CDTS .. cosseno polar de direção do elétron.
	   IZZ ... número atômico do átomo onde ocorreu a dispersão.
	   ISH ... camada de elétrons atômica que foi ionizada.
	*/

	/*if (imprimiu==0){
		printf("\n\nGCOA2\n\n");
		imprimiu++;
	}*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double D2 = 1.4142135623731e0;
	double D1 = 1.0e0 / D2;
	double D12 = 0.5e0;

	static const int NOCO = 512;

	double RN[NOCO];
	double PAC[NOCO];

	double EK, EK2, EKS, EK1, TAUMIN, TAUM2, A1, A2, S0, AUX, PZOMC, RNI, TAU, CDT1, S, TST, A, XQC, AF;
	double FPZ, FPZMAX, T, B1, B2, Q2;
	int I2, ISHELL, I3, JO;

	EK = E * d_RREV;
	EK2 = EK + EK + 1.0e0;
	EKS = EK * EK;
	EK1 = EKS - EK2 - 1.0e0;
	TAUMIN = 1.0e0 / EK2;
	TAUM2 = TAUMIN * TAUMIN;
	A1 = log(EK2);
	A2 = A1 + 2.0e0 * EK * (1.0e0 + EK) * TAUM2;

	if (E > 5.0e6)
		goto L4;

	// Função de dispersão incoerente para theta=d_PI

	S0 = 0.0e0;
	for (int I = 1; I <= dg_CGCO_.NOSCCO[M - 1]; I++)
	{
		if (dg_CGCO_.UICO[I - 1][M - 1] < E)
		{
			AUX = E * (E - dg_CGCO_.UICO[I - 1][M - 1]) * 2.0e0;
			PZOMC = dg_CGCO_.FJ0[I - 1][M - 1] * (AUX - d_REV * dg_CGCO_.UICO[I - 1][M - 1]) / (d_REV * sqrt(AUX + AUX + pow(dg_CGCO_.UICO[I - 1][M - 1], 2)));
			if (PZOMC > 0.0e0)
				RNI = 1.0e0 - 0.5e0 * exp(D12 - pow((D1 + D2 * PZOMC), 2));
			else
				RNI = 0.5e0 * exp(D12 - pow((D1 - D2 * PZOMC), 2));
			S0 = S0 + dg_CGCO_.FCO[I - 1][M - 1] * RNI;
		}
	}

	// Amostragem de tau.

L1:;
	if (d_rand2_(1.0e0) * A2 < A1)
		TAU = pow(TAUMIN, d_rand2_(2.0e0));
	else
		TAU = sqrt(1.0e0 + d_rand2_(3.0e0) * (TAUM2 - 1.0e0));
	CDT1 = (1.0e0 - TAU) / (EK * TAU);

	// Função de dispersão incoerente.
	S = 0.0e0;
	for (int I = 1; I <= dg_CGCO_.NOSCCO[M - 1]; I++)
	{
		if (dg_CGCO_.UICO[I - 1][M - 1] < E)
		{
			AUX = E * (E - dg_CGCO_.UICO[I - 1][M - 1]) * CDT1;
			PZOMC = dg_CGCO_.FJ0[I - 1][M - 1] * (AUX - d_REV * dg_CGCO_.UICO[I - 1][M - 1]) / (d_REV * sqrt(AUX + AUX + pow(dg_CGCO_.UICO[I - 1][M - 1], 2)));
			if (PZOMC > 0.0e0)
				RN[I - 1] = 1.0e0 - 0.5e0 * exp(D12 - pow((D1 + D2 * PZOMC), 2));
			else
				RN[I - 1] = 0.5e0 * exp(D12 - pow((D1 - D2 * PZOMC), 2));
			S = S + dg_CGCO_.FCO[I - 1][M - 1] * RN[I - 1];
			PAC[I - 1] = S;
		}
		else
		{
			PAC[I - 1] = S;
		}
	}

	// Funcao de Rejeição
	TST = S * (1.0e0 + TAU * (EK1 + TAU * (EK2 + TAU * EKS))) / (EKS * TAU * (1.0e0 + TAU * TAU));
	if (d_rand2_(4.0e0) * S0 > TST)
		goto L1;
	CDT = 1.0e0 - CDT1;

	// Escudo do elétron alvo.

L2:
	TST = S * d_rand2_(5.0e0);
	// Busca Binaria
	if (TST < PAC[1 - 1])
	{
		ISHELL = 1;
	}
	else
	{
		ISHELL = 1;
		JO = dg_CGCO_.NOSCCO[M - 1] + 1;
	L3:;
		I2 = (ISHELL + JO) / 2;
		if (TST > PAC[I2 - 1])
			ISHELL = I2;
		else
			JO = I2;
		if (JO - ISHELL > 1)
			goto L3;
		ISHELL = ISHELL + 1;
	}

	// Momento projetado do elétron alvo.
	A = d_rand2_(6.0e0) * RN[ISHELL - 1];
	if (A < 0.5e0)
		PZOMC = (D1 - sqrt(D12 - log(A + A))) / (D2 * dg_CGCO_.FJ0[ISHELL - 1][M - 1]);
	else
		PZOMC = (sqrt(D12 - log(2.0e0 - A - A)) - D1) / (D2 * dg_CGCO_.FJ0[ISHELL - 1][M - 1]);
	if (PZOMC < -1.0e0)
		goto L2;

	// Rejeição F(EP).
	XQC = 1.0e0 + TAU * (TAU - 2.0e0 * CDT);
	AF = sqrt(XQC) * (1.0e0 + TAU * (TAU - CDT) / XQC);

	if (AF > 0.0e0)
		FPZMAX = 1.0e0 + AF * 0.2e0;
	else
		FPZMAX = 1.0e0 - AF * 0.2e0;

	FPZ = 1.0e0 + AF * fmax(fmin(PZOMC, 0.2e0), -0.2e0);
	if (d_rand2_(7.0e0) * FPZMAX > FPZ)
		goto L2;

	// Energia do fóton espalhado
	T = pow(PZOMC, 2);
	B1 = 1.0e0 - T * TAU * TAU;
	B2 = 1.0e0 - T * TAU * CDT;
	if (PZOMC > 0.0e0)
		EP = E * (TAU / B1) * (B2 + sqrt(fabs(B2 * B2 - B1 * (1.0e0 - T))));
	else
		EP = E * (TAU / B1) * (B2 - sqrt(fabs(B2 * B2 - B1 * (1.0e0 - T))));
	goto L6;

	// Sem alargamento Doppler para E maior que 5 MeV.

L4:;
	if (d_rand2_(8.0e0) * A2 < A1)
		TAU = pow(TAUMIN, d_rand2_(9.0e0));
	else
		TAU = sqrt(1.0e0 + d_rand2_(10.0e0) * (TAUM2 - 1.0e0));

	// Funcao de Rejeição
	TST = (1.0e0 + TAU * (EK1 + TAU * (EK2 + TAU * EKS))) / (EKS * TAU * (1.0e0 + TAU * TAU));
	if (d_rand2_(11.0e0) > TST)
		goto L4;
	EP = TAU * E;
	CDT = 1.0e0 - (1.0e0 - TAU) / (EK * TAU);

	// Camada eletrônica alvo.
	TST = d_rand2_(12.0e0);

	// busca Binaria
	if (TST < dg_CGCO_.PTRSH[1 - 1][M - 1])
	{
		ISHELL = 1;
	}
	else
	{
		ISHELL = 1;
		JO = dg_CGCO_.NOSCCO[M - 1] + 1;
	L5:;
		I3 = (ISHELL + JO) / 2;
		if (TST > dg_CGCO_.PTRSH[I3 - 1][M - 1])
			ISHELL = I3;
		else
			JO = I3;
		if (JO - ISHELL > 1)
			goto L5;
		ISHELL = ISHELL + 1;
	}

	if (EP > (E - dg_CGCO_.UICO[ISHELL - 1][M - 1]))
		goto L4;

L6:;
	DE = E - EP;
	if (dg_CGCO_.KSCO[ISHELL - 1][M - 1] < 17)
	{
		if (dg_CGCO_.UICO[ISHELL - 1][M - 1] > dg_CECUTR_.ECUTR[M - 1])
		{
			ES = DE - dg_CGCO_.UICO[ISHELL - 1][M - 1];
		}
		else
		{
			ES = DE;
		}
	}
	else
	{
		ES = DE;
	}

	Q2 = E * E + EP * (EP - 2.0e0 * E * CDT);
	if (Q2 > 1.0e-12)
		CDTS = (E - EP * CDT) / sqrt(Q2);
	else
		CDTS = 1.0e0;

	IZZ = dg_CGCO_.KZCO[ISHELL - 1][M - 1];
	ISH = dg_CGCO_.KSCO[ISHELL - 1][M - 1];
}

__device__ void d_ebraa2_(double &E, double &DE, double &CDT, int &M)
{
	/*
		Amostragem aleatória da direção inicial dos fótons bremss, relativa
	na direção do projétil.
	Ajuste numérico/interpolação de funções de forma de onda parcial dadas por
	Kissel, Quarles e Pratt; ANDT 28(1993)381.

	Parâmetros de entrada:
	M ..... material onde o projétil se move.
	E ..... energia cinética do projétil.
	DE .... energia do fóton emitido.
	Parâmetro de saída:
	CDT ... cosseno do ângulo de emissão polar.
	*/

	/*if (imprimiu==0){
		printf("\n\nEBRAA2\n\n");
		imprimiu++;
	}*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double BETA, RK, P10, P11, P1, P20, P21, P2, BETAP;
	int IE, IE1, IET, IK;

	// Parametros de distribuição

	BETA = sqrt(E * (E + d_TREV)) / (E + d_REV);

	// Uma distribuição dipolar pura é usada para E>500 keV.

	if (E > 500.0e3)
	{
		CDT = 2.0e0 * d_rand2_(1.0e0) - 1.0e0;
		if (d_rand2_(2.0e0) > 0.75e0)
		{
			if (CDT < 0.0e0)
			{
				CDT = -pow((-CDT), 0.333333333333333e0);
			}
			else
			{
				CDT = pow(CDT, 0.333333333333333e0);
			}
		}
		CDT = (CDT + BETA) / (1.0e0 + BETA * CDT);
		return;
	}

	if (BETA > dg_CBRANG_.BET[6 - 1])
	{
		IE = 6;
		goto L20;
	}
	if (BETA < dg_CBRANG_.BET[1 - 1])
	{
		IE = 1;
		goto L20;
	}
	IE = 1;
	IE1 = 6;
L10:;
	IET = (IE + IE1) / 2;
	if (BETA > dg_CBRANG_.BET[IET - 1])
	{
		IE = IET;
	}
	else
	{
		IE1 = IET;
	}
	if ((IE1 - IE) > 1)
		goto L10;
L20:;

	RK = 1.0e0 + 20.0e0 * DE / E;
	IK = min(int(RK), 20);

	P10 = dg_CBRANG_.BP1[1 - 1][IK - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP1[2 - 1][IK - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP1[3 - 1][IK - 1][IE - 1][M - 1] + BETA * dg_CBRANG_.BP1[4 - 1][IK - 1][IE - 1][M - 1]));
	P11 = dg_CBRANG_.BP1[1 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP1[2 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP1[3 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * dg_CBRANG_.BP1[4 - 1][IK + 1 - 1][IE - 1][M - 1]));
	P1 = P10 + (RK - IK) * (P11 - P10);

	P20 = dg_CBRANG_.BP2[1 - 1][IK - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP2[2 - 1][IK - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP2[3 - 1][IK - 1][IE - 1][M - 1] + BETA * dg_CBRANG_.BP2[4 - 1][IK - 1][IE - 1][M - 1]));
	P21 = dg_CBRANG_.BP2[1 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP2[2 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * (dg_CBRANG_.BP2[3 - 1][IK + 1 - 1][IE - 1][M - 1] + BETA * dg_CBRANG_.BP2[4 - 1][IK + 1 - 1][IE - 1][M - 1]));
	P2 = P20 + (RK - IK) * (P21 - P20);

	// Amostragem das distribuições dipolo transformadas por Lorentz.

	P1 = fmin(exp(P1) / BETA, 1.0e0);
	BETAP = fmin(fmax(BETA * (1.0e0 + P2 / BETA), 0.0e0), 0.999999999e0);

	if (d_rand2_(3.0e0) < P1)
	{
	L1:;
		CDT = 2.0e0 * d_rand2_(4.0e0) - 1.0e0;
		if ((2.0e0 * d_rand2_(5.0e0)) > (1.0e0 + CDT * CDT))
			goto L1;
	}
	else
	{
	L2:;
		CDT = 2.0e0 * d_rand2_(4.0e0) - 1.0e0;
		if (d_rand2_(5.0e0) > 1.0e0 - CDT * CDT)
			goto L2;
	}
	CDT = (CDT + BETAP) / (1.0e0 + BETAP * CDT);
}

__device__ void d_esia2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IZZ, int &ISH)
{

	/*
		Amostragem aleatória de ionização de camada interna por impacto de elétrons.

		Modelo Sternheimer-Liljequist GOS.

		Argumentos de entrada:
		E ....... energia do elétron (eV).
		M ....... material onde os elétrons se propagam.
		DELTA ... Correção do efeito de densidade de Fermi.
		Argumentos de saída:
		DE ...... perda de energia (eV).
		EP ...... energia do elétron espalhado (eV).
		CDT ..... cosseno do ângulo de dispersão polar.
		ES ...... energia do elétron secundário emitido (eV).
		CDTS .... cosseno polar de direção do elétron secundário.
		IZZ ..... número atômico do elemento onde a ionização ocorreu.
		ISH ..... camada de elétrons atômicos que foi ionizada.
	*/

	/*if (imprimiu==0){
		printf("\n\nESIA2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double PK, TST, UK, WK, WTHR, WM, WKP, QKP, EE, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV;
	int IO, JO, IT, IOSC, JE;

	// Ponto da rede de energia

	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;
	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Seleção do oscilador ativo.

	TST = d_rand2_(2.0e0);

	// Busca binaria
	IO = 1;
	JO = dg_CESIAC_.NESI[M - 1] + 1;
L1:;
	IT = (IO + JO) / 2;
	if (TST > dg_CESIAC_.ESIAC[IT - 1][JE - 1][M - 1])
	{
		IO = IT;
	}
	else
	{
		JO = IT;
	}
	if ((JO - IO) > 1)
		goto L1;

	IOSC = dg_CESIAC_.IESI[IO - 1][M - 1];
	IZZ = dg_CEIN_.KZ[IOSC - 1][M - 1];
	ISH = dg_CEIN_.KS[IOSC - 1][M - 1];
	UK = dg_CEIN_.UI[IOSC - 1][M - 1];
	WK = dg_CEIN_.WRI[IOSC - 1][M - 1];

	if (UK > 1.0e-3)
	{
		WTHR = UK;
	}
	else
	{
		WTHR = WK;
	}

	if (E < (WTHR + 1.0e-6))
	{
		DE = UK;
		EP = E - DE;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	 conchas internas são variadas para produzir um limiar suave.
	*/

	WM = 3.0e0 * WK - 2.0e0 * UK;
	if (E > WM)
	{
		WKP = WK;
		QKP = UK;
	}
	else
	{
		WKP = (E + 2.0e0 * UK) / 3.0e0;
		QKP = UK * (E / WM);
		WM = E;
	}
	EE = E + UK;
	WCMAX = 0.5e0 * EE;
	WDMAX = fmin(WM, WCMAX);

	// Constantes

	RB = E + d_TREV;
	GAM = 1.0e0 + E * d_RREV;
	GAM2 = GAM * GAM;
	BETA2 = (GAM2 - 1.0e0) / GAM2;
	AMOL = pow(((GAM - 1.0e0) / GAM), 2);
	CPS = E * RB;
	CP = sqrt(CPS);

	// Seções transversais parciais do oscilador ativo.

	// Excitações distantes.

	CPPS = (E - WKP) * (E - WKP + d_TREV);
	CPP = sqrt(CPPS);

	if (WKP > (1.0e-6 * E))
	{
		QM = sqrt(pow((CP - CPP), 2) + d_REV * d_REV) - d_REV;
	}
	else
	{
		QM = pow(WKP, 2) / (BETA2 * d_TREV);
		QM = QM * (1.0e0 - QM * d_RTREV);
	}

	if (QM < QKP)
	{
		RWKP = 1.0e0 / WKP;
		XHDL = log(QKP * (QM + d_TREV) / (QM * (QKP + d_TREV))) * RWKP;
		XHDT = fmax(log(GAM2) - BETA2 - DELTA, 0.0e0) * RWKP;
		F0 = (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR) / pow((WM - UK), 2);
		XHDL = F0 * XHDL;
		XHDT = F0 * XHDT;
	}
	else
	{
		XHDL = 0.0e0;
		XHDT = 0.0e0;
	}

	// Colisoes Proximas
	RCL = WTHR / EE;
	RL1 = 1.0e0 - RCL;
	RRL1 = 1.0e0 / RL1;
	XHC = (AMOL * (0.5e0 - RCL) + 1.0e0 / RCL - RRL1 + (1.0e0 - AMOL) * log(RCL * RRL1)) / EE;

	XHTOT = XHC + XHDL + XHDT;

	if (XHTOT < 1.0e-35)
	{
		DE = UK;
		EP = E - DE;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		return;
	}

	// Amostragem de variáveis ​​de estado final.

	TST = d_rand2_(3.0e0) * XHTOT;

	// Colisão aproximada dura.

	TS1 = XHC;
	if (TST < TS1)
	{
		A = 5.0e0 * AMOL;
		ARCL = A * 0.5e0 * RCL;
	L2:;
		FB = (1.0e0 + ARCL) * d_rand2_(4.0e0);
		if (FB < 1.0e0)
		{
			RK = RCL / (1.0e0 - FB * (1.0e0 - (RCL + RCL)));
		}
		else
		{
			RK = RCL + (FB - 1.0e0) * (0.5e0 - RCL) / ARCL;
		}
		RK2 = RK * RK;
		RKF = RK / (1.0e0 - RK);
		PHI = 1.0e0 + pow(RKF, 2) - RKF + AMOL * (RK2 + RKF);
		if (d_rand2_(5.0e0) * (1.0e0 + A * RK2) > PHI)
			goto L2;
		// Energia e ângulo de espalhamento (elétron primário).
		DE = RK * EE;
		EP = E - DE;
		CDT = sqrt(EP * RB / (E * (RB - DE)));
		// Energia e ângulo de emissão do raio delta.
		ES = DE - UK;
		CDTS = sqrt(DE * RB / (E * (DE + d_TREV)));
		return;
	}

	// Interação longitudinal dura distante.
	TS1 = TS1 + XHDL;
	DE = WM - sqrt(pow((WM - WTHR), 2) - d_rand2_(7.0e0) * (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR));
	EP = E - DE;

	if (TST < TS1)
	{
		QS = QM / (1.0e0 + QM * d_RTREV);
		Q = QS / (pow(((QS / QKP) * (1.0e0 + QKP * d_RTREV)), d_rand2_(6.0e0)) - (QS * d_RTREV));
		QTREV = Q * (Q + d_TREV);
		CDT = (CPPS + CPS - QTREV) / (2.0e0 * CP * CPP);
		if (CDT > 1.0e0)
			CDT = 1.0e0;

		// Energia e ângulo de emissão do raio delta.
		ES = DE - UK;
		CDTS = 0.5e0 * (WKP * (E + RB - WKP) + QTREV) / sqrt(CPS * QTREV);
		if (CDTS > 1.0e0)
			CDTS = 1.0e0;
		return;
	}

	// Interação transversal distante difícil.

	CDT = 1.0e0;

	// Energia e ângulo de emissão do raio delta.

	ES = DE - UK;
	CDTS = 1.0e0;
}

__device__ void d_relax2_(int &IZ, int &IS)
{
	/*
	Esta sub-rotina simula o relaxamento de um átomo ionizado
	o elemento IZ com uma vaga no shell IS (o shell K ou um L
	subcamada C ou M). Essa vacância inicial é preenchida por elétrons de
	camadas externas através de transições radiativas e não radiativas, que
	pode produzir vagas adicionais.

	Usamos a seguinte notação para designar as transições possíveis:
	* Radiativo: IS0-IS1 (um elétron da camada IS1 preenche o
	vacância na casca IS0, deixando um buraco na casca IS1).
	* Não radiativo: IS0-IS1-IS2 (um elétron da camada IS1 preenche
	a vacância na camada IS0, e a energia liberada é tomada
	afastado por um elétron na camada IS2; este processo deixa dois
	, nas carcaças IS1 e IS2).
	A cascata de desexcitação (ou seja, o conjunto de transições que ocorrem para
	uma determinada vaga inicial) é amostrado a partir das probabilidades de transição
	contido na Biblioteca de Dados Atômicos Avaliados de Livermore (EADL). O
	A energia C da radiação emitida em cada transição é lida do
	Base de dados C PENELOPE.

	A simulação da cascata de desexcitação é descontinuada
	quando as conchas K, L e M estiverem cheias ou quando não houver
	energia suficiente para produzir radiação 'ativa' (com energia maior que
	EAB). A energia de excitação do íon residual é assumida como sendo
	depositado localmente. Desconsideramos a emissão e transporte de soft
	raios-x e elétrons lentos, cujas energias são menores do que a ligação
	energia C da camada N1 do elemento mais pesado no meio. este
	estabelece um limite inferior para o intervalo de energia que pode ser coberto pelo
	programa de simulação C de forma consistente.

	Os dados de desexcitação para os elementos carregados são armazenados no
	Bloco /CRELAX/, em uma forma projetada para minimizar a quantidade de memória
	e para facilitar a amostragem aleatória. As quantidades em comum
	bloco são os seguintes:
	IFIRST(99,16) ... dados de desexcitação para uma vaga no shell IS de
	o elemento IZ começa na posição K=IFIRST(IZ,IS) no
	matrizes de armazenamento. Os valores permitidos para IS são de 1 a 16 (K shell
	subcamadas L, M e N).
	ILAST(99,16) ... os dados de desexcitação para uma vaga no shell
	IS do elemento IZ termina na posição K=ILAST(IZ,IS) no
	matrizes de armazenamento.
	IS1(K), IS2(K) ... shells que estão ativos na transição (veja o
	código do rótulo do shell abaixo). Para transições radiativas, IS2(K)=0.
	P(K) ... probabilidade relativa para a transição IS-IS1(K)-IS2(K).
	ET(K) ... energia da partícula secundária emitida na transição.
	F(K), IAL(K) ... valores de corte e alias (método de amostragem de Walker).

  ---------------------------------------------------------------------
  Label code IS for electron shells:
	  1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
	  2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
	  3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
	  4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
	  5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
	  6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
	  7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
	  8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
	  9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
	 10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  ---------------------------------------------------------------------
	*/

	/*	if (imprimiu==0){
		 printf("\n\nRELAX2\n\n");
		 imprimiu++;
	 }*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double PTIM[256];
	int ISV[256];

	int NV, ISP, KF, KL, K1, NVIS, IS1K, IS2K, KPARS;
	double PAGE0, RN, TST, WS, US, VS, DF, SDTS;

	// Inicializacao

	if ((IZ < 3) || (IS > 16))
	{
		return;
	}

	if (dg_CADATA_.EB[IS - 1][IZ - 1] < dg_CECUTR_.ECUTR[dg_TRACK_mod_[index].MAT - 1]) ////Se a energia de ionização da casca for menor que ECUTR, a cascata não é seguido.
		return;

	NV = 1;
	ISV[1 - 1] = IS;
	PAGE0 = dg_TRACK_mod_[index].PAGE;
	PTIM[1 - 1] = dg_TRACK_mod_[index].PAGE;

	// Proxima transição

L1:;
	ISP = ISV[NV - 1];
	dg_TRACK_mod_[index].PAGE = PTIM[NV - 1];
	KF = dg_CRELAX_.IFIRST[ISP - 1][IZ - 1];
	KL = dg_CRELAX_.ILAST[ISP - 1][IZ - 1];
	NV = NV - 1;

	if (KL > KF)
	{
		// Algoritmo de amostragem de Walker.
		RN = d_rand2_(1.0e0) * (KL - KF + 1);
		K1 = int(RN);
		TST = RN - K1;
		if (TST > dg_CRELAX_.F[KF + K1 - 1])
		{
			dg_CRELAX_.KS[trackIndex] = dg_CRELAX_.IS0[KF + K1 - 1];
		}
		else
		{
			dg_CRELAX_.KS[trackIndex] = KF + K1;
		}
	}
	else
	{
		dg_CRELAX_.KS[trackIndex] = KF;
	}

	/*
	Se MODER=0, o controle é retornado ao programa chamador após
	determinando a primeira transição, KS. Útil para testar o aleatório
	amostragem. Para operação normal, podemos comentar o seguinte
	declaração .
	*/

	if (dg_CRELAX_.MODER == 0)
		return;
	// Se LAGE=true, a idade da partícula é registrada.
	if (dg_TRACK_mod_[index].LAGE)
	{
		NVIS = 1;
		if (NV > 1)
		{ // Várias vagas no shell ativo?
			for (int ISC = 1; ISC <= NV; ISC++)
			{
				if (ISV[ISC - 1] == ISP)
					NVIS = NVIS + 1;
			}
		}
		dg_TRACK_mod_[index].PAGE = dg_TRACK_mod_[index].PAGE - (dg_CADATA_.ALW[ISP - 1][IZ - 1] / (NVIS)) * log(d_rand2_(2.0e0));
	}

	// radiação fluorescente

	IS1K = dg_CRELAX_.IS1[dg_CRELAX_.KS[trackIndex] - 1];
	IS2K = dg_CRELAX_.IS2[dg_CRELAX_.KS[trackIndex] - 1];
	if (IS2K == 0)
	{
		KPARS = 2;
		if (IS1K < 17)
		{
			if (dg_CADATA_.EB[IS1K - 1][IZ - 1] > dg_CECUTR_.ECUTR[dg_TRACK_mod_[index].MAT - 1])
			{
				NV = NV + 1;
				ISV[NV - 1] = IS1K;
				PTIM[NV - 1] = dg_TRACK_mod_[index].PAGE;
			}
		}
	}
	else
	{
		KPARS = 1;
		if (IS1K < 17)
		{
			if (dg_CADATA_.EB[IS1K - 1][IZ - 1] > dg_CECUTR_.ECUTR[dg_TRACK_mod_[index].MAT - 1])
			{
				NV = NV + 1;
				ISV[NV - 1] = IS1K;
				PTIM[NV - 1] = dg_TRACK_mod_[index].PAGE;
			}
		}
		if (IS2K < 17)
		{
			if (dg_CADATA_.EB[IS2K - 1][IZ - 1] > dg_CECUTR_.ECUTR[dg_TRACK_mod_[index].MAT - 1])
			{
				NV = NV + 1;
				ISV[NV - 1] = IS2K;
				PTIM[NV - 1] = dg_TRACK_mod_[index].PAGE;
			}
		}
	}

	/*
	A partícula emitida é armazenada na pilha secundária quando
	sua energia ET(K) é maior que EABS.
	*/

	if (dg_CRELAX_.ET[dg_CRELAX_.KS[trackIndex] - 1] > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][KPARS - 1])
	{
		// Direção inicial (isotrópica).
		WS = -1.0e0 + 2.0e0 * d_rand2_(2.0e0);
		SDTS = sqrt(1.0e0 - WS * WS);
		DF = d_TWOPI * d_rand2_(3.0e0);
		US = cos(DF) * SDTS;
		VS = sin(DF) * SDTS;
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = IZ * 1000000 + ISP * 10000 + IS1K * 100 + IS2K;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		d_stores2_(dg_CRELAX_.ET[dg_CRELAX_.KS[trackIndex] - 1], dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, KPARS, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}

	// Existem vagas não preenchidas nas conchas internas?

	if (NV > 0)
		goto L1;
	dg_TRACK_mod_[index].PAGE = PAGE0;
}

__device__ void d_stores2_(double &EI, double &XI, double &YI, double &ZI, double &UI, double &VI, double &WI, double &WGHTI, int &KPARI, int *ILBI, int &IPOLI)
{

	/*
		Esta sub-rotina armazena o estado inicial de uma nova partícula secundária
		na pilha secundária. Os valores de entrada são:
		EI ........... energia inicial.
		XI, YI, ZI ... coordenadas da posição inicial.
		UI, VI, WI ... cossenos de direção inicial.
		WGHTI ........ peso (=1 na simulação analógica).
		KPARI ........ tipo de partícula (1: elétron, 2: fóton,
		3: pósitron).
		ILBI(5) ...... rótulos de partículas.
		IPOLI ........ sinalizador de polarização.

		O parâmetro NMS fixa o tamanho da pilha secundária (ou seja,
		número máximo de partículas que podem ser armazenadas). Se este número for
		excedido, uma mensagem de aviso é impressa na unidade 26. Quando a memória
		O armazenamento de é esgotado, cada nova partícula secundária é armazenada no
		posição do elétron secundário menos energético ou fóton já
		produzido, que é assim descartado.
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double EM = 1.0e35;
	int IE, IG, IP;

	//atomicAdd2(&d_contStores, 1);
	//printf("Contagem de Stores %d", d_contStores);

	if (KPARI == 1)
	{
		if (dg_nTRACKS_.nSECTRACK_E == pilhaPart * 100)
		{
			// printf("pilha eletrons estourou\n");
			for (int I = 0; I < pilhaPart * 100; I++)
			{
				if (dg_SECTRACK_E_[I].E < EM)
				{
					EM = dg_SECTRACK_E_[I].E;
					IE = I;
				}
			}
		}
		else
		{
			IE = atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, 1);
			//atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, 1);

			//printf("Index nSecTrack_E %d", dg_nTRACKS_.nSECTRACK_E);

			//IE = dg_nTRACKS_.nSECTRACK_E;
			/*if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, 1);
			}
			else
			{
				dg_nTRACKS_.nSECTRACK_E = dg_nTRACKS_.nSECTRACK_E + 1;
			}*/
		}

		if (EI == 0.0e0){
			printf("Energia ZERADA Index %d , ICOL = %d \n\n\n", index, ILBI[3 - 1]);
		}

		dg_SECTRACK_E_[IE].E = EI;
		dg_SECTRACK_E_[IE].X = XI;
		dg_SECTRACK_E_[IE].Y = YI;
		dg_SECTRACK_E_[IE].Z = ZI;
		dg_SECTRACK_E_[IE].U = UI;
		dg_SECTRACK_E_[IE].V = VI;
		dg_SECTRACK_E_[IE].W = WI;
		dg_SECTRACK_E_[IE].WGHT = WGHTI;
		dg_SECTRACK_E_[IE].KPAR = KPARI;
		dg_SECTRACK_E_[IE].IBODY = dg_TRACK_mod_[index].IBODY;
		dg_SECTRACK_E_[IE].MAT = dg_TRACK_mod_[index].MAT;

		dg_SECTRACK_E_[IE].ILB[1 - 1] = ILBI[1 - 1];
		dg_SECTRACK_E_[IE].ILB[2 - 1] = ILBI[2 - 1];
		dg_SECTRACK_E_[IE].ILB[3 - 1] = ILBI[3 - 1];
		dg_SECTRACK_E_[IE].ILB[4 - 1] = ILBI[4 - 1];
		dg_SECTRACK_E_[IE].ILB[5 - 1] = ILBI[5 - 1];

		dg_SECTRACK_E_[IE].N = dg_TRACK_mod_[index].N;
		dg_SECTRACK_E_[IE].INDEX = dg_TRACK_mod_[index].INDEX;
		dg_SECTRACK_E_[IE].IEXIT = 0;
		dg_SECTRACK_E_[IE].STEP = 11;
		dg_SECTRACK_E_[IE].CROSS = false;
		


		if (IPOLI == 1)
		{
			dg_SECTRACK_E_[IE].SP1 = dg_TRACK_mod_[index].SP1;
			dg_SECTRACK_E_[IE].SP2 = dg_TRACK_mod_[index].SP2;
			dg_SECTRACK_E_[IE].SP3 = dg_TRACK_mod_[index].SP3;
			dg_SECTRACK_E_[IE].IPOL = IPOLI;
		}
		else
		{
			dg_SECTRACK_E_[IE].SP1 = 0.0e0;
			dg_SECTRACK_E_[IE].SP2 = 0.0e0;
			dg_SECTRACK_E_[IE].SP3 = 0.0e0;
			dg_SECTRACK_E_[IE].IPOL = 0;
		}
		dg_SECTRACK_E_[IE].PAGE = dg_TRACK_mod_[index].PAGE;
	}
	else if (KPARI == 2)
	{
		if (dg_nTRACKS_.nSECTRACK_G == pilhaPart)
		{
			// printf("pilha fotons estourou\n");
			for (int I = 0; I < pilhaPart; I++)
			{
				if (dg_SECTRACK_G_[I].E < EM)
				{
					EM = dg_SECTRACK_G_[I].E;
					IG = I;
				}
			}
		}
		else
		{
			IG = dg_nTRACKS_.nSECTRACK_G;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nSECTRACK_G, 1);
			}
			else
			{
				dg_nTRACKS_.nSECTRACK_G = dg_nTRACKS_.nSECTRACK_G + 1;
			}
		}

		dg_SECTRACK_G_[IG].E = EI;
		dg_SECTRACK_G_[IG].X = XI;
		dg_SECTRACK_G_[IG].Y = YI;
		dg_SECTRACK_G_[IG].Z = ZI;
		dg_SECTRACK_G_[IG].U = UI;
		dg_SECTRACK_G_[IG].V = VI;
		dg_SECTRACK_G_[IG].W = WI;
		dg_SECTRACK_G_[IG].WGHT = WGHTI;
		dg_SECTRACK_G_[IG].KPAR = KPARI;
		dg_SECTRACK_G_[IG].IBODY = dg_TRACK_mod_[index].IBODY;
		dg_SECTRACK_G_[IG].MAT = dg_TRACK_mod_[index].MAT;
		dg_SECTRACK_G_[IG].ILB[1 - 1] = ILBI[1 - 1];
		dg_SECTRACK_G_[IG].ILB[2 - 1] = ILBI[2 - 1];
		dg_SECTRACK_G_[IG].ILB[3 - 1] = ILBI[3 - 1];
		dg_SECTRACK_G_[IG].ILB[4 - 1] = ILBI[4 - 1];
		dg_SECTRACK_G_[IG].ILB[5 - 1] = ILBI[5 - 1];

		dg_SECTRACK_G_[IG].N = dg_TRACK_mod_[index].N;
		dg_SECTRACK_G_[IG].INDEX = dg_TRACK_mod_[index].INDEX;
		dg_SECTRACK_G_[IG].IEXIT = dg_TRACK_mod_[index].IEXIT;
		dg_SECTRACK_E_[IG].STEP = 11;
		dg_SECTRACK_E_[IE].CROSS = false;

		if (IPOLI == 1)
		{
			dg_SECTRACK_G_[IG].SP1 = dg_TRACK_mod_[index].SP1;
			dg_SECTRACK_G_[IG].SP2 = dg_TRACK_mod_[index].SP2;
			dg_SECTRACK_G_[IG].SP3 = dg_TRACK_mod_[index].SP3;
			dg_SECTRACK_G_[IG].IPOL = IPOLI;
		}
		else
		{
			dg_SECTRACK_G_[IG].SP1 = 0.0e0;
			dg_SECTRACK_G_[IG].SP2 = 0.0e0;
			dg_SECTRACK_G_[IG].SP3 = 0.0e0;
			dg_SECTRACK_G_[IG].IPOL = 0;
		}
		dg_SECTRACK_G_[IG].PAGE = dg_TRACK_mod_[index].PAGE;
	}
	else if (KPARI == 3)
	{
		if (dg_nTRACKS_.nSECTRACK_P == pilhaPart)
		{
			// printf("pilha positrons estourou\n");
			for (int I = 0; I < pilhaPart; I++)
			{
				if (dg_SECTRACK_P_[I].E < EM)
				{
					EM = dg_SECTRACK_P_[I].E;
					IP = I;
				}
			}
		}
		else
		{
			IP = dg_nTRACKS_.nSECTRACK_P;
			if (d_cudaCapability >= 6)
			{
				atomicAdd2(&dg_nTRACKS_.nSECTRACK_P, 1);
			}
			else
			{
				dg_nTRACKS_.nSECTRACK_P = dg_nTRACKS_.nSECTRACK_P + 1;
			}
		}

		dg_SECTRACK_P_[IP].E = EI;
		dg_SECTRACK_P_[IP].X = XI;
		dg_SECTRACK_P_[IP].Y = YI;
		dg_SECTRACK_P_[IP].Z = ZI;
		dg_SECTRACK_P_[IP].U = UI;
		dg_SECTRACK_P_[IP].V = VI;
		dg_SECTRACK_P_[IP].W = WI;
		dg_SECTRACK_P_[IP].WGHT = WGHTI;
		dg_SECTRACK_P_[IP].KPAR = KPARI;
		dg_SECTRACK_P_[IP].IBODY = dg_TRACK_mod_[index].IBODY;
		dg_SECTRACK_P_[IP].MAT = dg_TRACK_mod_[index].MAT;
		dg_SECTRACK_P_[IP].ILB[1 - 1] = ILBI[1 - 1];
		dg_SECTRACK_P_[IP].ILB[2 - 1] = ILBI[2 - 1];
		dg_SECTRACK_P_[IP].ILB[3 - 1] = ILBI[3 - 1];
		dg_SECTRACK_P_[IP].ILB[4 - 1] = ILBI[4 - 1];
		dg_SECTRACK_P_[IP].ILB[5 - 1] = ILBI[5 - 1];

		dg_SECTRACK_P_[IP].N = dg_TRACK_mod_[index].N;
		dg_SECTRACK_P_[IP].INDEX = dg_TRACK_mod_[index].INDEX;
		dg_SECTRACK_P_[IP].IEXIT = dg_TRACK_mod_[index].IEXIT;
		dg_SECTRACK_E_[IP].STEP = 11;
		dg_SECTRACK_E_[IE].CROSS = false;


		if (IPOLI == 1)
		{
			dg_SECTRACK_P_[IP].SP1 = dg_TRACK_mod_[index].SP1;
			dg_SECTRACK_P_[IP].SP2 = dg_TRACK_mod_[index].SP2;
			dg_SECTRACK_P_[IP].SP3 = dg_TRACK_mod_[index].SP3;
			dg_SECTRACK_P_[IP].IPOL = IPOLI;
		}
		else
		{
			dg_SECTRACK_P_[IP].SP1 = 0.e0;
			dg_SECTRACK_P_[IP].SP2 = 0.e0;
			dg_SECTRACK_P_[IP].SP3 = 0.e0;
			dg_SECTRACK_P_[IP].IPOL = 0;
		}
		dg_SECTRACK_P_[IP].PAGE = dg_TRACK_mod_[index].PAGE;
	}
}

//__device__ void d_stores2_(double& EI, double& XI, double& YI, double& ZI, double& UI, double& VI, double& WI, double& WGHTI, int& KPARI, int* ILBI, int& IPOLI) {

/*
	Esta sub-rotina armazena o estado inicial de uma nova partícula secundária
	na pilha secundária. Os valores de entrada são:
	EI ........... energia inicial.
	XI, YI, ZI ... coordenadas da posição inicial.
	UI, VI, WI ... cossenos de direção inicial.
	WGHTI ........ peso (=1 na simulação analógica).
	KPARI ........ tipo de partícula (1: elétron, 2: fóton,
	3: pósitron).
	ILBI(5) ...... rótulos de partículas.
	IPOLI ........ sinalizador de polarização.

	O parâmetro NMS fixa o tamanho da pilha secundária (ou seja,
	número máximo de partículas que podem ser armazenadas). Se este número for
	excedido, uma mensagem de aviso é impressa na unidade 26. Quando a memória
	O armazenamento de é esgotado, cada nova partícula secundária é armazenada no
	posição do elétron secundário menos energético ou fóton já
	produzido, que é assim descartado.
*/
/*	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double EME, EMG;
	int IE, IG, IS;

	if (dg_SECST_.NSEC < NMS) {
		dg_SECST_.NSEC = dg_SECST_.NSEC + 1;
		IS = dg_SECST_.NSEC;
	}
	else {
		if (dg_CERSEC_.IERSEC == 0) {
			printf("   *** WARNING: (STORES) not enough storage for secondaries.\n EABS(KPAR,MAT) or the parameter NMS should be enlarged\n");
			dg_CERSEC_.IERSEC = 1;
		}
		dg_SECST_.NSEC = NMS;
		EME = 1.0e35;
		EMG = 1.0e35;
		IE = 0;
		IG = 0;
		for (int I = 1; I <= NMS; I++) {
			if (dg_SECST_.KS[I - 1] == 1) {
				if (dg_SECST_.ES[I - 1] < EME) {
					EME = dg_SECST_.ES[I - 1];
					IE = I;
				}
			}
			else if (dg_SECST_.KS[I - 1] == 2) {
				if (dg_SECST_.ES[I - 1] < EMG) {
					EMG = dg_SECST_.ES[I - 1];
					IG = I;
				}
			}
		}

		if (IE > 0) {
			IS = IE;
		}
		else if (IG > 0) {
			IS = IG;
		}
		else {
			printf("   *** Not enough storage for secondary positrons./n       JOB ABORTED.");
			IS = 0;
			return;
			//exit(1);
			//assert(0);
		}
	}

	dg_SECST_.ES[IS - 1] = EI;
	dg_SECST_.XS[IS - 1] = XI;
	dg_SECST_.YS[IS - 1] = YI;
	dg_SECST_.ZS[IS - 1] = ZI;
	dg_SECST_.US[IS - 1] = UI;
	dg_SECST_.VS[IS - 1] = VI;
	dg_SECST_.WS[IS - 1] = WI;
	dg_SECST_.WGHTS[IS - 1] = WGHTI;
	dg_SECST_.KS[IS - 1] = KPARI;
	dg_SECST_.IBODYS[IS - 1] = dg_TRACK_mod_[index].IBODY;
	dg_SECST_.MS[IS - 1] = dg_TRACK_mod_[index].MAT;
	dg_SECST_.ILBS[IS - 1][1 - 1] = ILBI[1 - 1];
	dg_SECST_.ILBS[IS - 1][2 - 1] = ILBI[2 - 1];
	dg_SECST_.ILBS[IS - 1][3 - 1] = ILBI[3 - 1];
	dg_SECST_.ILBS[IS - 1][4 - 1] = ILBI[4 - 1];
	dg_SECST_.ILBS[IS - 1][5 - 1] = ILBI[5 - 1];

	if (IPOLI == 1) {
		dg_SECST_.SP1S[IS - 1] = dg_TRACK_mod_[index].SP1;
		dg_SECST_.SP2S[IS - 1] = dg_TRACK_mod_[index].SP2;
		dg_SECST_.SP3S[IS - 1] = dg_TRACK_mod_[index].SP3;
		dg_SECST_.IPOLS[IS - 1] = IPOLI;
	}
	else {
		dg_SECST_.SP1S[IS - 1] = 0.e0;ste
		dg_SECST_.SP2S[IS - 1] = 0.e0;
		dg_SECST_.SP3S[IS - 1] = 0.e0;
		dg_SECST_.IPOLS[IS - 1] = 0;
	}

	dg_SECST_.PAGES[IS - 1] = dg_TRACK_mod_[index].PAGE;


}*/

__device__ void d_eeld2_(double &RNDC, double &RMU)
{

	/*
	Simulação de eventos elásticos rígidos de elétrons. Seções transversais de
	a base de dados numérica ELSEPA.

	Valor do argumento C:
	RNDC ... valor de corte do número aleatório uniforme
	(apenas eventos difíceis são simulados).
	RMU .... deflexão angular amostrada, =(1-CDT)/2.

	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	static const int NP = 128;
	static const int NPM1 = NP - 1;

	double PK, RU, RR, PP, XX, AA, BB, D;
	int ITN, JE, I, J, K;

	// Ponto da rede de energia

	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;

	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Ponto
	RU = RNDC + d_rand2_(2.0e0) * (1.0e0 - RNDC);

	// Selection of the interval (binary search in a restricted interval).
	ITN = (int)(RU * NPM1 + 1);
	I = dg_CEELDB_.ITLE[dg_TRACK_mod_[index].MAT - 1][JE - 1][ITN - 1];
	J = dg_CEELDB_.ITUE[dg_TRACK_mod_[index].MAT - 1][JE - 1][ITN - 1];
	if ((J - I) < 2)
		goto L2;
L1:;
	K = (I + J) / 2;
	if (RU > dg_CEELDB_.PSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][K - 1])
	{
		I = K;
	}
	else
	{
		J = K;
	}

	if ((J - I) > 1)
		goto L1;

	// Amostragem da distribuição cumulativa inversa racional.
L2:;
	PP = dg_CEELDB_.PSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
	RR = RU - PP;
	if (RR > 1.0e-16)
	{
		XX = dg_CEELDB_.XSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		AA = dg_CEELDB_.ASE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		BB = dg_CEELDB_.BSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
		D = dg_CEELDB_.PSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I + 1 - 1] - PP;
		RMU = XX + ((1.0e0 + AA + BB) * D * RR / (D * D + (AA * D + BB * RR) * RR)) * (dg_CEELDB_.XSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I + 1 - 1] - XX);
	}
	else
	{
		RMU = dg_CEELDB_.XSE[dg_TRACK_mod_[index].MAT - 1][JE - 1][I - 1];
	}
}

__device__ void d_eina2_(double &E, double &DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int &M, int &IOSC)
{

	/*
	Amostragem aleatória de colisões inelásticas duras de elétrons.

	Modelo C Sternheimer-Liljequist GOS.

	Argumentos de entrada:
	E ....... energia do elétron (eV).
	M ....... material onde os elétrons se propagam.
	DELTA ... Correção do efeito de densidade de Fermi.
	Argumentos de saída:
	DE ...... perda de energia (eV).
	EP ...... energia do elétron espalhado (eV).
	CDT ..... cosseno do ângulo de dispersão polar.
	ES ...... energia do elétron secundário emitido (eV).
	CDTS .... cosseno polar de direção do elétron secundário.
	IOSC .... índice do oscilador que foi 'ionizado'.

	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	bool LDIST;

	static const int NO = 512;

	double WCCM, PK, TST, UK, WK, WTHR, WM, WKP, QKP, EE, WCMAX, WDMAX, RB, GAM, GAM2, BETA2, AMOL, CPS, CP;
	double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0, RCL, RL1, RRL1, XHC, XHTOT, TS1, A, ARCL, FB, RK, RK2, RKF, PHI;
	double QS, Q, QTREV;
	int IO, JO, IT, JE;

	WCCM = dg_PENELOPE_mod_.WCC[M - 1];

	if (WCCM > E)
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}
	// Ponto da rede de energia
	PK = (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP[dg_CEGRID_.KE[trackIndex] - 1]) * dg_CEGRID_.DLFC;
	if (d_rand2_(1.0e0) < PK)
	{
		JE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		JE = dg_CEGRID_.KE[trackIndex];
	}

	// Seleção do oscilador ativo.
	TST = d_rand2_(2.0e0);
	// Busca binaria
	IO = 1;
	JO = dg_CEINAC_.NEIN[M - 1] + 1;
L1:;
	IT = (IO + JO) / 2;
	if (TST > dg_CEINAC_.EINAC[IT - 1][JE - 1][M - 1])
	{
		IO = IT;
	}
	else
	{
		JO = IT;
	}
	if (JO - IO > 1)
		goto L1;
	IOSC = dg_CEINAC_.IEIN[IO - 1][M - 1];
	UK = dg_CEIN_.UI[IOSC - 1][M - 1];
	WK = dg_CEIN_.WRI[IOSC - 1][M - 1];

	if (UK > 1.0e-3)
	{
		WTHR = fmax(WCCM, UK);
	}
	else
	{
		WTHR = fmax(WCCM, WK);
	}

	if (E < WTHR + 1.0e-6)
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}

	/*
	Truque: A energia de ressonância e a energia de recuo de corte de
	conchas internas são variadas para produzir um limiar suave.
	*/

	LDIST = true;
	if (UK > 1.0e-3)
	{
		WM = 3.0e0 * WK - 2.0e0 * UK;
		if (E > WM)
		{
			WKP = WK;
			QKP = UK;
		}
		else
		{
			WKP = (E + 2.0e0 * UK) / 3.0e0;
			QKP = UK * (E / WM);
			WM = E;
		}
		if (WCCM > WM)
			LDIST = false;

		EE = E + UK;
		WCMAX = 0.5e0 * EE;
		WDMAX = fmin(WM, WCMAX);
		if (WTHR > WDMAX)
			LDIST = false;
	}
	else
	{
		if (WCCM > WK)
			LDIST = false;
		WKP = WK;
		QKP = WK;
		WM = E;
		EE = E;
		WCMAX = 0.5e0 * EE;
		WDMAX = WKP + 1.0e0;
	}

	// Constantes

	RB = E + d_TREV;
	GAM = 1.0e0 + E * d_RREV;
	GAM2 = GAM * GAM;
	BETA2 = (GAM2 - 1.0e0) / GAM2;
	AMOL = pow(((GAM - 1.0e0) / GAM), 2);
	CPS = E * RB;
	CP = sqrt(CPS);

	// Seções transversais parciais do oscilador ativo.
	// Excitalçoes Distantes
	if (LDIST)
	{
		CPPS = (E - WKP) * (E - WKP + d_TREV);
		CPP = sqrt(CPPS);
		if (WKP > 1.0e-6 * E)
		{
			QM = sqrt(pow((CP - CPP), 2) + d_REV * d_REV) - d_REV;
		}
		else
		{
			QM = pow(WKP, 2) / (BETA2 * d_TREV);
			QM = QM * (1.0e0 - QM * d_RTREV);
		}
		if (QM < QKP)
		{
			RWKP = 1.0e0 / WKP;
			XHDL = log(QKP * (QM + d_TREV) / (QM * (QKP + d_TREV))) * RWKP;
			XHDT = fmax(log(GAM2) - BETA2 - DELTA, 0.0e0) * RWKP;
			if (UK > 1.0e-3)
			{
				F0 = (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR) / pow((WM - UK), 2);
				XHDL = F0 * XHDL;
				XHDT = F0 * XHDT;
			}
		}
		else
		{
			XHDL = 0.0e0;
			XHDT = 0.0e0;
		}
	}
	else
	{
		QM = 0.0e0; // Definido para evitar avisos de compilação.
		CPP = 0.0e0;
		CPPS = 0.0e0;
		XHDL = 0.0e0;
		XHDT = 0.0e0;
	}

	// Colisoes Fechadas
	RCL = WTHR / EE;

	if (RCL < 0.5e0)
	{
		RL1 = 1.0e0 - RCL;
		RRL1 = 1.0e0 / RL1;
		XHC = (AMOL * (0.5e0 - RCL) + 1.0e0 / RCL - RRL1 + (1.0e0 - AMOL) * log(RCL * RRL1)) / EE;
	}
	else
	{
		XHC = 0.0e0;
	}

	XHTOT = XHC + XHDL + XHDT;
	if (XHTOT < 1.0e-35)
	{
		DE = 0.0e0;
		EP = E;
		CDT = 1.0e0;
		ES = 0.0e0;
		CDTS = 0.0e0;
		IOSC = NO;
		return;
	}

	// Amostragem de variáveis ​​de estado final.

	TST = d_rand2_(3.0e0) * XHTOT;

	// Colisão fechada dura

	TS1 = XHC;
	if (TST < TS1)
	{
		A = 5.0e0 * AMOL;
		ARCL = A * 0.5e0 * RCL;
	L2:;
		FB = (1.0e0 + ARCL) * d_rand2_(4.0e0);
		if (FB < 1.0e0)
		{
			RK = RCL / (1.0e0 - FB * (1.0e0 - (RCL + RCL)));
		}
		else
		{
			RK = RCL + (FB - 1.0e0) * (0.5e0 - RCL) / ARCL;
		}
		RK2 = RK * RK;
		RKF = RK / (1.0e0 - RK);
		PHI = 1.0e0 + pow(RKF, 2) - RKF + AMOL * (RK2 + RKF);
		if (d_rand2_(5.0e0) * (1.0e0 + A * RK2) > PHI)
			goto L2;
		// Energia e ângulo de espalhamento (elétron primário).
		DE = RK * EE;
		EP = E - DE;
		CDT = sqrt(EP * RB / (E * (RB - DE)));
		// Energia e ângulo de emissão do raio delta.
		if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
		{
			if (UK > dg_CECUTR_.ECUTR[M - 1])
			{
				ES = DE - UK;
			}
			else
			{
				ES = DE;
			}
		}
		else
		{
			ES = DE;
		}
		CDTS = sqrt(DE * RB / (E * (DE + d_TREV)));
		return;
	}

	// Interação longitudinal dura distante.
	TS1 = TS1 + XHDL;
	if (UK > 1.0e-3)
	{
		DE = WM - sqrt(pow((WM - WTHR), 2) - d_rand2_(7.0e0) * (WDMAX - WTHR) * (WM + WM - WDMAX - WTHR));
	}
	else
	{
		DE = WKP;
	}
	EP = E - DE;
	if (TST < TS1)
	{
		QS = QM / (1.0e0 + QM * d_RTREV);
		Q = QS / (pow(((QS / QKP) * (1.0e0 + QKP * d_RTREV)), d_rand2_(6.0e0)) - (QS * d_RTREV));
		QTREV = Q * (Q + d_TREV);
		CDT = (CPPS + CPS - QTREV) / (2.0e0 * CP * CPP);
		if (CDT > 1.0e0)
		{
			CDT = 1.0e0;
		}
		// Energia e ângulo de emissão do raio delta.
		if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
		{
			ES = DE - UK; // Apenas conchas internas.
		}
		else
		{
			ES = DE;
		}
		CDTS = 0.5e0 * (WKP * (E + RB - WKP) + QTREV) / sqrt(CPS * QTREV);
		if (CDTS > 1.0e0)
			CDTS = 1.0e0;
		return;
	}

	// Interação transversal distante difícil.
	CDT = 1.0e0;
	// Energia e ângulo de emissão do raio delta.
	if (dg_CEIN_.KS[IOSC - 1][M - 1] < 17)
	{
		if (UK > dg_CECUTR_.ECUTR[M - 1])
		{
			ES = DE - UK; // Apenas conchas internas.
		}
		else
		{
			ES = DE;
		}
	}
	else
	{
		ES = DE;
	}
	CDTS = 1.0e0;
}

__device__ void d_ebra2_(double &E, double &W, int &M)
{

	/*
	   Simulação da emissão de bremsstrahlung por elétrons ou pósitrons em
	   material M.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	static const int NBW = 32;
	int IE, I, J, K;
	double PT, W1, W2, DW, B, A, PMAX;

	if (dg_PENELOPE_mod_.WCR[M - 1] > E)
	{
		W = 0.0e0;
		return;
	}

	// Seleção do ponto da rede de energia.

	if (d_rand2_(1.0e0) < dg_CEGRID_.XEK[trackIndex])
	{
		IE = dg_CEGRID_.KE[trackIndex] + 1;
	}
	else
	{
		IE = dg_CEGRID_.KE[trackIndex];
	}
	// Ponteiro
L1:;
	PT = dg_CEBR_.PBCUT[IE - 1][M - 1] + d_rand2_(2.0e0) * (dg_CEBR_.PACB[NBW - 1][IE - 1][M - 1] - dg_CEBR_.PBCUT[IE - 1][M - 1]);

	// Pesquisa binária do intervalo W.
	I = 1;
	J = NBW;
L2:;
	K = (I + J) / 2;
	if (PT > dg_CEBR_.PACB[K - 1][IE - 1][M - 1])
	{
		I = K;
	}
	else
	{
		J = K;
	}

	if ((J - I) > 1)
	{
		//	printf("loop ebra2 gotoL2\n");
		//	printf("PT= %f, PACB= %f, I= %d, J= %d, K= %d\n\n", PT,  dg_CEBR_.PACB[K-1][IE-1][M-1], I, J, K );
		goto L2;
	}

	// Amostragem da energia do fóton (método de rejeição).
	W1 = dg_CEBR_.WB[I - 1];
	W2 = dg_CEBR_.WB[I + 1 - 1];
	DW = W2 - W1;
	B = dg_CEBR_.DPDFB[I - 1][IE - 1][M - 1] / DW;
	A = dg_CEBR_.PDFB[I - 1][IE - 1][M - 1] - B * W1;
	if (W1 < dg_CEBR_.WBCUT[IE - 1][M - 1])
	{
		W1 = dg_CEBR_.WBCUT[IE - 1][M - 1];
	}
	if (W2 < W1)
	{
		printf(" **** WARNING: EBR. Conflicting end-point values.\n");
		W = W1;
		return;
	}
	PMAX = fmax(A + B * W1, A + B * W2);
L3:;
	W = W1 * pow((W2 / W1), d_rand2_(3.0e0));
	if ((d_rand2_(4.0e0) * PMAX) > (A + B * W))
	{
		//	printf("loop ebra2 gotoL3\n");
		goto L3;
	}
	W = W * E;
	if (W < dg_PENELOPE_mod_.WCR[M - 1])
	{
		//	printf("loop ebra2 gotoL1\n");
		goto L1;
	}
}

__device__ void d_fimdet2_(int &N, int &ID, double &DSEF)
{

	/*
	Distribuição de fluência de partículas dentro do
	detector. Apenas colisões discretas.

	*/
	/*if (imprimiu==0){
		printf("\n\nFIMDET2\n\n");
		imprimiu++;
	}*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int IE;

	if (dg_CIMDET_.LLE[ID - 1] == 1)
	{
		IE = (int)(1.0e0 + (log(dg_TRACK_mod_[index].E) - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}
	else
	{
		IE = (int)(1.0e0 + (dg_TRACK_mod_[index].E - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}

	if ((IE > 0) && (IE <= dg_CIMDET_.NE[ID - 1]))
	{
		if (N != dg_CIMDET_.LFLT[IE - 1][ID - 1])
		{
			dg_CIMDET_.FLT[IE - 1][ID - 1] = dg_CIMDET_.FLT[IE - 1][ID - 1] + dg_CIMDET_.FLTP[IE - 1][ID - 1];
			dg_CIMDET_.FLT2[IE - 1][ID - 1] = dg_CIMDET_.FLT2[IE - 1][ID - 1] + pow(dg_CIMDET_.FLTP[IE - 1][ID - 1], 2);
			dg_CIMDET_.FLTP[IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT * DSEF;
			dg_CIMDET_.LFLT[IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.FLTP[IE - 1][ID - 1] = dg_CIMDET_.FLTP[IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT * DSEF;
		}

		if (N != dg_CIMDET_.LFLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1])
		{
			dg_CIMDET_.FLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1];
			dg_CIMDET_.FLP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + pow(dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1], 2);
			dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT * DSEF;
			dg_CIMDET_.LFLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT * DSEF;
		}
	}
}

__device__ void d_fimdes2_(int &N, int &ID, double &EI, double &DECSD, double &DSEF)
{

	/*
	Distribuição de fluência de partículas dentro do
	detector C. Desaceleração contínua
	*/

	/*	if (imprimiu==0){
			printf("\n\nFINDES2\n\n");
			imprimiu++;
		}*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int IEI, IEF;

	double EIC, EF, FACT, EA, EB, TLBIN;

	if (EI < dg_CIMDET_.EL[ID - 1])
		return;
	EIC = fmin(EI, dg_CIMDET_.EU[ID - 1]);
	EF = fmax(EI - DECSD, dg_CIMDET_.EL[ID - 1]);

	if (EF > EIC)
		return;

	if (dg_CIMDET_.LLE[ID - 1] == 1)
	{
		IEI = (int)(1.0e0 + (log(EIC) - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
		IEF = (int)(1.0e0 + (log(EF) - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}
	else
	{
		IEI = (int)(1.0e0 + (EIC - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
		IEF = (int)(1.0e0 + (EF - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}

	FACT = DSEF / DECSD;

	for (int IE = IEF; IE <= IEI; IE++)
	{
		EA = fmax(EF, dg_CIMDET_.ET[IE - 1][ID - 1]);
		EB = fmin(EIC, dg_CIMDET_.ET[IE + 1 - 1][ID - 1]);
		TLBIN = (EB - EA) * FACT;
		if (N != dg_CIMDET_.LFLT[IE - 1][ID - 1])
		{
			dg_CIMDET_.FLT[IE - 1][ID - 1] = dg_CIMDET_.FLT[IE - 1][ID - 1] + dg_CIMDET_.FLTP[IE - 1][ID - 1];
			dg_CIMDET_.FLT2[IE - 1][ID - 1] = dg_CIMDET_.FLT2[IE - 1][ID - 1] + pow(dg_CIMDET_.FLTP[IE - 1][ID - 1], 2);
			dg_CIMDET_.FLTP[IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT * TLBIN;
			dg_CIMDET_.LFLT[IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.FLTP[IE - 1][ID - 1] = dg_CIMDET_.FLTP[IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT * TLBIN;
		}

		if (N != dg_CIMDET_.LFLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1])
		{
			dg_CIMDET_.FLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1];
			dg_CIMDET_.FLP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + pow(dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1], 2);
			dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT * TLBIN;
			dg_CIMDET_.LFLP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.FLPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT * TLBIN;
		}
	}
}

__global__ void g_knock2_G(int size)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double CDT, DF, STS, SS, ES, EP, CDTS;
	double DFS, US, VS, WS, ECDT, CONS, EE, CDTE, CDTP;
	int IZA, ISA, IEFF;
	if ((index < size) && (dg_TRACK_mod_[index].STEP == 7))
	{

		//if ((dg_TRACK_mod_[index].IEXIT == 0) && (dg_wSHOWERS_[trackIndex].CROSS == false))
		//{

			dg_TRACK_mod_[index].STEP = 8; //g_showers_step5_G/E/P
			// Fotons KPAR=2

			// L2000:;

			STS = dg_CJUMP0_[trackIndex].ST * d_rand2_(1.0e0);
			SS = dg_CJUMP0_[trackIndex].P[1 - 1];
			if (SS > STS)
				goto L2100;

			SS = SS + dg_CJUMP0_[trackIndex].P[2 - 1];
			if (SS > STS)
				goto L2200;

			SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
			if (SS > STS)
				goto L2300;

			SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
			if (SS > STS)
				goto L2400;

			SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
			if (SS > STS)
				goto L2800;

			// Espalhamento Rayleigh ICOL=1

		L2100:;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			d_graa2_(dg_TRACK_mod_[index].E, CDT, IEFF, dg_TRACK_mod_[index].MAT);

			/*
			Interação delta. Introduzido para corrigir o uso de um
			limite superior do coeficiente de atenuação de Rayleigh.
			*/

			if (IEFF == 0)
			{
				dg_wSHOWERS_[trackIndex].ICOL = 7;
				return;
			}
			dg_wSHOWERS_[trackIndex].ICOL = 1;
			if (dg_TRACK_mod_[index].IPOL == 1)
			{
				double wCONS = 0.0e0;
				d_dirpol2_(CDT, DF, wCONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}
			else
			{
				DF = d_TWOPI * d_rand2_(2.0e0);
				d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}

			dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
			dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
			dg_TRACK_mod_[index].ILB[3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
			return;

			// Espalhamento Compton ICOL=2

		L2200:;
			dg_wSHOWERS_[trackIndex].ICOL = 2;
			d_gcoa2_(dg_TRACK_mod_[index].E, dg_wSHOWERS_[trackIndex].DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
			US = dg_TRACK_mod_[index].U;
			VS = dg_TRACK_mod_[index].V;
			WS = dg_TRACK_mod_[index].W;
			DF = -1.0e0;
			if ((IZA > 0) && (ISA < 17))
			{
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				d_relax2_(IZA, ISA);
			}

			// Nova direção e energia
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
			{
				if (dg_TRACK_mod_[index].IPOL == 1)
				{
					ECDT = dg_TRACK_mod_[index].E * d_RREV * (1.0e0 - CDT);
					CONS = ECDT * ECDT / (1.0e0 + ECDT);
					d_dirpol2_(CDT, DF, CONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
				}
				else
				{
					DF = d_TWOPI * d_rand2_(3.0e0);
					d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
				}
				dg_TRACK_mod_[index].E = EP;
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E;
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			// Electron Compton - particula secundaria
			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				if (DF < -0.5e0)
					DF = d_TWOPI * d_rand2_(4.0e0);
				DFS = DF + d_PI;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wKPARP = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
			dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
			dg_TRACK_mod_[index].ILB[3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
			return;

			// Absorção Fotoeletrica ICOL=3

		L2300:;
			dg_wSHOWERS_[trackIndex].ICOL = 3;
			d_gpha2_(ES, IZA, ISA);
			/*
			Interação delta. Introduzido para corrigir o uso de um
			limite superior do coeficiente de atenuação fotoelétrica.
			*/

			if (IZA == 0)
			{
				dg_wSHOWERS_[trackIndex].ICOL = 7;
				dg_wSHOWERS_[trackIndex].DE = 0.0e0;
				return;
			}

			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				d_sauter2_(ES, CDTS);
				DFS = d_TWOPI * d_rand2_(5.0e0);
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wKPARP = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			if (ISA < 17)
			{
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				d_relax2_(IZA, ISA);
			}
			dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E;
			dg_TRACK_mod_[index].E = 0.0e0;
			return;

			// Produção de pares eletron-positron ICOL=4

		L2400:;

			dg_wSHOWERS_[trackIndex].ICOL = 4;
			d_gppa2_(EE, CDTE, EP, CDTP, IZA, ISA);
			dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E;
			// Eletron
			if (EE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				DF = d_TWOPI * d_rand2_(6.0e0);
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTE, DF, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wKPARP = 1;
				d_stores2_(EE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			// Positron
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
			{
				DF = d_TWOPI * d_rand2_(7.0e0);
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTP, DF, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wKPARP = 3;
				d_stores2_(EP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
				// O pósitron carrega uma energia 'latente' de 1022 keV.
				dg_wSHOWERS_[trackIndex].DE = dg_wSHOWERS_[trackIndex].DE - d_TREV;
			}
			else
			{
				d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
			}
			dg_TRACK_mod_[index].E = 0.0e0;

			// Relaxamento atômico após a produção de tripletos.
			if (ISA < 17)
			{
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				d_relax2_(IZA, ISA);
			}
			return;

			// Mecanismo fictício auxiliar (ICOL=8).

		L2800:;

			dg_wSHOWERS_[trackIndex].ICOL = 8;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			d_gaux2_();
			return;
		//}
	}
}

__device__ void d_knock2_G(double &DE, int &ICOL)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double CDT, DF, STS, SS, ES, EP, CDTS;
	double DFS, US, VS, WS, ECDT, CONS, EE, CDTE, CDTP;
	int IZA, ISA, IEFF;

	// Fotons KPAR=2

	// L2000:;

	STS = dg_CJUMP0_[trackIndex].ST * d_rand2_(1.0e0);
	SS = dg_CJUMP0_[trackIndex].P[1 - 1];
	if (SS > STS)
		goto L2100;

	SS = SS + dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L2200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L2300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L2400;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L2800;

	// Espalhamento Rayleigh ICOL=1

L2100:;
	DE = 0.0e0;
	d_graa2_(dg_TRACK_mod_[index].E, CDT, IEFF, dg_TRACK_mod_[index].MAT);

	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação de Rayleigh.
	*/

	if (IEFF == 0)
	{
		ICOL = 7;
		return;
	}
	ICOL = 1;
	if (dg_TRACK_mod_[index].IPOL == 1)
	{
		double wCONS = 0.0e0;
		d_dirpol2_(CDT, DF, wCONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DF = d_TWOPI * d_rand2_(2.0e0);
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}

	dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
	dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
	dg_TRACK_mod_[index].ILB[3 - 1] = ICOL;
	return;

	// Espalhamento Compton ICOL=2

L2200:;
	ICOL = 2;
	d_gcoa2_(dg_TRACK_mod_[index].E, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	US = dg_TRACK_mod_[index].U;
	VS = dg_TRACK_mod_[index].V;
	WS = dg_TRACK_mod_[index].W;
	DF = -1.0e0;
	if ((IZA > 0) && (ISA < 17))
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	// Nova direção e energia
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		if (dg_TRACK_mod_[index].IPOL == 1)
		{
			ECDT = dg_TRACK_mod_[index].E * d_RREV * (1.0e0 - CDT);
			CONS = ECDT * ECDT / (1.0e0 + ECDT);
			d_dirpol2_(CDT, DF, CONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
		}
		else
		{
			DF = d_TWOPI * d_rand2_(3.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
		}
		dg_TRACK_mod_[index].E = EP;
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	// Electron Compton - particula secundaria
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		if (DF < -0.5e0)
			DF = d_TWOPI * d_rand2_(4.0e0);
		DFS = DF + d_PI;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
	dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
	dg_TRACK_mod_[index].ILB[3 - 1] = ICOL;
	return;

	// Absorção Fotoeletrica ICOL=3

L2300:;
	ICOL = 3;
	d_gpha2_(ES, IZA, ISA);
	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação fotoelétrica.
	*/

	if (IZA == 0)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		d_sauter2_(ES, CDTS);
		DFS = d_TWOPI * d_rand2_(5.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	if (ISA < 17)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}
	DE = dg_TRACK_mod_[index].E;
	dg_TRACK_mod_[index].E = 0.0e0;
	return;

	// Produção de pares eletron-positron ICOL=4

L2400:;

	ICOL = 4;
	d_gppa2_(EE, CDTE, EP, CDTP, IZA, ISA);
	DE = dg_TRACK_mod_[index].E;
	// Eletron
	if (EE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DF = d_TWOPI * d_rand2_(6.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTE, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(EE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Positron
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		DF = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTP, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 3;
		d_stores2_(EP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
		// O pósitron carrega uma energia 'latente' de 1022 keV.
		DE = DE - d_TREV;
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
	}
	dg_TRACK_mod_[index].E = 0.0e0;

	// Relaxamento atômico após a produção de tripletos.
	if (ISA < 17)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L2800:;

	ICOL = 8;
	DE = 0.0e0;
	d_gaux2_();
	return;
}

__device__ void d_knock2_E(double &DE, int &ICOL)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS;
	int IOSC, IZA, ISA;

	//	L1000:;
	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		goto L1100;

	// Eletrons
	// Dobradiça, evento suave artificial (ICOL=1).

	ICOL = 1;
	dg_CJUMP1_[trackIndex].MHINGE = 1;

	// Perda de Energia

	if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
	{
		DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
		dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
		if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
		{
			DE = dg_PENELOPE_mod_.E0STEP[trackIndex];
			dg_TRACK_mod_[index].E = 0.0e0;
			return;
		}
		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
		if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
			return;
		dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
	}
	else
	{
		DE = 0.0e0;
	}

	// Deflexão Angular
	if (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
	if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
		return;
	// 1º e 2º momentos da distribuição angular.
	EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
	EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
	// Amostragem de um histograma de duas barras com esses momentos.
	PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
	PDEN = 1.0e0 - 2.0e0 * EMU1;
	PMU0 = PNUM / PDEN;
	PA = PDEN + PMU0;
	RND = d_rand2_(2.0e0);

	if (RND < PA)
	{
		CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
	}
	else
	{
		CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
	}
	DF = d_TWOPI * d_rand2_(3.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	return;

	// Evento duro

L1100:;

	dg_CJUMP1_[trackIndex].MHINGE = 0;
	// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (dg_CJUMP1_[trackIndex].KDELTA == 1)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	// Amostragem aleatória do tipo de interação.
	STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
	SS = dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L1200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L1300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L1400;

	SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
	if (SS > STS)
		goto L1500;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L1800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL = 7;
	DE = 0.0e0;
	return;

	// Colisão Elastica Dura ICOL=2

L1200:;
	ICOL = 2;
	if (dg_TRACK_mod_[index].E >= dg_CELSEP_.EELMAX[dg_TRACK_mod_[index].MAT - 1])
	{
		TRNDC = dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		TA = exp(dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		TB = dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eela2_(TA, TB, TRNDC, RMU);
	}
	else
	{
		// Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC = dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eeld2_(TRNDC, RMU);
	}

	CDT = 1.0e0 - (RMU + RMU);
	DF = d_TWOPI * d_rand2_(5.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	DE = 0.0e0;
	return;

	// Colisão Dura inelastica (ICOL=3)

L1300:;
	ICOL = 3;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_eina2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(6.0e0);
	// Raio Delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia e direção.
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	return;

	// Emissão bremsstrahlung dura (ICOL=4).
L1400:;
	ICOL = 4;
	d_ebra2_(dg_TRACK_mod_[index].E, DE, dg_TRACK_mod_[index].MAT);

	// fóton de Bremsstrahlung.

	if (DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		d_ebraa2_(dg_TRACK_mod_[index].E, DE, CDTS, dg_TRACK_mod_[index].MAT);
		DFS = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia
	dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
	if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DE = dg_TRACK_mod_[index].E + DE;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Ionização de uma casca interna (ICOL=5).

L1500:;
	ICOL = 5;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_esia2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	// relaxamento atomico
	if (IZA > 2)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(8.0e0);
	// raio delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}

	// Nova energia e direção
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L1800:;
	ICOL = 8;
	DE = 0.0e0;
	d_eaux2_();
	return;
}

__global__ void g_knock2_E(int size)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS;
	int IOSC, IZA, ISA;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 7))
	{
		int trackIndex = dg_TRACK_mod_[index].INDEX;
		//if ((dg_TRACK_mod_[index].IEXIT == 0) && (dg_wSHOWERS_[trackIndex].CROSS == false))
		//{


			dg_TRACK_mod_[index].STEP = 8; //g_showers_step5_G/E/P
			//	L1000:;
			if (dg_CJUMP1_[trackIndex].MHINGE == 1)
				goto L1100;

			// Eletrons
			// Dobradiça, evento suave artificial (ICOL=1).

			dg_wSHOWERS_[trackIndex].ICOL = 1;
			dg_CJUMP1_[trackIndex].MHINGE = 1;

			// Perda de Energia

			if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
			{
				dg_wSHOWERS_[trackIndex].DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
				dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - dg_wSHOWERS_[trackIndex].DE;
				if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
				{
					dg_wSHOWERS_[trackIndex].DE = dg_PENELOPE_mod_.E0STEP[trackIndex];
					dg_TRACK_mod_[index].E = 0.0e0;
					return;
				}
				dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
				if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
					return;
				dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			}

			// Deflexão Angular
			if (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
			{
				dg_CJUMP0_[trackIndex].T1 = exp(dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
				dg_CJUMP0_[trackIndex].T2 = exp(dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
			}
			else
			{
				dg_CJUMP0_[trackIndex].T1 = 0.0e0;
				dg_CJUMP0_[trackIndex].T2 = 0.0e0;
			}
			if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
				return;
			// 1º e 2º momentos da distribuição angular.
			EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
			EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
			// Amostragem de um histograma de duas barras com esses momentos.
			PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
			PDEN = 1.0e0 - 2.0e0 * EMU1;
			PMU0 = PNUM / PDEN;
			PA = PDEN + PMU0;
			RND = d_rand2_(2.0e0);

			if (RND < PA)
			{
				CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
			}
			else
			{
				CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
			}
			DF = d_TWOPI * d_rand2_(3.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			return;

			// Evento duro

		L1100:;

			dg_CJUMP1_[trackIndex].MHINGE = 0;
			// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
			if (dg_CJUMP1_[trackIndex].KDELTA == 1)
			{
				dg_wSHOWERS_[trackIndex].ICOL = 7;
				dg_wSHOWERS_[trackIndex].DE = 0.0e0;
				return;
			}

			// Amostragem aleatória do tipo de interação.
			STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
			STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
			SS = dg_CJUMP0_[trackIndex].P[2 - 1];
			if (SS > STS)
				goto L1200;

			SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
			if (SS > STS)
				goto L1300;

			SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
			if (SS > STS)
				goto L1400;

			SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
			if (SS > STS)
				goto L1500;

			SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
			if (SS > STS)
				goto L1800;

			/*
			Uma interação delta (ICOL=7) pode ocorrer quando o total
			A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
			*/

			dg_wSHOWERS_[trackIndex].ICOL = 7;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			return;

			// Colisão Elastica Dura ICOL=2

		L1200:;
			dg_wSHOWERS_[trackIndex].ICOL = 2;
			if (dg_TRACK_mod_[index].E >= dg_CELSEP_.EELMAX[dg_TRACK_mod_[index].MAT - 1])
			{
				TRNDC = dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				TA = exp(dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
				TB = dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				d_eela2_(TA, TB, TRNDC, RMU);
			}
			else
			{
				// Implementacao do modelo alternativo utilzando  ELSEPA database
				TRNDC = dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				d_eeld2_(TRNDC, RMU);
			}

			CDT = 1.0e0 - (RMU + RMU);
			DF = d_TWOPI * d_rand2_(5.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			return;

			// Colisão Dura inelastica (ICOL=3)

		L1300:;
			dg_wSHOWERS_[trackIndex].ICOL = 3;
			DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
			d_eina2_(dg_TRACK_mod_[index].E, DELTA, dg_wSHOWERS_[trackIndex].DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

			//Ângulos de espalhamento (elétron primário).
			DF = d_TWOPI * d_rand2_(6.0e0);
			// Raio Delta
			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				DFS = DF + d_PI;
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			// Nova energia e direção.
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				dg_TRACK_mod_[index].E = EP;
				d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E;
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			return;

			// Emissão bremsstrahlung dura (ICOL=4).
		L1400:;
			dg_wSHOWERS_[trackIndex].ICOL = 4;
			d_ebra2_(dg_TRACK_mod_[index].E, dg_wSHOWERS_[trackIndex].DE, dg_TRACK_mod_[index].MAT);

			// fóton de Bremsstrahlung.

			if (dg_wSHOWERS_[trackIndex].DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
			{
				d_ebraa2_(dg_TRACK_mod_[index].E, dg_wSHOWERS_[trackIndex].DE, CDTS, dg_TRACK_mod_[index].MAT);
				DFS = d_TWOPI * d_rand2_(7.0e0);
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 2;
				d_stores2_(dg_wSHOWERS_[trackIndex].DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			// Nova energia
			dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - dg_wSHOWERS_[trackIndex].DE;
			if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E + dg_wSHOWERS_[trackIndex].DE;
				dg_TRACK_mod_[index].E = 0.0e0;
			}
			return;

			// Ionização de uma casca interna (ICOL=5).

		L1500:;
			dg_wSHOWERS_[trackIndex].ICOL = 5;
			DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
			d_esia2_(dg_TRACK_mod_[index].E, DELTA, dg_wSHOWERS_[trackIndex].DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
			// relaxamento atomico
			if (IZA > 2)
			{
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				d_relax2_(IZA, ISA);
			}

			//Ângulos de espalhamento (elétron primário).
			DF = d_TWOPI * d_rand2_(8.0e0);
			// raio delta
			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				DFS = DF + d_PI;
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}

			// Nova energia e direção
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				dg_TRACK_mod_[index].E = EP;
				d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E;
				dg_TRACK_mod_[index].E = 0.0e0;
			}
			return;

			// Mecanismo fictício auxiliar (ICOL=8).

		L1800:;
			dg_wSHOWERS_[trackIndex].ICOL = 8;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			d_eaux2_();
			return;
		//}
	}
}

__device__ void d_knock2_P(double &DE, int &ICOL)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS, E1, CDT1, E2, CDT2;
	int IOSC, IZA, ISA;

	// L3000:;

	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		goto L3100;

	// Dobradiça, evento suave artificial (ICOL=1).

	ICOL = 1;
	dg_CJUMP1_[trackIndex].MHINGE = 1;

	// Perda de Energia

	if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
	{
		DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
		dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
		if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
		{
			d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
			DE = dg_PENELOPE_mod_.E0STEP[trackIndex] + d_TREV;
			dg_TRACK_mod_[index].E = 0.0e0;
			return;
		}
		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
		if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
			return;

		dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
	}
	else
	{
		DE = 0.0e0;
	}

	// Deflexão Angular
	if (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
	if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
		return;
	// 1º e 2º momentos da distribuição angular.
	EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
	EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
	// Amostragem de um histograma de duas barras com esses momentos.
	PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
	PDEN = 1.0e0 - 2.0e0 * EMU1;
	PMU0 = PNUM / PDEN;
	PA = PDEN + PMU0;
	RND = d_rand2_(2.0e0);

	if (RND < PA)
	{
		CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
	}
	else
	{
		CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
	}
	DF = d_TWOPI * d_rand2_(3.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	return;

	// Evento duro

L3100:;

	dg_CJUMP1_[trackIndex].MHINGE = 0;
	// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (dg_CJUMP1_[trackIndex].KDELTA == 1)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	// Amostragem aleatória do tipo de interação.
	STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
	SS = dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L3200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L3300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L3400;

	SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
	if (SS > STS)
		goto L3500;

	SS = SS + dg_CJUMP0_[trackIndex].P[6 - 1];
	if (SS > STS)
		goto L3600;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L3800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL = 7;
	DE = 0.0e0;
	return;

	// Colisão Elastica Dura ICOL=2

L3200:;
	ICOL = 2;
	if (dg_TRACK_mod_[index].E >= dg_CELSEP_.PELMAX[dg_TRACK_mod_[index].MAT - 1])
	{
		TRNDC = dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		TA = exp(dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		TB = dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eela2_(TA, TB, TRNDC, RMU);
	}
	else
	{
		// Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC = dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_peld2_(TRNDC, RMU);
	}

	CDT = 1.0e0 - (RMU + RMU);
	DF = d_TWOPI * d_rand2_(5.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	DE = 0.0e0;
	return;

	// Colisão Dura inelastica (ICOL=3)

L3300:;
	ICOL = 3;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_pina2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

	//Ângulos de espalhamento (positron primário).
	DF = d_TWOPI * d_rand2_(6.0e0);
	// Raio Delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia e direção.
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
		DE = dg_TRACK_mod_[index].E + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	return;

	// Emissão bremsstrahlung dura (ICOL=4).
L3400:;
	ICOL = 4;
	d_ebra2_(dg_TRACK_mod_[index].E, DE, dg_TRACK_mod_[index].MAT);

	// fóton de Bremsstrahlung.

	if (DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		d_ebraa2_(dg_TRACK_mod_[index].E, DE, CDTS, dg_TRACK_mod_[index].MAT);
		DFS = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia
	dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
	if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
		DE = dg_TRACK_mod_[index].E + DE + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Ionização de uma casca interna (ICOL=5).

L3500:;
	ICOL = 5;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_psia2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	// relaxamento atomico
	if (IZA > 2)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(8.0e0);
	// raio delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}

	// Nova energia e direção
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
		DE = dg_TRACK_mod_[index].E + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

L3600:; // Aniquilação de positron em voo

	ICOL = 6;
	d_pana2_(dg_TRACK_mod_[index].E, E1, CDT1, E2, CDT2, dg_TRACK_mod_[index].MAT);
	DF = d_TWOPI * d_rand2_(9.0e0);
	if (E1 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDT1, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(E1, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	if (E2 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		DF = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDT2, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(E2, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	DE = dg_TRACK_mod_[index].E + d_TREV;
	dg_TRACK_mod_[index].E = 0.0e0;
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L3800:;
	ICOL = 8;
	DE = 0.0e0;
	d_paux2_();
	return;
}

__global__ void g_knock2_P(int size)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS, E1, CDT1, E2, CDT2;
	int IOSC, IZA, ISA;

	if (index < size)
	{

		if (dg_TRACK_mod_[index].IEXIT == 0)
		{

			// L3000:;

			if (dg_CJUMP1_[trackIndex].MHINGE == 1)
				goto L3100;

			// Dobradiça, evento suave artificial (ICOL=1).

			dg_wSHOWERS_[trackIndex].ICOL = 1;
			dg_CJUMP1_[trackIndex].MHINGE = 1;

			// Perda de Energia

			if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
			{
				dg_wSHOWERS_[trackIndex].DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
				dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - dg_wSHOWERS_[trackIndex].DE;
				if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
				{
					d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
					dg_wSHOWERS_[trackIndex].DE = dg_PENELOPE_mod_.E0STEP[trackIndex] + d_TREV;
					dg_TRACK_mod_[index].E = 0.0e0;
					return;
				}
				dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
				if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
					return;

				dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			}

			// Deflexão Angular
			if (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
			{
				dg_CJUMP0_[trackIndex].T1 = exp(dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
				dg_CJUMP0_[trackIndex].T2 = exp(dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
			}
			else
			{
				dg_CJUMP0_[trackIndex].T1 = 0.0e0;
				dg_CJUMP0_[trackIndex].T2 = 0.0e0;
			}
			if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
				return;
			// 1º e 2º momentos da distribuição angular.
			EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
			EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
			// Amostragem de um histograma de duas barras com esses momentos.
			PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
			PDEN = 1.0e0 - 2.0e0 * EMU1;
			PMU0 = PNUM / PDEN;
			PA = PDEN + PMU0;
			RND = d_rand2_(2.0e0);

			if (RND < PA)
			{
				CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
			}
			else
			{
				CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
			}
			DF = d_TWOPI * d_rand2_(3.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			return;

			// Evento duro

		L3100:;

			dg_CJUMP1_[trackIndex].MHINGE = 0;
			// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
			if (dg_CJUMP1_[trackIndex].KDELTA == 1)
			{
				dg_wSHOWERS_[trackIndex].ICOL = 7;
				dg_wSHOWERS_[trackIndex].DE = 0.0e0;
				return;
			}

			// Amostragem aleatória do tipo de interação.
			STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
			STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
			SS = dg_CJUMP0_[trackIndex].P[2 - 1];
			if (SS > STS)
				goto L3200;

			SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
			if (SS > STS)
				goto L3300;

			SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
			if (SS > STS)
				goto L3400;

			SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
			if (SS > STS)
				goto L3500;

			SS = SS + dg_CJUMP0_[trackIndex].P[6 - 1];
			if (SS > STS)
				goto L3600;

			SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
			if (SS > STS)
				goto L3800;

			/*
			Uma interação delta (ICOL=7) pode ocorrer quando o total
			A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
			*/

			dg_wSHOWERS_[trackIndex].ICOL = 7;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			return;

			// Colisão Elastica Dura ICOL=2

		L3200:;
			dg_wSHOWERS_[trackIndex].ICOL = 2;
			if (dg_TRACK_mod_[index].E >= dg_CELSEP_.PELMAX[dg_TRACK_mod_[index].MAT - 1])
			{
				TRNDC = dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				TA = exp(dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
				TB = dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				d_eela2_(TA, TB, TRNDC, RMU);
			}
			else
			{
				// Implementacao do modelo alternativo utilzando  ELSEPA database
				TRNDC = dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
				d_peld2_(TRNDC, RMU);
			}

			CDT = 1.0e0 - (RMU + RMU);
			DF = d_TWOPI * d_rand2_(5.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			return;

			// Colisão Dura inelastica (ICOL=3)

		L3300:;
			dg_wSHOWERS_[trackIndex].ICOL = 3;
			DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
			d_pina2_(dg_TRACK_mod_[index].E, DELTA, dg_wSHOWERS_[trackIndex].DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

			//Ângulos de espalhamento (positron primário).
			DF = d_TWOPI * d_rand2_(6.0e0);
			// Raio Delta
			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				DFS = DF + d_PI;
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			// Nova energia e direção.
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
			{
				dg_TRACK_mod_[index].E = EP;
				d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}
			else
			{
				d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E + d_TREV;
				dg_TRACK_mod_[index].E = 0.0e0;
			}

			return;

			// Emissão bremsstrahlung dura (ICOL=4).
		L3400:;
			dg_wSHOWERS_[trackIndex].ICOL = 4;
			d_ebra2_(dg_TRACK_mod_[index].E, dg_wSHOWERS_[trackIndex].DE, dg_TRACK_mod_[index].MAT);

			// fóton de Bremsstrahlung.

			if (dg_wSHOWERS_[trackIndex].DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
			{
				d_ebraa2_(dg_TRACK_mod_[index].E, dg_wSHOWERS_[trackIndex].DE, CDTS, dg_TRACK_mod_[index].MAT);
				DFS = d_TWOPI * d_rand2_(7.0e0);
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 2;
				d_stores2_(dg_wSHOWERS_[trackIndex].DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			// Nova energia
			dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - dg_wSHOWERS_[trackIndex].DE;
			if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
			{
				d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E + dg_wSHOWERS_[trackIndex].DE + d_TREV;
				dg_TRACK_mod_[index].E = 0.0e0;
			}
			return;

			// Ionização de uma casca interna (dg_wSHOWERS_[trackIndex].ICOL=5).

		L3500:;
			dg_wSHOWERS_[trackIndex].ICOL = 5;
			DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
			d_psia2_(dg_TRACK_mod_[index].E, DELTA, dg_wSHOWERS_[trackIndex].DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
			// relaxamento atomico
			if (IZA > 2)
			{
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				d_relax2_(IZA, ISA);
			}

			//Ângulos de espalhamento (elétron primário).
			DF = d_TWOPI * d_rand2_(8.0e0);
			// raio delta
			if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
			{
				DFS = DF + d_PI;
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDTS, DFS, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 1;
				d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}

			// Nova energia e direção
			if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
			{
				dg_TRACK_mod_[index].E = EP;
				d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
			}
			else
			{
				d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
				dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E + d_TREV;
				dg_TRACK_mod_[index].E = 0.0e0;
			}
			return;

		L3600:; // Aniquilação de positron em voo

			dg_wSHOWERS_[trackIndex].ICOL = 6;
			d_pana2_(dg_TRACK_mod_[index].E, E1, CDT1, E2, CDT2, dg_TRACK_mod_[index].MAT);
			DF = d_TWOPI * d_rand2_(9.0e0);
			if (E1 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
			{
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDT1, DF, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 2;
				d_stores2_(E1, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			if (E2 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
			{
				DF = DF + d_PI;
				US = dg_TRACK_mod_[index].U;
				VS = dg_TRACK_mod_[index].V;
				WS = dg_TRACK_mod_[index].W;
				d_direct2_(CDT2, DF, US, VS, WS);
				dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
				dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
				dg_CHIST_.ILBA[trackIndex][3 - 1] = dg_wSHOWERS_[trackIndex].ICOL;
				dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
				dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
				int wvar = 2;
				d_stores2_(E2, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
			}
			dg_wSHOWERS_[trackIndex].DE = dg_TRACK_mod_[index].E + d_TREV;
			dg_TRACK_mod_[index].E = 0.0e0;
			return;

			// Mecanismo fictício auxiliar (ICOL=8).

		L3800:;
			dg_wSHOWERS_[trackIndex].ICOL = 8;
			dg_wSHOWERS_[trackIndex].DE = 0.0e0;
			d_paux2_();
			return;
		}
	}
}

__device__ void d_knock2_(double &DE, int &ICOL)
{

	/*
	Simulação de dobradiças aleatórias e eventos de interação difícil.

	Argumentos de saída:
	DE ..... energia depositada pela partícula no material. Isto é
	geralmente igual à diferença entre as energias
	antes e depois da interação.
	ICOL ... tipo de interação sofrida pela partícula.
	*/

	/*	if (imprimiu==0){
			printf("\n\nKNOCK2\n\n");
			imprimiu++;
		}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	double EMU1, EMU2, PNUM, PDEN, PMU0, PA, RND, CDT, DF, STNOW, STS, SS, TRNDC, TA, TB, RMU, DELTA, ES, EP, CDTS;
	double DFS, US, VS, WS, ECDT, CONS, EE, CDTE, CDTP, E1, CDT1, E2, CDT2;
	int IOSC, IZA, ISA, IEFF;

	if (dg_TRACK_mod_[index].KPAR == 1)
		goto L1000;
	else if (dg_TRACK_mod_[index].KPAR == 2)
		goto L2000;
	else if (dg_TRACK_mod_[index].KPAR == 3)
		goto L3000;
	else
	{
		printf("   KNOCK: Incorrect particle type.\n");
		// exit(0);
		return;
	}

L1000:;
	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		goto L1100;

	// Eletrons
	// Dobradiça, evento suave artificial (ICOL=1).

	ICOL = 1;
	dg_CJUMP1_[trackIndex].MHINGE = 1;

	// Perda de Energia

	if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
	{
		DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
		dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
		if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
		{
			DE = dg_PENELOPE_mod_.E0STEP[trackIndex];
			dg_TRACK_mod_[index].E = 0.0e0;
			return;
		}
		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
		if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
			return;
		dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
	}
	else
	{
		DE = 0.0e0;
	}

	// Deflexão Angular
	if (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
	if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
		return;
	// 1º e 2º momentos da distribuição angular.
	EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
	EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
	// Amostragem de um histograma de duas barras com esses momentos.
	PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
	PDEN = 1.0e0 - 2.0e0 * EMU1;
	PMU0 = PNUM / PDEN;
	PA = PDEN + PMU0;
	RND = d_rand2_(2.0e0);

	if (RND < PA)
	{
		CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
	}
	else
	{
		CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
	}
	DF = d_TWOPI * d_rand2_(3.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	return;

	// Evento duro

L1100:;

	dg_CJUMP1_[trackIndex].MHINGE = 0;
	// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (dg_CJUMP1_[trackIndex].KDELTA == 1)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	// Amostragem aleatória do tipo de interação.
	STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
	SS = dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L1200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L1300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L1400;

	SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
	if (SS > STS)
		goto L1500;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L1800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL = 7;
	DE = 0.0e0;
	return;

	// Colisão Elastica Dura ICOL=2

L1200:;
	ICOL = 2;
	if (dg_TRACK_mod_[index].E >= dg_CELSEP_.EELMAX[dg_TRACK_mod_[index].MAT - 1])
	{
		TRNDC = dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.RNDCE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		TA = exp(dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.AE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		TB = dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.BE[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eela2_(TA, TB, TRNDC, RMU);
	}
	else
	{
		// Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC = dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCED[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eeld2_(TRNDC, RMU);
	}

	CDT = 1.0e0 - (RMU + RMU);
	DF = d_TWOPI * d_rand2_(5.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	DE = 0.0e0;
	return;

	// Colisão Dura inelastica (ICOL=3)

L1300:;
	ICOL = 3;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_eina2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(6.0e0);
	// Raio Delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia e direção.
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	return;

	// Emissão bremsstrahlung dura (ICOL=4).
L1400:;
	ICOL = 4;
	d_ebra2_(dg_TRACK_mod_[index].E, DE, dg_TRACK_mod_[index].MAT);

	// fóton de Bremsstrahlung.

	if (DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		d_ebraa2_(dg_TRACK_mod_[index].E, DE, CDTS, dg_TRACK_mod_[index].MAT);
		DFS = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia
	dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
	if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DE = dg_TRACK_mod_[index].E + DE;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Ionização de uma casca interna (ICOL=5).

L1500:;
	ICOL = 5;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_esia2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	// relaxamento atomico
	if (IZA > 2)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(8.0e0);
	// raio delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}

	// Nova energia e direção
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L1800:;
	ICOL = 8;
	DE = 0.0e0;
	d_eaux2_();
	return;

	// Fotons KPAR=2

L2000:;

	STS = dg_CJUMP0_[trackIndex].ST * d_rand2_(1.0e0);
	SS = dg_CJUMP0_[trackIndex].P[1 - 1];
	if (SS > STS)
		goto L2100;

	SS = SS + dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L2200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L2300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L2400;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L2800;

	// Espalhamento Rayleigh ICOL=1

L2100:;
	DE = 0.0e0;
	d_graa2_(dg_TRACK_mod_[index].E, CDT, IEFF, dg_TRACK_mod_[index].MAT);

	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação de Rayleigh.
	*/

	if (IEFF == 0)
	{
		ICOL = 7;
		return;
	}
	ICOL = 1;
	if (dg_TRACK_mod_[index].IPOL == 1)
	{
		double wCONS = 0.0e0;
		d_dirpol2_(CDT, DF, wCONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		DF = d_TWOPI * d_rand2_(2.0e0);
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}

	dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
	dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
	dg_TRACK_mod_[index].ILB[3 - 1] = ICOL;
	return;

	// Espalhamento Compton ICOL=2

L2200:;
	ICOL = 2;
	d_gcoa2_(dg_TRACK_mod_[index].E, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	US = dg_TRACK_mod_[index].U;
	VS = dg_TRACK_mod_[index].V;
	WS = dg_TRACK_mod_[index].W;
	DF = -1.0e0;
	if ((IZA > 0) && (ISA < 17))
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	// Nova direção e energia
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		if (dg_TRACK_mod_[index].IPOL == 1)
		{
			ECDT = dg_TRACK_mod_[index].E * d_RREV * (1.0e0 - CDT);
			CONS = ECDT * ECDT / (1.0e0 + ECDT);
			d_dirpol2_(CDT, DF, CONS, dg_TRACK_mod_[index].SP1, dg_TRACK_mod_[index].SP2, dg_TRACK_mod_[index].SP3, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
		}
		else
		{
			DF = d_TWOPI * d_rand2_(3.0e0);
			d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
		}
		dg_TRACK_mod_[index].E = EP;
	}
	else
	{
		DE = dg_TRACK_mod_[index].E;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	// Electron Compton - particula secundaria
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		if (DF < -0.5e0)
			DF = d_TWOPI * d_rand2_(4.0e0);
		DFS = DF + d_PI;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	dg_TRACK_mod_[index].ILB[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
	dg_TRACK_mod_[index].ILB[2 - 1] = dg_TRACK_mod_[index].KPAR;
	dg_TRACK_mod_[index].ILB[3 - 1] = ICOL;
	return;

	// Absorção Fotoeletrica ICOL=3

L2300:;
	ICOL = 3;
	d_gpha2_(ES, IZA, ISA);
	/*
	Interação delta. Introduzido para corrigir o uso de um
	limite superior do coeficiente de atenuação fotoelétrica.
	*/

	if (IZA == 0)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		d_sauter2_(ES, CDTS);
		DFS = d_TWOPI * d_rand2_(5.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	if (ISA < 17)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}
	DE = dg_TRACK_mod_[index].E;
	dg_TRACK_mod_[index].E = 0.0e0;
	return;

	// Produção de pares eletron-positron ICOL=4

L2400:;

	ICOL = 4;
	d_gppa2_(EE, CDTE, EP, CDTP, IZA, ISA);
	DE = dg_TRACK_mod_[index].E;
	// Eletron
	if (EE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DF = d_TWOPI * d_rand2_(6.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTE, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 1;
		d_stores2_(EE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Positron
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		DF = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTP, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wKPARP = 3;
		d_stores2_(EP, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wKPARP, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
		// O pósitron carrega uma energia 'latente' de 1022 keV.
		DE = DE - d_TREV;
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
	}
	dg_TRACK_mod_[index].E = 0.0e0;

	// Relaxamento atômico após a produção de tripletos.
	if (ISA < 17)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L2800:;

	ICOL = 8;
	DE = 0.0e0;
	d_gaux2_();
	return;

	// Positrons KPAR=3

L3000:;

	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		goto L3100;

	// Dobradiça, evento suave artificial (ICOL=1).

	ICOL = 1;
	dg_CJUMP1_[trackIndex].MHINGE = 1;

	// Perda de Energia

	if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
	{
		DE = dg_PENELOPE_mod_.DESOFT[trackIndex];
		dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
		if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
		{
			d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
			DE = dg_PENELOPE_mod_.E0STEP[trackIndex] + d_TREV;
			dg_TRACK_mod_[index].E = 0.0e0;
			return;
		}
		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_PENELOPE_mod_.E0STEP[trackIndex] - dg_PENELOPE_mod_.SSOFT[trackIndex] * (dg_CJUMP0_[trackIndex].DST - dg_CJUMP0_[trackIndex].DSR);
		if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
			return;

		dg_CEGRID_.XEL[trackIndex] = log(dg_PENELOPE_mod_.E0STEP[trackIndex]);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
	}
	else
	{
		DE = 0.0e0;
	}

	// Deflexão Angular
	if (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
	if (dg_CJUMP0_[trackIndex].T1 < 1.0e-20)
		return;
	// 1º e 2º momentos da distribuição angular.
	EMU1 = 0.5e0 * (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T1));
	EMU2 = EMU1 - (1.0e0 - exp(-dg_CJUMP0_[trackIndex].DST * dg_CJUMP0_[trackIndex].T2)) / 6.0e0;
	// Amostragem de um histograma de duas barras com esses momentos.
	PNUM = 2.0e0 * EMU1 - 3.0e0 * EMU2;
	PDEN = 1.0e0 - 2.0e0 * EMU1;
	PMU0 = PNUM / PDEN;
	PA = PDEN + PMU0;
	RND = d_rand2_(2.0e0);

	if (RND < PA)
	{
		CDT = 1.0e0 - 2.0e0 * PMU0 * (RND / PA);
	}
	else
	{
		CDT = 1.0e0 - 2.0e0 * (PMU0 + (1.0e0 - PMU0) * ((RND - PA) / (1.0e0 - PA)));
	}
	DF = d_TWOPI * d_rand2_(3.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	return;

	// Evento duro

L3100:;

	dg_CJUMP1_[trackIndex].MHINGE = 0;
	// Uma interação delta (ICOL=7) ocorre quando o máximo comprimento de passo permitido é excedido.
	if (dg_CJUMP1_[trackIndex].KDELTA == 1)
	{
		ICOL = 7;
		DE = 0.0e0;
		return;
	}

	// Amostragem aleatória do tipo de interação.
	STNOW = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	STS = fmax(STNOW, dg_CJUMP0_[trackIndex].ST) * d_rand2_(4.0e0);
	SS = dg_CJUMP0_[trackIndex].P[2 - 1];
	if (SS > STS)
		goto L3200;

	SS = SS + dg_CJUMP0_[trackIndex].P[3 - 1];
	if (SS > STS)
		goto L3300;

	SS = SS + dg_CJUMP0_[trackIndex].P[4 - 1];
	if (SS > STS)
		goto L3400;

	SS = SS + dg_CJUMP0_[trackIndex].P[5 - 1];
	if (SS > STS)
		goto L3500;

	SS = SS + dg_CJUMP0_[trackIndex].P[6 - 1];
	if (SS > STS)
		goto L3600;

	SS = SS + dg_CJUMP0_[trackIndex].P[8 - 1];
	if (SS > STS)
		goto L3800;

	/*
	Uma interação delta (ICOL=7) pode ocorrer quando o total
	A probabilidade de interação por unidade de comprimento do caminho, ST, é maior que STNOW.
	*/

	ICOL = 7;
	DE = 0.0e0;
	return;

	// Colisão Elastica Dura ICOL=2

L3200:;
	ICOL = 2;
	if (dg_TRACK_mod_[index].E >= dg_CELSEP_.PELMAX[dg_TRACK_mod_[index].MAT - 1])
	{
		TRNDC = dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.RNDCP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		TA = exp(dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.AP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		TB = dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.BP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_eela2_(TA, TB, TRNDC, RMU);
	}
	else
	{
		// Implementacao do modelo alternativo utilzando  ELSEPA database
		TRNDC = dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CELSEP_.RNDCPD[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
		d_peld2_(TRNDC, RMU);
	}

	CDT = 1.0e0 - (RMU + RMU);
	DF = d_TWOPI * d_rand2_(5.0e0);
	d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	DE = 0.0e0;
	return;

	// Colisão Dura inelastica (ICOL=3)

L3300:;
	ICOL = 3;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_pina2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IOSC);

	//Ângulos de espalhamento (positron primário).
	DF = d_TWOPI * d_rand2_(6.0e0);
	// Raio Delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia e direção.
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
		DE = dg_TRACK_mod_[index].E + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}

	return;

	// Emissão bremsstrahlung dura (ICOL=4).
L3400:;
	ICOL = 4;
	d_ebra2_(dg_TRACK_mod_[index].E, DE, dg_TRACK_mod_[index].MAT);

	// fóton de Bremsstrahlung.

	if (DE > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		d_ebraa2_(dg_TRACK_mod_[index].E, DE, CDTS, dg_TRACK_mod_[index].MAT);
		DFS = d_TWOPI * d_rand2_(7.0e0);
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(DE, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	// Nova energia
	dg_TRACK_mod_[index].E = dg_TRACK_mod_[index].E - DE;
	if (dg_TRACK_mod_[index].E < dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]); // Aniquilação em repouso.
		DE = dg_TRACK_mod_[index].E + DE + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

	// Ionização de uma casca interna (ICOL=5).

L3500:;
	ICOL = 5;
	DELTA = dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.DEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex];
	d_psia2_(dg_TRACK_mod_[index].E, DELTA, DE, EP, CDT, ES, CDTS, dg_TRACK_mod_[index].MAT, IZA, ISA);
	// relaxamento atomico
	if (IZA > 2)
	{
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		d_relax2_(IZA, ISA);
	}

	//Ângulos de espalhamento (elétron primário).
	DF = d_TWOPI * d_rand2_(8.0e0);
	// raio delta
	if (ES > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][1 - 1])
	{
		DFS = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDTS, DFS, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 1;
		d_stores2_(ES, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}

	// Nova energia e direção
	if (EP > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][3 - 1])
	{
		dg_TRACK_mod_[index].E = EP;
		d_direct2_(CDT, DF, dg_TRACK_mod_[index].U, dg_TRACK_mod_[index].V, dg_TRACK_mod_[index].W);
	}
	else
	{
		d_panar2_(dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1]);
		DE = dg_TRACK_mod_[index].E + d_TREV;
		dg_TRACK_mod_[index].E = 0.0e0;
	}
	return;

L3600:; // Aniquilação de positron em voo

	ICOL = 6;
	d_pana2_(dg_TRACK_mod_[index].E, E1, CDT1, E2, CDT2, dg_TRACK_mod_[index].MAT);
	DF = d_TWOPI * d_rand2_(9.0e0);
	if (E1 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDT1, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(E1, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	if (E2 > dg_PENELOPE_mod_.EABS[dg_TRACK_mod_[index].MAT - 1][2 - 1])
	{
		DF = DF + d_PI;
		US = dg_TRACK_mod_[index].U;
		VS = dg_TRACK_mod_[index].V;
		WS = dg_TRACK_mod_[index].W;
		d_direct2_(CDT2, DF, US, VS, WS);
		dg_CHIST_.ILBA[trackIndex][1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
		dg_CHIST_.ILBA[trackIndex][2 - 1] = dg_TRACK_mod_[index].KPAR;
		dg_CHIST_.ILBA[trackIndex][3 - 1] = ICOL;
		dg_CHIST_.ILBA[trackIndex][4 - 1] = 0;
		dg_CHIST_.ILBA[trackIndex][5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
		int wvar = 2;
		d_stores2_(E2, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, wvar, dg_CHIST_.ILBA[trackIndex], d_wIPOLI);
	}
	DE = dg_TRACK_mod_[index].E + d_TREV;
	dg_TRACK_mod_[index].E = 0.0e0;
	return;

	// Mecanismo fictício auxiliar (ICOL=8).

L3800:;
	ICOL = 8;
	DE = 0.0e0;
	d_paux2_();
	return;
}

__device__ void d_eela2_(double &A, double &B, double &RNDC, double &RMU)
{

	/*
	Simulação de eventos elásticos duros. Modelo Wentzel modificado.

	Argumentos de entrada:
	A, B ... parâmetros de distribuição angular.
	RNDC ... probabilidade de corte.
	Valores de saída :
	RMU .... deflexão angular, =(1-CDT)/2.
	*/

	/*	if (imprimiu==0){
			printf("\n\nEELA2\n\n");
			imprimiu++;
		}*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double A1, B1, RMUAV, RND0, RND, RNDMB, BB, RMUC, PW, RNDRC;

	A1 = A + 1.0e0;

	if (B >= 0.0e0)
	{
		// Caso I

		RMUAV = A * A1 * log(A1 / A) - A;
		B1 = 1.0e0 - B;
		RND0 = B1 * A1 * RMUAV / (A + RMUAV);
		RND = RNDC + d_rand2_(1.0e0) * (1.0e0 - RNDC);

		if (RND < RND0)
		{
			RMU = RND * A / (B1 * A1 - RND);
		}
		else if (RND > RND0 + B)
		{
			RNDMB = RND - B;
			RMU = RNDMB * A / (B1 * A1 - RNDMB);
		}
		else
		{
			RMU = RMUAV;
		}
	}
	else
	{
		// Caso II
		BB = -B;
		B1 = 1.0e0 - BB;
		RMUC = RNDC * A / (B1 * A1 - RNDC);
		PW = B1 * A * (1.0e0 - RMUC) / (A + RMUC);
		if (d_rand2_(2.0e0) * (BB + PW) < BB)
		{
			RMU = 0.5e0 * (1.0e0 + sqrt(d_rand2_(3.0e0)));
		}
		else
		{
			RNDRC = d_rand2_(3.0e0) * (1.0e0 - RMUC);
			RMU = (A * RNDRC + A1 * RMUC) / (A1 - RNDRC);
		}
	}
}

__device__ void d_direct2_(double &CDT, double &DF, double &U, double &V, double &W)
{

	/*
	Esta sub-rotina calcula os novos cossenos de direção da partícula
	Velocidade após uma colisão com determinado espalhamento polar e azimutal
	ângulos.

	Entrada: U,V,W ... cossenos da direção inicial.
	CDT ..... cosseno do ângulo de dispersão polar.
	DF ...... ângulo de espalhamento azimutal (rad).

	Saída: U,V,W ... novos cossenos de direção.
	CDT e DF permanecem inalterados

	*/

	//	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double UV, UVW, FNORM, SDT, SDTSDF, SDTCDF, SUV, UN, VN;

	// Garante a Normalizacao

	UV = U * U + V * V;
	UVW = UV + W * W;

	if (fabs(UVW - 1.0e0) > 1.0e-13)
	{
		FNORM = 1.0e0 / sqrt(UVW);
		U = FNORM * U;
		V = FNORM * V;
		W = FNORM * W;
		UV = U * U + V * V;
	}

	// Calculando a nova direção

	if ((1.0e0 - fabs(CDT)) > 1.0e-8)
	{
		SDT = sqrt(1.0e0 - CDT * CDT);
	}
	else
	{
		SDT = sqrt(2.0e0 * (1.0e0 - fabs(CDT)));
	}

	if (SDT < 1.0e-13)
	{
		if (CDT < 0.0e0)
		{
			U = -U;
			V = -V;
			W = -W;
		}
	}
	else
	{
		SDTSDF = SDT * sin(DF);
		SDTCDF = SDT * cos(DF);
		if (UV > 1.0e-26)
		{
			SUV = sqrt(UV);
			UN = U / SUV;
			VN = V / SUV;
			U = U * CDT + (UN * W * SDTCDF - VN * SDTSDF);
			V = V * CDT + (VN * W * SDTCDF + UN * SDTSDF);
			W = W * CDT - SUV * SDTCDF;
		}
		else
		{
			if (W > 0.0e0)
			{
				U = SDTCDF;
				V = SDTSDF;
				W = CDT;
			}
			else
			{
				U = -SDTCDF;
				V = -SDTSDF;
				W = -CDT;
			}
		}
	}
	/*	if (imprimiu==0){
			printf("\n\nDIRECT2\n\n");
			imprimiu++;
		}*/
}

__device__ void d_panar2_(double &ECUT)
{

	/*Simulação da aniquilação de pósitrons em repouso. quanta de aniquilação
	são armazenados na pilha secundária somente quando ECUT é menor que d_REV.*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double US, VS, WS, CDT1, DF;

	int ILBA[5];

	if (d_REV < ECUT)
		return;

	US = dg_TRACK_mod_[index].U;
	VS = dg_TRACK_mod_[index].V;
	WS = dg_TRACK_mod_[index].W;
	CDT1 = -1.0e0 + 2.0e0 * d_rand2_(1.0e0);
	DF = d_TWOPI * d_rand2_(2.0e0);

	d_direct2_(CDT1, DF, US, VS, WS);

	ILBA[1 - 1] = dg_TRACK_mod_[index].ILB[1 - 1] + 1;
	ILBA[2 - 1] = 3;
	ILBA[3 - 1] = 6;
	ILBA[4 - 1] = 0;
	ILBA[5 - 1] = dg_TRACK_mod_[index].ILB[5 - 1];
	int WKPARP = 2;
	d_stores2_(d_REV, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, WKPARP, ILBA, d_wIPOLI);
	US = -US;
	VS = -VS;
	WS = -WS;
	d_stores2_(d_REV, dg_TRACK_mod_[index].X, dg_TRACK_mod_[index].Y, dg_TRACK_mod_[index].Z, US, VS, WS, dg_TRACK_mod_[index].WGHT, WKPARP, ILBA, d_wIPOLI);
}

__device__ void d_pimfp2_(int &IEND)
{
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para
	ações de pósitrons com a energia atual no material M.
	*/

	/*if (imprimiu==0){
		 printf("\n\nPIMFP2\n\n");
		 imprimiu++;
	 }*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	dg_CJUMP0_[trackIndex].P[2 - 1] = exp(dg_CPIMFP_.SPHEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPHEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPHEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[3 - 1] = exp(dg_CPIMFP_.SPHIN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPHIN[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPHIN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[4 - 1] = exp(dg_CPIMFP_.SPHBR[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPHBR[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPHBR[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[5 - 1] = exp(dg_CPIMFP_.SPISI[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPISI[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPISI[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[6 - 1] = exp(dg_CPIMFP_.SPAN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPAN[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPAN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[8 - 1] = 0.0e0;

	if (IEND == 1)
	{
		return;
	}

	if (dg_CPIMFP_.W1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].W1 = exp(dg_CPIMFP_.W1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.W1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.W1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].W2 = exp(dg_CPIMFP_.W2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.W2P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.W2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].W1 = 0.0e0;
		dg_CJUMP0_[trackIndex].W2 = 0.0e0;
	}

	if (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T1P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.T2P[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
}

__device__ void d_eimfp2_(int &IEND)
{
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para
	ações de eletrons com a energia atual no material M.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	if (dg_TRACK_mod_[index].MAT <= 0)
		printf("MATERIAL %d\n", dg_TRACK_mod_[index].MAT);

	if (dg_CEGRID_.KE[trackIndex] <= 0)
	{
		printf("CEGRID KE %d\n", dg_CEGRID_.KE[trackIndex]);
		printf("mATERIAL %d\n", dg_TRACK_mod_[index].MAT);
		printf("Energia %lf\n", dg_TRACK_mod_[index].E);
		printf("Index %d\n", index);
		printf("trackIndex %d\n", trackIndex);
		printf("CROSS %d\n", dg_wSHOWERS_[trackIndex].CROSS);
	}

	dg_CJUMP0_[trackIndex].P[2 - 1] = exp(dg_CEIMFP_.SEHEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SEHEL[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SEHEL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[3 - 1] = exp(dg_CEIMFP_.SEHIN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SEHIN[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SEHIN[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[4 - 1] = exp(dg_CEIMFP_.SEHBR[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SEHBR[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SEHBR[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[5 - 1] = exp(dg_CEIMFP_.SEISI[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SEISI[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SEISI[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[8 - 1] = 0.0e0;

	if (IEND == 1)
	{
		return;
	}

	if (dg_CEIMFP_.W1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].W1 = exp(dg_CEIMFP_.W1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.W1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.W1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].W2 = exp(dg_CEIMFP_.W2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.W2E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.W2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].W1 = 0.0e0;
		dg_CJUMP0_[trackIndex].W2 = 0.0e0;
	}

	if (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] > -78.3e0)
	{
		dg_CJUMP0_[trackIndex].T1 = exp(dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T1E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
		dg_CJUMP0_[trackIndex].T2 = exp(dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.T2E[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	else
	{
		dg_CJUMP0_[trackIndex].T1 = 0.0e0;
		dg_CJUMP0_[trackIndex].T2 = 0.0e0;
	}
}

__device__ void d_gimfp2_()
{
	/*
	Esta sub-rotina calcula os caminhos livres médios inversos para interações
	de fótons com a energia atual no material M.
	*/

	/*if (imprimiu==0){
		printf("\n\nGIMFP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	dg_CJUMP0_[trackIndex].P[1 - 1] = dg_CGIMFP_.SGRA[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1];
	dg_CJUMP0_[trackIndex].P[2 - 1] = exp(dg_CGIMFP_.SGCO[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CGIMFP_.SGCO[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CGIMFP_.SGCO[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	dg_CJUMP0_[trackIndex].P[3 - 1] = dg_CGIMFP_.SGPH[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1];

	if (dg_TRACK_mod_[index].E < 1.023e6)
	{
		dg_CJUMP0_[trackIndex].P[4 - 1] = 0.0e0;
	}
	else
	{
		dg_CJUMP0_[trackIndex].P[4 - 1] = exp(dg_CGIMFP_.SGPP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CGIMFP_.SGPP[dg_CEGRID_.KE[trackIndex] + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CGIMFP_.SGPP[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1]) * dg_CEGRID_.XEK[trackIndex]);
	}
	dg_CJUMP0_[trackIndex].P[8 - 1] = 0.0e0;
}

__device__ void d_start2_()
{

	/*
		Esta sub-rotina força o próximo evento a ser um evento soft artificial.
		Deve ser chamado quando uma nova trilha de partículas (primária ou secundária) é
		é iniciado e quando cruza uma interface.
	*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	if ((dg_TRACK_mod_[index].E < dg_CEGRID_.EMIN) || (dg_TRACK_mod_[index].E > 0.99999999e0 * dg_CEGRID_.EU))
	{
		printf("   *** Energy out of range. KPAR = %2d, E = %.5E eV\n", dg_TRACK_mod_[index].KPAR, dg_TRACK_mod_[index].E);
		for (int J = 1; J <= 5; J++)
		{
			printf("       ILB%d = %d\n", J, dg_TRACK_mod_[index].ILB[J - 1]);
		}
		printf("       EMIN = %.5E eV, EMAX = %.5E eV\n", dg_CEGRID_.EL, dg_CEGRID_.EU);
		printf("       Check the values of EABS(KPAR,M) and EMAX.\n");
		// exit(0);
		return;
	}
	dg_CJUMP1_[trackIndex].MHINGE = 0;
	dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E + 1.0e30;
	dg_CJUMP1_[trackIndex].ELAST2 = dg_CJUMP1_[trackIndex].ELAST1;

	/*	if (imprimiu==0){
		printf("\n\nSTART2\n\n");
		imprimiu++;
	}*/
}

__device__ void d_jump2_G(double &DSMAX, double &DS)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;

	// Fotons

	if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
	{
		dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
		d_gimfp2_();
		dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[1 - 1] + dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	}

	DS = -log(d_rand2_(1.0e0)) / dg_CJUMP0_[trackIndex].ST;
}

__global__ void g_jump2_G(int size)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 4))
	{
		int trackIndex = dg_TRACK_mod_[index].INDEX;
		dg_TRACK_mod_[index].STEP = 5; //g_step2
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{

			// Fotons

			if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
			{
				dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
				d_gimfp2_();
				dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
				dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[1 - 1] + dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
			}

			dg_wSHOWERS_[trackIndex].DS = -log(d_rand2_(1.0e0)) / dg_CJUMP0_[trackIndex].ST;
		}
	}
	
}

__device__ void d_jump2_E(double &DSMAX, double &DS)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
	{
		if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
		{
			dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
			dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
			dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			int wvar = 1;
			d_eimfp2_(wvar);
			dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		}
		DS = dg_CJUMP0_[trackIndex].DSR;
		return;
	}
	dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
	if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
	{
		dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
		int wvar = 2;
		d_eimfp2_(wvar);
		dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
		dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
	}

	// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
	dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

	/*
	Interações de parada suave.
	KSOFTI=1, parada suave está ativa,
	KSOFTI=0, a parada suave não está ativa.
	*/

	if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
	{
		dg_CJUMP1_[trackIndex].KSOFTI = 1;
		/*
		O comprimento máximo do passo, DSMAXP, é determinado em termos do
		valor DSMAX de entrada (que é especificado pelo usuário) e a média
		caminho livre para interações difíceis (1/ST).
		*/
		DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
		if (DSMAXP > DSMC)
		{
			DSMAXP = DSMC;
		}
		else if (DSMAXP < 1.0e-8)
		{
			DSMAXP = DSMC;
		}

		// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
		DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

		// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

		EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
		VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
		FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
		FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
		EDEM = EDE0 * FSEDE;
		VDEM = VDE0 * FSVDE;
		W21 = VDEM / EDEM;

		if (EDEM > 9.0e0 * W21)
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
		}
		else if (EDEM > 3.0e0 * W21)
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
		}
		else
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
		}

		XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		KE1 = (int)XE1;
		XEK1 = XE1 - KE1;
		STLWR = exp(dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SETOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
		dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
	}
	else
	{
		dg_CJUMP1_[trackIndex].KSOFTI = 0;
		dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
		dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
	}

	/*
	Dispersão elástica suave.
	KSOFTE=1, dispersão suave está ativa,
	KSOFTE=0, dispersão suave não está ativa.
	*/

	if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
	{
		dg_CJUMP1_[trackIndex].KSOFTE = 1;
	}
	else
	{
		dg_CJUMP1_[trackIndex].KSOFTE = 0;
	}

	/*
	Interações delta.
	KDELTA=0, segue-se uma interação difícil,
	KDELTA=1, segue-se uma interação delta.
	*/

	dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
	if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
	{
		dg_CJUMP1_[trackIndex].KDELTA = 0;
	}
	else
	{
		dg_CJUMP0_[trackIndex].DST = DSMAXP;
		dg_CJUMP1_[trackIndex].KDELTA = 1;
	}

	if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
	{
		dg_CJUMP1_[trackIndex].MHINGE = 1;
		DS = dg_CJUMP0_[trackIndex].DST;
	}
	else
	{
		DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
		dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - DS;
		if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
		{
			if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
			{
				dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
				dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
			}
			else
			{
				EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
				VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
				FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				EDE = EDE0 * FSEDE;
				VDE = VDE0 * FSVDE;

				// Geração de valores aleatórios DE com média EDE e variância VDE.
				SIGMA = sqrt(VDE);
				if (SIGMA < 0.333333333e0 * EDE)
				{
					// Distribuição gaussiana truncada.
					dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
				}
				else
				{
					RU = d_rand2_(4.0e0);
					EDE2 = EDE * EDE;
					VDE3 = 3.0e0 * VDE;
					if (EDE2 < VDE3)
					{
						PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
						if (RU < PNULL)
						{
							dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
							dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
							if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
							{
								dg_CJUMP1_[trackIndex].MHINGE = 1;
								DS = dg_CJUMP0_[trackIndex].DST;
							}
							else
							{
								dg_CJUMP1_[trackIndex].KSOFTI = 0;
							}
							return;
						}
						else
						{
							// Distribuição Uniforme
							dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
						}
					}
					else
					{
						dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
					}
				}
				dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
			}
		}
	}
}

__global__ void g_jump2_E(int size)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 4))
	{

		int trackIndex = dg_TRACK_mod_[index].INDEX;
		//if (dg_TRACK_mod_[index].IEXIT == 0)
		//{
		dg_TRACK_mod_[index].STEP = 5; //g_step2

			if (dg_CJUMP1_[trackIndex].MHINGE == 1)
			{
				if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
				{
					dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
					if (dg_TRACK_mod_[index].E == 0.0e0)
						printf("MHINGE valor energia %f\n", dg_TRACK_mod_[index].E);
					dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
					dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
					dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
					int wvar = 1;
					d_eimfp2_(wvar);
					dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
				}
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DSR;
				return;
			}
			dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
			if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
			{
				dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
				if (dg_TRACK_mod_[index].E == 0.0e0)
					printf("MHINGE valor energia %f\n", dg_TRACK_mod_[index].E);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
				int wvar = 2;
				d_eimfp2_(wvar);
				dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
				dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
			}

			// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
			dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
			DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

			/*
			Interações de parada suave.
			KSOFTI=1, parada suave está ativa,
			KSOFTI=0, a parada suave não está ativa.
			*/

			if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
			{
				dg_CJUMP1_[trackIndex].KSOFTI = 1;
				/*
				O comprimento máximo do passo, DSMAXP, é determinado em termos do
				valor DSMAX de entrada (que é especificado pelo usuário) e a média
				caminho livre para interações difíceis (1/ST).
				*/
				DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
				if (DSMAXP > DSMC)
				{
					DSMAXP = DSMC;
				}
				else if (DSMAXP < 1.0e-8)
				{
					DSMAXP = DSMC;
				}

				// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
				DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

				// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

				EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
				VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
				FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				EDEM = EDE0 * FSEDE;
				VDEM = VDE0 * FSVDE;
				W21 = VDEM / EDEM;

				if (EDEM > 9.0e0 * W21)
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
				}
				else if (EDEM > 3.0e0 * W21)
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
				}
				else
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
				}

				XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				KE1 = (int)XE1;
				XEK1 = XE1 - KE1;
				STLWR = exp(dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SETOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
				dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
			}
			else
			{
				dg_CJUMP1_[trackIndex].KSOFTI = 0;
				dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
				dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
			}

			/*
			Dispersão elástica suave.
			KSOFTE=1, dispersão suave está ativa,
			KSOFTE=0, dispersão suave não está ativa.
			*/

			if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
			{
				dg_CJUMP1_[trackIndex].KSOFTE = 1;
			}
			else
			{
				dg_CJUMP1_[trackIndex].KSOFTE = 0;
			}

			/*
			Interações delta.
			KDELTA=0, segue-se uma interação difícil,
			KDELTA=1, segue-se uma interação delta.
			*/

			dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
			if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
			{
				dg_CJUMP1_[trackIndex].KDELTA = 0;
			}
			else
			{
				dg_CJUMP0_[trackIndex].DST = DSMAXP;
				dg_CJUMP1_[trackIndex].KDELTA = 1;
			}

			if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
			{
				dg_CJUMP1_[trackIndex].MHINGE = 1;
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST;
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
				dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - dg_wSHOWERS_[trackIndex].DS;
				if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
				{
					if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
					{
						dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
						dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
					}
					else
					{
						EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
						VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
						FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
						FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
						EDE = EDE0 * FSEDE;
						VDE = VDE0 * FSVDE;

						// Geração de valores aleatórios DE com média EDE e variância VDE.
						SIGMA = sqrt(VDE);
						if (SIGMA < 0.333333333e0 * EDE)
						{
							// Distribuição gaussiana truncada.
							dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
						}
						else
						{
							RU = d_rand2_(4.0e0);
							EDE2 = EDE * EDE;
							VDE3 = 3.0e0 * VDE;
							if (EDE2 < VDE3)
							{
								PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
								if (RU < PNULL)
								{
									dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
									dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
									if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
									{
										dg_CJUMP1_[trackIndex].MHINGE = 1;
										dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST;
									}
									else
									{
										dg_CJUMP1_[trackIndex].KSOFTI = 0;
									}
									return;
								}
								else
								{
									// Distribuição Uniforme
									dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
								}
							}
							else
							{
								dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
							}
						}
						dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
					}
				}
			}
		//}
	}
}

__device__ void d_jump2_P(double &DSMAX, double &DS)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if (dg_CJUMP1_[trackIndex].MHINGE == 1)
	{
		if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
		{
			dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
			dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			dg_CEGRID_.KE[trackIndex] = (int)(dg_CEGRID_.XE[trackIndex]);
			dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			int wvar = 1;
			d_pimfp2_(wvar);
			dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		}
		DS = dg_CJUMP0_[trackIndex].DSR;
		return;
	}

	dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
	if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
	{
		dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
		dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
		dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
		int wvar = 2;
		d_pimfp2_(wvar);
		dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
		dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
	}

	// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
	dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
	DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

	/*
	Interações de parada suave.
	KSOFTI=1, parada suave está ativa,
	KSOFTI=0, a parada suave não está ativa.
	*/

	if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
	{
		dg_CJUMP1_[trackIndex].KSOFTI = 1;
		/*
		O comprimento máximo do passo, DSMAXP, é determinado em termos do
		valor DSMAX de entrada (que é especificado pelo usuário) e a média
		caminho livre para interações difíceis (1/ST).
		*/
		DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
		if (DSMAXP > DSMC)
		{
			DSMAXP = DSMC;
		}
		else if (DSMAXP < 1.0e-8)
		{
			DSMAXP = DSMC;
		}

		// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
		DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

		// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

		EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
		VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
		FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
		FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
		EDEM = EDE0 * FSEDE;
		VDEM = VDE0 * FSVDE;
		W21 = VDEM / EDEM;

		if (EDEM > 9.0e0 * W21)
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
		}
		else if (EDEM > 3.0e0 * W21)
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
		}
		else
		{
			ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
		}

		XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
		KE1 = (int)XE1;
		XEK1 = XE1 - KE1;
		STLWR = exp(dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPTOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
		dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
	}
	else
	{
		dg_CJUMP1_[trackIndex].KSOFTI = 0;
		dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
		dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
	}

	/*
	Dispersão elástica suave.
	KSOFTE=1, dispersão suave está ativa,
	KSOFTE=0, dispersão suave não está ativa.
	*/

	if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
	{
		dg_CJUMP1_[trackIndex].KSOFTE = 1;
	}
	else
	{
		dg_CJUMP1_[trackIndex].KSOFTE = 0;
	}

	/*
	Interações delta.
	KDELTA=0, segue-se uma interação difícil,
	KDELTA=1, segue-se uma interação delta.
	*/

	dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
	if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
	{
		dg_CJUMP1_[trackIndex].KDELTA = 0;
	}
	else
	{
		dg_CJUMP0_[trackIndex].DST = DSMAXP;
		dg_CJUMP1_[trackIndex].KDELTA = 1;
	}

	if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
	{
		dg_CJUMP1_[trackIndex].MHINGE = 1;
		DS = dg_CJUMP0_[trackIndex].DST;
	}
	else
	{
		DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
		dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - DS;
		if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
		{
			if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
			{
				dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
				dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
			}
			else
			{
				EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
				VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
				FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				EDE = EDE0 * FSEDE;
				VDE = VDE0 * FSVDE;

				// Geração de valores aleatórios DE com média EDE e variância VDE.
				SIGMA = sqrt(VDE);
				if (SIGMA < 0.333333333e0 * EDE)
				{
					// Distribuição gaussiana truncada.
					dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
				}
				else
				{
					RU = d_rand2_(4.0e0);
					EDE2 = EDE * EDE;
					VDE3 = 3.0e0 * VDE;
					if (EDE2 < VDE3)
					{
						PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
						if (RU < PNULL)
						{
							dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
							dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
							if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
							{
								dg_CJUMP1_[trackIndex].MHINGE = 1;
								DS = dg_CJUMP0_[trackIndex].DST;
							}
							else
							{
								dg_CJUMP1_[trackIndex].KSOFTI = 0;
							}
							return;
						}
						else
						{
							// Distribuição Uniforme
							dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
						}
					}
					else
					{
						dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
					}
				}
				dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
			}
		}
	}
}

__global__ void g_jump2_P(int size)
{

	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if (index < size)
	{
		if (dg_TRACK_mod_[index].IEXIT == 0)
		{

			if (dg_CJUMP1_[trackIndex].MHINGE == 1)
			{
				if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
				{
					dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
					dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
					dg_CEGRID_.KE[trackIndex] = (int)(dg_CEGRID_.XE[trackIndex]);
					dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
					int wvar = 1;
					d_pimfp2_(wvar);
					dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
				}
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DSR;
				return;
			}

			dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
			if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
			{
				dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
				int wvar = 2;
				d_pimfp2_(wvar);
				dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
				dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
			}

			// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
			dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
			DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

			/*
			Interações de parada suave.
			KSOFTI=1, parada suave está ativa,
			KSOFTI=0, a parada suave não está ativa.
			*/

			if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
			{
				dg_CJUMP1_[trackIndex].KSOFTI = 1;
				/*
				O comprimento máximo do passo, DSMAXP, é determinado em termos do
				valor DSMAX de entrada (que é especificado pelo usuário) e a média
				caminho livre para interações difíceis (1/ST).
				*/
				DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
				if (DSMAXP > DSMC)
				{
					DSMAXP = DSMC;
				}
				else if (DSMAXP < 1.0e-8)
				{
					DSMAXP = DSMC;
				}

				// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
				DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

				// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

				EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
				VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
				FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
				EDEM = EDE0 * FSEDE;
				VDEM = VDE0 * FSVDE;
				W21 = VDEM / EDEM;

				if (EDEM > 9.0e0 * W21)
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
				}
				else if (EDEM > 3.0e0 * W21)
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
				}
				else
				{
					ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
				}

				XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				KE1 = (int)XE1;
				XEK1 = XE1 - KE1;
				STLWR = exp(dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPTOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
				dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
			}
			else
			{
				dg_CJUMP1_[trackIndex].KSOFTI = 0;
				dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
				dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
			}

			/*
			Dispersão elástica suave.
			KSOFTE=1, dispersão suave está ativa,
			KSOFTE=0, dispersão suave não está ativa.
			*/

			if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
			{
				dg_CJUMP1_[trackIndex].KSOFTE = 1;
			}
			else
			{
				dg_CJUMP1_[trackIndex].KSOFTE = 0;
			}

			/*
			Interações delta.
			KDELTA=0, segue-se uma interação difícil,
			KDELTA=1, segue-se uma interação delta.
			*/

			dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
			if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
			{
				dg_CJUMP1_[trackIndex].KDELTA = 0;
			}
			else
			{
				dg_CJUMP0_[trackIndex].DST = DSMAXP;
				dg_CJUMP1_[trackIndex].KDELTA = 1;
			}

			if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
			{
				dg_CJUMP1_[trackIndex].MHINGE = 1;
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST;
			}
			else
			{
				dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
				dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - dg_wSHOWERS_[trackIndex].DS;
				if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
				{
					if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
					{
						dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
						dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
					}
					else
					{
						EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
						VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
						FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
						FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
						EDE = EDE0 * FSEDE;
						VDE = VDE0 * FSVDE;

						// Geração de valores aleatórios DE com média EDE e variância VDE.
						SIGMA = sqrt(VDE);
						if (SIGMA < 0.333333333e0 * EDE)
						{
							// Distribuição gaussiana truncada.
							dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
						}
						else
						{
							RU = d_rand2_(4.0e0);
							EDE2 = EDE * EDE;
							VDE3 = 3.0e0 * VDE;
							if (EDE2 < VDE3)
							{
								PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
								if (RU < PNULL)
								{
									dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
									dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
									if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
									{
										dg_CJUMP1_[trackIndex].MHINGE = 1;
										dg_wSHOWERS_[trackIndex].DS = dg_CJUMP0_[trackIndex].DST;
									}
									else
									{
										dg_CJUMP1_[trackIndex].KSOFTI = 0;
									}
									return;
								}
								else
								{
									// Distribuição Uniforme
									dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
								}
							}
							else
							{
								dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
							}
						}
						dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
					}
				}
			}
		}
	}
}

__device__ void d_jump2_(double &DSMAX, double &DS)
{
	/*
	Cálculo do caminho livre do ponto de partida até a posição
	do próximo evento e das probabilidades de ocorrência de diferentes
	eventos .

	Argumentos :
	DSMAX .... comprimento máximo permitido do passo (entrada),
	DS ....... comprimento do segmento (saída).

	Saída , através do módulo PENELOPE_mod:
	E0STEP ... energia no início do segmento,
	DESOFT ... perda de energia devido a interações suaves ao longo da etapa,
	SSOFT .... poder de parada devido a interações suaves,
	= DESOFT/passo_comprimento.
	*/

	/*	if (imprimiu==0){
		printf("\n\nJUMP2\n\n");
		imprimiu++;
	}*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int trackIndex = dg_TRACK_mod_[index].INDEX;
	double DSMAXP, DSMC, EDE0, VDE0, FSEDE, FSVDE, EDEM, VDEM, W21, ELOWER, XE1, XEK1, STLWR, EDE, VDE, SIGMA;
	double RU, EDE2, VDE3, PNULL;
	int KE1;

	if (dg_TRACK_mod_[index].KPAR == 1)
	{ // eletrons
		if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		{
			if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
			{
				dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
				int wvar = 1;
				d_eimfp2_(wvar);
				dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
			}
			DS = dg_CJUMP0_[trackIndex].DSR;
			return;
		}
		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
		if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
		{
			dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
			dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
			dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			int wvar = 2;
			d_eimfp2_(wvar);
			dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
			dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		}

		// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
		dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
		DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

		/*
		Interações de parada suave.
		KSOFTI=1, parada suave está ativa,
		KSOFTI=0, a parada suave não está ativa.
		*/

		if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
		{
			dg_CJUMP1_[trackIndex].KSOFTI = 1;
			/*
			O comprimento máximo do passo, DSMAXP, é determinado em termos do
			valor DSMAX de entrada (que é especificado pelo usuário) e a média
			caminho livre para interações difíceis (1/ST).
			*/
			DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
			if (DSMAXP > DSMC)
			{
				DSMAXP = DSMC;
			}
			else if (DSMAXP < 1.0e-8)
			{
				DSMAXP = DSMC;
			}

			// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
			DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

			// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

			EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
			VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
			FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
			FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
			EDEM = EDE0 * FSEDE;
			VDEM = VDE0 * FSVDE;
			W21 = VDEM / EDEM;

			if (EDEM > 9.0e0 * W21)
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
			}
			else if (EDEM > 3.0e0 * W21)
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
			}
			else
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
			}

			XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			KE1 = (int)XE1;
			XEK1 = XE1 - KE1;
			STLWR = exp(dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CEIMFP_.SETOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CEIMFP_.SETOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
			dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
		}
		else
		{
			dg_CJUMP1_[trackIndex].KSOFTI = 0;
			dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
			dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
		}

		/*
		Dispersão elástica suave.
		KSOFTE=1, dispersão suave está ativa,
		KSOFTE=0, dispersão suave não está ativa.
		*/

		if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
		{
			dg_CJUMP1_[trackIndex].KSOFTE = 1;
		}
		else
		{
			dg_CJUMP1_[trackIndex].KSOFTE = 0;
		}

		/*
		Interações delta.
		KDELTA=0, segue-se uma interação difícil,
		KDELTA=1, segue-se uma interação delta.
		*/

		dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
		if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
		{
			dg_CJUMP1_[trackIndex].KDELTA = 0;
		}
		else
		{
			dg_CJUMP0_[trackIndex].DST = DSMAXP;
			dg_CJUMP1_[trackIndex].KDELTA = 1;
		}

		if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
		{
			dg_CJUMP1_[trackIndex].MHINGE = 1;
			DS = dg_CJUMP0_[trackIndex].DST;
		}
		else
		{
			DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
			dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - DS;
			if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
			{
				if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
				{
					dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
					dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
				}
				else
				{
					EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
					VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
					FSEDE = fmax(1.0e0 - dg_CEIMFP_.DW1EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
					FSVDE = fmax(1.0e0 - dg_CEIMFP_.DW2EL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
					EDE = EDE0 * FSEDE;
					VDE = VDE0 * FSVDE;

					// Geração de valores aleatórios DE com média EDE e variância VDE.
					SIGMA = sqrt(VDE);
					if (SIGMA < 0.333333333e0 * EDE)
					{
						// Distribuição gaussiana truncada.
						dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
					}
					else
					{
						RU = d_rand2_(4.0e0);
						EDE2 = EDE * EDE;
						VDE3 = 3.0e0 * VDE;
						if (EDE2 < VDE3)
						{
							PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
							if (RU < PNULL)
							{
								dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
								dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
								if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
								{
									dg_CJUMP1_[trackIndex].MHINGE = 1;
									DS = dg_CJUMP0_[trackIndex].DST;
								}
								else
								{
									dg_CJUMP1_[trackIndex].KSOFTI = 0;
								}
								return;
							}
							else
							{
								// Distribuição Uniforme
								dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
							}
						}
						else
						{
							dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
						}
					}
					dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
				}
			}
		}
		return;
	}
	else if (dg_TRACK_mod_[index].KPAR == 3)
	{ // Positrons

		if (dg_CJUMP1_[trackIndex].MHINGE == 1)
		{
			if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
			{
				dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
				dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
				dg_CEGRID_.KE[trackIndex] = (int)(dg_CEGRID_.XE[trackIndex]);
				dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
				int wvar = 1;
				d_pimfp2_(wvar);
				dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
			}
			DS = dg_CJUMP0_[trackIndex].DSR;
			return;
		}

		dg_PENELOPE_mod_.E0STEP[trackIndex] = dg_TRACK_mod_[index].E;
		if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST2)
		{
			dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
			dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
			dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			int wvar = 2;
			d_pimfp2_(wvar);
			dg_CJUMP1_[trackIndex].ELAST2 = dg_TRACK_mod_[index].E;
			dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
		}

		// Caminho livre médio rígido inverso (probabilidade de interação por unidade comprimento do caminho).
		dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[5 - 1] + dg_CJUMP0_[trackIndex].P[6 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
		DSMAXP = dg_CSPGEO_.DSMAX[dg_TRACK_mod_[index].IBODY - 1];

		/*
		Interações de parada suave.
		KSOFTI=1, parada suave está ativa,
		KSOFTI=0, a parada suave não está ativa.
		*/

		if (dg_CJUMP0_[trackIndex].W1 > 1.0e-20)
		{
			dg_CJUMP1_[trackIndex].KSOFTI = 1;
			/*
			O comprimento máximo do passo, DSMAXP, é determinado em termos do
			valor DSMAX de entrada (que é especificado pelo usuário) e a média
			caminho livre para interações difíceis (1/ST).
			*/
			DSMC = 4.0e0 / dg_CJUMP0_[trackIndex].ST;
			if (DSMAXP > DSMC)
			{
				DSMAXP = DSMC;
			}
			else if (DSMAXP < 1.0e-8)
			{
				DSMAXP = DSMC;
			}

			// O valor de DSMAXP é randomizado para eliminar artefatos de dose no final da primeira etapa.
			DSMAXP = (0.5e0 + d_rand2_(1.0e0) * 0.5e0) * DSMAXP;

			// Limite superior para a probabilidade de interação ao longo da etapa (incluindo straggling de energia suave).

			EDE0 = dg_CJUMP0_[trackIndex].W1 * DSMAXP;
			VDE0 = dg_CJUMP0_[trackIndex].W2 * DSMAXP;
			FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
			FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
			EDEM = EDE0 * FSEDE;
			VDEM = VDE0 * FSVDE;
			W21 = VDEM / EDEM;

			if (EDEM > 9.0e0 * W21)
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + 3.0e0 * sqrt(VDEM)), dg_CEGRID_.EMIN);
			}
			else if (EDEM > 3.0e0 * W21)
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - (EDEM + sqrt(3.0e0 * VDEM)), dg_CEGRID_.EMIN);
			}
			else
			{
				ELOWER = fmax(dg_TRACK_mod_[index].E - 1.5e0 * (EDEM + W21), dg_CEGRID_.EMIN);
			}

			XE1 = 1.0e0 + (log(ELOWER) - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			KE1 = (int)XE1;
			XEK1 = XE1 - KE1;
			STLWR = exp(dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1] + (dg_CPIMFP_.SPTOT[KE1 + 1 - 1][dg_TRACK_mod_[index].MAT - 1] - dg_CPIMFP_.SPTOT[KE1 - 1][dg_TRACK_mod_[index].MAT - 1]) * XEK1);
			dg_CJUMP0_[trackIndex].ST = fmax(dg_CJUMP0_[trackIndex].ST, STLWR);
		}
		else
		{
			dg_CJUMP1_[trackIndex].KSOFTI = 0;
			dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
			dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
		}

		/*
		Dispersão elástica suave.
		KSOFTE=1, dispersão suave está ativa,
		KSOFTE=0, dispersão suave não está ativa.
		*/

		if (dg_CJUMP0_[trackIndex].T1 > 1.0e-20)
		{
			dg_CJUMP1_[trackIndex].KSOFTE = 1;
		}
		else
		{
			dg_CJUMP1_[trackIndex].KSOFTE = 0;
		}

		/*
		Interações delta.
		KDELTA=0, segue-se uma interação difícil,
		KDELTA=1, segue-se uma interação delta.
		*/

		dg_CJUMP0_[trackIndex].DST = -log(d_rand2_(2.0e0)) / dg_CJUMP0_[trackIndex].ST;
		if (dg_CJUMP0_[trackIndex].DST < DSMAXP)
		{
			dg_CJUMP1_[trackIndex].KDELTA = 0;
		}
		else
		{
			dg_CJUMP0_[trackIndex].DST = DSMAXP;
			dg_CJUMP1_[trackIndex].KDELTA = 1;
		}

		if (dg_CJUMP1_[trackIndex].KSOFTE + dg_CJUMP1_[trackIndex].KSOFTI == 0)
		{
			dg_CJUMP1_[trackIndex].MHINGE = 1;
			DS = dg_CJUMP0_[trackIndex].DST;
		}
		else
		{
			DS = dg_CJUMP0_[trackIndex].DST * d_rand2_(3.0e0);
			dg_CJUMP0_[trackIndex].DSR = dg_CJUMP0_[trackIndex].DST - DS;
			if (dg_CJUMP1_[trackIndex].KSOFTI == 1)
			{
				if (dg_CJUMP0_[trackIndex].DST < 1.0e-8)
				{
					dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_CJUMP0_[trackIndex].W1;
					dg_PENELOPE_mod_.DESOFT[trackIndex] = dg_PENELOPE_mod_.SSOFT[trackIndex] * dg_CJUMP0_[trackIndex].DST;
				}
				else
				{
					EDE0 = dg_CJUMP0_[trackIndex].W1 * dg_CJUMP0_[trackIndex].DST;
					VDE0 = dg_CJUMP0_[trackIndex].W2 * dg_CJUMP0_[trackIndex].DST;
					FSEDE = fmax(1.0e0 - dg_CPIMFP_.DW1PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
					FSVDE = fmax(1.0e0 - dg_CPIMFP_.DW2PL[dg_CEGRID_.KE[trackIndex] - 1][dg_TRACK_mod_[index].MAT - 1] * EDE0, 0.75e0);
					EDE = EDE0 * FSEDE;
					VDE = VDE0 * FSVDE;

					// Geração de valores aleatórios DE com média EDE e variância VDE.
					SIGMA = sqrt(VDE);
					if (SIGMA < 0.333333333e0 * EDE)
					{
						// Distribuição gaussiana truncada.
						dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + d_rndg32_() * SIGMA;
					}
					else
					{
						RU = d_rand2_(4.0e0);
						EDE2 = EDE * EDE;
						VDE3 = 3.0e0 * VDE;
						if (EDE2 < VDE3)
						{
							PNULL = (VDE3 - EDE2) / (VDE3 + 3.0e0 * EDE2);
							if (RU < PNULL)
							{
								dg_PENELOPE_mod_.DESOFT[trackIndex] = 0.0e0;
								dg_PENELOPE_mod_.SSOFT[trackIndex] = 0.0e0;
								if (dg_CJUMP1_[trackIndex].KSOFTE == 0)
								{
									dg_CJUMP1_[trackIndex].MHINGE = 1;
									DS = dg_CJUMP0_[trackIndex].DST;
								}
								else
								{
									dg_CJUMP1_[trackIndex].KSOFTI = 0;
								}
								return;
							}
							else
							{
								// Distribuição Uniforme
								dg_PENELOPE_mod_.DESOFT[trackIndex] = 1.5e0 * (EDE + VDE / EDE) * (RU - PNULL) / (1.0e0 - PNULL);
							}
						}
						else
						{
							dg_PENELOPE_mod_.DESOFT[trackIndex] = EDE + (2.0e0 * RU - 1.0e0) * sqrt(VDE3);
						}
					}
					dg_PENELOPE_mod_.SSOFT[trackIndex] = dg_PENELOPE_mod_.DESOFT[trackIndex] / dg_CJUMP0_[trackIndex].DST;
				}
			}
		}
		return;
	}
	else
	{ // Fotons

		if (dg_TRACK_mod_[index].E < dg_CJUMP1_[trackIndex].ELAST1)
		{
			dg_CEGRID_.XEL[trackIndex] = log(dg_TRACK_mod_[index].E);
			dg_CEGRID_.XE[trackIndex] = 1.0e0 + (dg_CEGRID_.XEL[trackIndex] - dg_CEGRID_.DLEMP1) * dg_CEGRID_.DLFC;
			dg_CEGRID_.KE[trackIndex] = (int)dg_CEGRID_.XE[trackIndex];
			dg_CEGRID_.XEK[trackIndex] = dg_CEGRID_.XE[trackIndex] - dg_CEGRID_.KE[trackIndex];
			d_gimfp2_();
			dg_CJUMP1_[trackIndex].ELAST1 = dg_TRACK_mod_[index].E;
			dg_CJUMP0_[trackIndex].ST = dg_CJUMP0_[trackIndex].P[1 - 1] + dg_CJUMP0_[trackIndex].P[2 - 1] + dg_CJUMP0_[trackIndex].P[3 - 1] + dg_CJUMP0_[trackIndex].P[4 - 1] + dg_CJUMP0_[trackIndex].P[8 - 1];
		}

		DS = -log(d_rand2_(1.0e0)) / dg_CJUMP0_[trackIndex].ST;
	}

	/*	if (imprimiu==0){
			printf("\n\nJUMP2\n\n");
			imprimiu++;
		}*/
}

__device__ double d_rndg32_()
{

	/*
	Esta função entrega valores aleatórios no intervalo (-3.0,3.0)
	amostrado de uma distribuição gaussiana truncada que tem média zero e
	variação da unidade. A amostragem é realizada pelo método RITA.
	*/

	/*if (imprimiu==0){
		printf("\n\nRNDG32\n\n");
		imprimiu++;
	}*/

	// Selection of the interval (Walker's aliasing).
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;
	double RN, TST, RR, D, resultado;
	int K, I;

	RN = d_rand2_(1.0e0) * dg_CRNDG3_.NPM1 + 1.0e0;
	K = int(RN);
	TST = RN - K;
	if (TST < dg_CRNDG3_.F[K - 1])
	{
		I = K;
		RR = TST;
		D = dg_CRNDG3_.F[K - 1];
	}
	else
	{
		I = dg_CRNDG3_.KA[K - 1];
		RR = TST - dg_CRNDG3_.F[K - 1];
		D = 1.0e0 - dg_CRNDG3_.F[K - 1];
	}

	// Amostragem da distribuição cumulativa inversa racional.
	if (RR > 1.0e-12)
	{
		resultado = dg_CRNDG3_.X[I - 1] + ((1.0e0 + dg_CRNDG3_.A[I - 1] + dg_CRNDG3_.B[I - 1]) * D * RR / (D * D + (dg_CRNDG3_.A[I - 1] * D + dg_CRNDG3_.B[I - 1] * RR) * RR)) * (dg_CRNDG3_.X[I + 1 - 1] - dg_CRNDG3_.X[I - 1]);
	}
	else
	{
		resultado = dg_CRNDG3_.X[I - 1] + d_rand2_(2.0e0) * (dg_CRNDG3_.X[I + 1 - 1] - dg_CRNDG3_.X[I - 1]);
	}

	return resultado;
}

__device__ void d_gcone2_(double &UF, double &VF, double &WF)
{

	/*
	Esta sub-rotina amostra uma direção aleatória uniformemente dentro de um cone
	 com eixo central na direção (THETA,PHI) e abertura ALPHA.
	Os parâmetros são inicializados chamando a sub-rotina GCONE0.
	*/
	//	int index = blockDim.x * blockIdx.x + threadIdx.x;

	double WT, DF, SUV, UT, VT, DXY, DXYZ, FNORM;

	// Defina uma direção relativa ao eixo z.
	WT = dg_CGCONE_.CAPER + (1.0e0 - dg_CGCONE_.CAPER) * d_rand2_(1.0e0);
	DF = d_TWOPI * d_rand2_(2.0e0);
	SUV = sqrt(1.0e0 - WT * WT);
	UT = SUV * cos(DF);
	VT = SUV * sin(DF);

	// Rotacao para a direção do eixo do feixe

	UF = dg_CGCONE_.CPCT * UT - dg_CGCONE_.SPHI * VT + dg_CGCONE_.CPST * WT;
	VF = dg_CGCONE_.SPCT * UT + dg_CGCONE_.CPHI * VT + dg_CGCONE_.SPST * WT;
	WF = -dg_CGCONE_.STHE * UT + dg_CGCONE_.CTHE * WT;

	// Normalizaçao
	DXY = UF * UF + VF * VF;
	DXYZ = DXY + WF * WF;
	if (fabs(DXYZ - 1.0e0) > 1.0e-14)
	{
		FNORM = 1.0e0 / sqrt(DXYZ);
		UF = FNORM * UF;
		VF = FNORM * VF;
		WF = FNORM * WF;
	}
	/*if (imprimiu==0){
		printf("\n\ngcone2\n\n");
		imprimiu++;
	}*/
}

__device__ void d_simdet2_(int &N, int &ID)
{

	/*Calcula espectros de detectores de impacto, grava e carrega arquivos de despejo,
	acumula arquivos de despejo de diferentes execuções e grava os resultados.*/
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int IE, IT;

	// Espectro de energia das partículas que entram.

	if (dg_CIMDET_.LLE[ID - 1] == 1)
	{
		IE = (int)(1.0e0 + (log(dg_TRACK_mod_[index].E) - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}
	else
	{
		IE = (int)(1.0e0 + (dg_TRACK_mod_[index].E - dg_CIMDET_.EL[ID - 1]) * dg_CIMDET_.RBSE[ID - 1]);
	}

	if (N != dg_CIMDET_.LEDEP[ID - 1])
	{
		dg_CIMDET_.EDEP[ID - 1] = dg_CIMDET_.EDEP[ID - 1] + dg_CIMDET_.EDEPP[ID - 1];
		dg_CIMDET_.EDEP2[ID - 1] = dg_CIMDET_.EDEP2[ID - 1] + pow(dg_CIMDET_.EDEPP[ID - 1], 2);
		dg_CIMDET_.EDEPP[ID - 1] = dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
		dg_CIMDET_.LEDEP[ID - 1] = N;
	}
	else
	{
		dg_CIMDET_.EDEPP[ID - 1] = dg_CIMDET_.EDEPP[ID - 1] + dg_TRACK_mod_[index].E * dg_TRACK_mod_[index].WGHT;
	}

	if ((IE > 0) && (IE <= dg_CIMDET_.NE[ID - 1]))
	{
		if (N != dg_CIMDET_.LDIT[IE - 1][ID - 1])
		{
			dg_CIMDET_.DIT[IE - 1][ID - 1] = dg_CIMDET_.DIT[IE - 1][ID - 1] + dg_CIMDET_.DITP[IE - 1][ID - 1];
			dg_CIMDET_.DIT2[IE - 1][ID - 1] = dg_CIMDET_.DIT2[IE - 1][ID - 1] + pow(dg_CIMDET_.DITP[IE - 1][ID - 1], 2);
			dg_CIMDET_.DITP[IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT;
			dg_CIMDET_.LDIT[IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.DITP[IE - 1][ID - 1] = dg_CIMDET_.DITP[IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT;
		}

		if (N != dg_CIMDET_.LDIP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1])
		{
			dg_CIMDET_.DIP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.DIP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_CIMDET_.DIPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1];
			dg_CIMDET_.DIP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.DIP2[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + pow(dg_CIMDET_.DIPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1], 2);
			dg_CIMDET_.DIPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_TRACK_mod_[index].WGHT;
			dg_CIMDET_.LDIP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = N;
		}
		else
		{
			dg_CIMDET_.DIPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] = dg_CIMDET_.DIPP[dg_TRACK_mod_[index].KPAR - 1][IE - 1][ID - 1] + dg_TRACK_mod_[index].WGHT;
		}

		// Distribuição de Idades das particulas
		if (dg_CIMDET_.NAGE[ID - 1] > 0)
		{
			if (dg_CIMDET_.LLAGE[ID - 1] == 1)
			{
				IT = (int)(1.0e0 + (log(dg_TRACK_mod_[index].PAGE) - dg_CIMDET_.AGEL[ID - 1]) * dg_CIMDET_.RBAGE[ID - 1]);
			}
			else
			{
				IT = (int)(1.0e0 + (dg_TRACK_mod_[index].PAGE - dg_CIMDET_.AGEL[ID - 1]) * dg_CIMDET_.RBAGE[ID - 1]);
			}
			if ((IT > 0) && (IT <= dg_CIMDET_.NAGE[ID - 1]))
			{
				if (N != dg_CIMDET_.LAGEA[IT - 1][ID - 1])
				{
					dg_CIMDET_.AGE[IT - 1][ID - 1] = dg_CIMDET_.AGE[IT - 1][ID - 1] + dg_CIMDET_.AGEP[IT - 1][ID - 1];
					dg_CIMDET_.AGE2[IT - 1][ID - 1] = dg_CIMDET_.AGE2[IT - 1][ID - 1] + pow(dg_CIMDET_.AGEP[IT - 1][ID - 1], 2);
					dg_CIMDET_.AGEP[IT - 1][ID - 1] = dg_TRACK_mod_[index].WGHT;
					dg_CIMDET_.LAGEA[IT - 1][ID - 1] = N;
				}
				else
				{
					dg_CIMDET_.AGEP[IT - 1][ID - 1] = dg_CIMDET_.AGEP[IT - 1][ID - 1] + dg_TRACK_mod_[index].WGHT;
				}
			}
		}
	}
}

__device__ void d_sdose2_(double &DEP, double &XD, double &YD, double &ZD, int &MATC, int &N)
{

	/*
		Registra a distribuição de dose dentro da caixa de dose, grava e carrega
		despeja arquivos, acumula arquivos de despejo de diferentes execuções e grava
		Resultados .
	*/
	// int index = blockDim.x * blockIdx.x + threadIdx.x;
	int I1, I2, I3;
	double RD;

	// Distribuição de Dose
	// printf("\nKDOSE: %d\n", dg_CDOSE1_.KDOSE);
	if (dg_CDOSE1_.KDOSE == 1)
	{ // Caixa
		if ((ZD > dg_CDOSE3_.DXL[3 - 1]) && (ZD < dg_CDOSE3_.DXU[3 - 1]))
		{
			I3 = (int)(1.0e0 + (ZD - dg_CDOSE3_.DXL[3 - 1]) * dg_CDOSE3_.RBDOSE[3 - 1]);
			if (N != dg_CDOSE2_.LDDOSE[I3 - 1])
			{
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_CDOSE2_.DDOSE[I3 - 1], dg_CDOSE2_.DDOSEP[I3 - 1]);
					atomicAdd2(&dg_CDOSE2_.DDOSE2[I3 - 1], pow(dg_CDOSE2_.DDOSEP[I3 - 1], 2));
				}
				else
				{
					dg_CDOSE2_.DDOSE[I3 - 1] = dg_CDOSE2_.DDOSE[I3 - 1] + dg_CDOSE2_.DDOSEP[I3 - 1];
					dg_CDOSE2_.DDOSE2[I3 - 1] = dg_CDOSE2_.DDOSE2[I3 - 1] + pow(dg_CDOSE2_.DDOSEP[I3 - 1], 2);
				}
				dg_CDOSE2_.DDOSEP[I3 - 1] = DEP * dg_PENELOPE_mod_.RDEN[MATC - 1];
				dg_CDOSE2_.LDDOSE[I3 - 1] = N;
			}
			else
			{
				if (d_cudaCapability >= 6)
				{
					atomicAdd2(&dg_CDOSE2_.DDOSEP[I3 - 1], DEP * dg_PENELOPE_mod_.RDEN[MATC - 1]);
				}
				else
				{
					dg_CDOSE2_.DDOSEP[I3 - 1] = dg_CDOSE2_.DDOSEP[I3 - 1] + DEP * dg_PENELOPE_mod_.RDEN[MATC - 1];
				}
			}

			if (((XD > dg_CDOSE3_.DXL[1 - 1]) && (XD < dg_CDOSE3_.DXU[1 - 1])) && ((YD > dg_CDOSE3_.DXL[2 - 1]) && (YD < dg_CDOSE3_.DXU[2 - 1])))
			{
				I1 = (int)(1.0e0 + (XD - dg_CDOSE3_.DXL[1 - 1]) * dg_CDOSE3_.RBDOSE[1 - 1]);
				I2 = (int)(1.0e0 + (YD - dg_CDOSE3_.DXL[2 - 1]) * dg_CDOSE3_.RBDOSE[2 - 1]);
				if (N != dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1])
				{
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1], dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1]);
						atomicAdd2(&dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1], pow(dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1], 2));
					}
					else
					{
						dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] + dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1];
						dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] + pow(dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1], 2);
					}
					dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = DEP;
					dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1] = N;
				}
				else
				{
					if (d_cudaCapability >= 6)
					{
						atomicAdd2(&dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1], DEP);
					}
					else
					{
						dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] + DEP;
					}
				}
			}
		}
	}
	else if (dg_CDOSE1_.KDOSE == 2)
	{ // Cilindro
		if ((ZD > dg_CDOSE3_.DXL[3 - 1]) && (ZD < dg_CDOSE3_.DXU[3 - 1]))
		{
			I3 = (int)(1.0e0 + (ZD - dg_CDOSE3_.DXL[3 - 1]) * dg_CDOSE3_.RBDOSE[3 - 1]);
			if (N != dg_CDOSE2_.LDDOSE[I3 - 1])
			{
				dg_CDOSE2_.DDOSE[I3 - 1] = dg_CDOSE2_.DDOSE[I3 - 1] + dg_CDOSE2_.DDOSEP[I3 - 1];
				dg_CDOSE2_.DDOSE2[I3 - 1] = dg_CDOSE2_.DDOSE2[I3 - 1] + pow(dg_CDOSE2_.DDOSEP[I3 - 1], 2);
				dg_CDOSE2_.DDOSEP[I3 - 1] = DEP * dg_PENELOPE_mod_.RDEN[MATC - 1];
				dg_CDOSE2_.LDDOSE[I3 - 1] = N;
			}
			else
			{
				dg_CDOSE2_.DDOSEP[I3 - 1] = dg_CDOSE2_.DDOSEP[I3 - 1] + DEP * dg_PENELOPE_mod_.RDEN[MATC - 1];
			}

			RD = sqrt(XD * XD + YD * YD);
			if (RD < dg_CDOSE3_.DXU[1 - 1])
			{
				I1 = (int)(1.0e0 + RD * dg_CDOSE3_.RBDOSE[1 - 1]);
				I2 = 1;
				if (N != dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1])
				{
					dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] + dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1];
					dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] + pow(dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1], 2);
					dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = DEP;
					dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1] = N;
				}
				else
				{
					dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] + DEP;
				}
			}
		}
	}
	else
	{ // Esfera
		RD = sqrt(XD * XD + YD * YD + ZD * ZD);
		if (RD < dg_CDOSE3_.DXU[1 - 1])
		{
			I1 = (int)(1.0e0 + RD * dg_CDOSE3_.RBDOSE[1 - 1]);
			I2 = 1;
			I3 = 1;
			if (N != dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1])
			{
				dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE[I3 - 1][I2 - 1][I1 - 1] + dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1];
				dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSE2[I3 - 1][I2 - 1][I1 - 1] + pow(dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1], 2);
				dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = DEP;
				dg_CDOSE1_.LDOSE[I3 - 1][I2 - 1][I1 - 1] = N;
			}
			else
			{
				dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] = dg_CDOSE1_.DOSEP[I3 - 1][I2 - 1][I1 - 1] + DEP;
			}
		}
	}
}

__device__ void d_cleans2_()
{

	// Esta sub-rotina inicializa a pilha secundária. Deve ser chamada antes de iniciar a simulação de cada pista primária.

	dg_SECST_.NSEC = 0;
	dg_nTRACKS_.nPRITRACK = 0;
	dg_nTRACKS_.nSECTRACK_E = 0;
	dg_nTRACKS_.nSECTRACK_G = 0;
	dg_nTRACKS_.nSECTRACK_P = 0;
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

/*__device__ int atomicAdd2(int* address, int val)
{
 unsigned long long int* address_as_ull =
 (unsigned long long int*)address;
 unsigned long long int old = *address_as_ull, assumed;
 do {
 /*assumed = old;
 old = atomicCAS(address_as_ull, assumed,
 __int_as_longlong(val +
 __longlong_as_int(assumed)));
 // Note: uses integer comparison to avoid hang in case of NaN (since NaN !=NaN)
 } while (assumed != old);
 return __longlong_as_int(old);
}*/

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

__global__ void g_cpyTrack_simular(int tipo, int size){

	int index =  blockDim.x * blockIdx.x + threadIdx.x;

	if (index < size){
		if (tipo == 1)
		{ // eletrons
			dg_TRACK_mod_[index].E = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].E;
			dg_TRACK_mod_[index].X = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].X;
			dg_TRACK_mod_[index].Y = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].Y;
			dg_TRACK_mod_[index].Z = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].Z;
			dg_TRACK_mod_[index].U = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].U;
			dg_TRACK_mod_[index].V = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].V;
			dg_TRACK_mod_[index].W = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].W;
			dg_TRACK_mod_[index].WGHT = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].WGHT;
			dg_TRACK_mod_[index].KPAR = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].KPAR;
			dg_TRACK_mod_[index].IBODY = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].IBODY;
			dg_TRACK_mod_[index].MAT = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].MAT;

			dg_TRACK_mod_[index].ILB[1 - 1] = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].ILB[1 - 1];
			dg_TRACK_mod_[index].ILB[2 - 1] = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].ILB[2 - 1];
			dg_TRACK_mod_[index].ILB[3 - 1] = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].ILB[3 - 1];
			dg_TRACK_mod_[index].ILB[4 - 1] = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].ILB[4 - 1];
			dg_TRACK_mod_[index].ILB[5 - 1] = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].ILB[5 - 1];

			dg_TRACK_mod_[index].N = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].N;
			dg_TRACK_mod_[index].INDEX = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].INDEX;
			dg_TRACK_mod_[index].IEXIT = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].IEXIT;
			dg_TRACK_mod_[index].STEP = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].STEP;

			if (d_wIPOLI == 1)
			{
				dg_TRACK_mod_[index].SP1 = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].SP1;
				dg_TRACK_mod_[index].SP2 =dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].SP2;
				dg_TRACK_mod_[index].SP3 = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].SP3;
				dg_TRACK_mod_[index].IPOL = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].IPOL;
			}
			else
			{
				dg_TRACK_mod_[index].SP1 = 0.0e0;
				dg_TRACK_mod_[index].SP2 = 0.0e0;
				dg_TRACK_mod_[index].SP3 = 0.0e0;
				dg_TRACK_mod_[index].IPOL = 0;
			}
			dg_TRACK_mod_[index].PAGE = dg_SECTRACK_E_[(dg_nTRACKS_.nSECTRACK_E-1) - index].PAGE;


		}
	}
}

__global__ void g_cpyTrack_complementar(int tipo, int size){

	int index =  blockDim.x * blockIdx.x + threadIdx.x;
	int i = 0;

	if ((index < size) && (dg_TRACK_mod_[index].STEP == 9)){
		if (tipo == 1){
		
				//copia para o vetor contabilizador novo
				//printf("%d\n\n\n\n", atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, -1));
				i = atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, -1);
				//printf(" Energia da particula sec eletron  %f\n\n\n",dg_SECTRACK_E_[atomicAdd2(&dg_nTRACKS_.nSECTRACK_E, -1)].E );
				atomicAdd2(&dg_nTRACKS_.nFINISH, 1);
				dg_TRACK_mod_[index].E = dg_SECTRACK_E_[i].E;
				dg_TRACK_mod_[index].X = dg_SECTRACK_E_[i].X;
				dg_TRACK_mod_[index].Y = dg_SECTRACK_E_[i].Y;
				dg_TRACK_mod_[index].Z = dg_SECTRACK_E_[i].Z;
				dg_TRACK_mod_[index].U = dg_SECTRACK_E_[i].U;
				dg_TRACK_mod_[index].V = dg_SECTRACK_E_[i].V;
				dg_TRACK_mod_[index].W = dg_SECTRACK_E_[i].W;
				dg_TRACK_mod_[index].WGHT = dg_SECTRACK_E_[i].WGHT;
				dg_TRACK_mod_[index].KPAR = dg_SECTRACK_E_[i].KPAR;
				dg_TRACK_mod_[index].IBODY = dg_SECTRACK_E_[i].IBODY;
				dg_TRACK_mod_[index].MAT = dg_SECTRACK_E_[i].MAT;

				dg_TRACK_mod_[index].ILB[1 - 1] = dg_SECTRACK_E_[i].ILB[1 - 1];
				dg_TRACK_mod_[index].ILB[2 - 1] = dg_SECTRACK_E_[i].ILB[2 - 1];
				dg_TRACK_mod_[index].ILB[3 - 1] = dg_SECTRACK_E_[i].ILB[3 - 1];
				dg_TRACK_mod_[index].ILB[4 - 1] = dg_SECTRACK_E_[i].ILB[4 - 1];
				dg_TRACK_mod_[index].ILB[5 - 1] = dg_SECTRACK_E_[i].ILB[5 - 1];

				dg_TRACK_mod_[index].N = dg_SECTRACK_E_[i].N;
				dg_TRACK_mod_[index].INDEX = dg_SECTRACK_E_[i].INDEX;
				dg_TRACK_mod_[index].IEXIT = dg_SECTRACK_E_[i].IEXIT;
				dg_TRACK_mod_[index].STEP = dg_SECTRACK_E_[i].STEP;


				if (d_wIPOLI == 1)
				{
					dg_TRACK_mod_[index].SP1 = dg_SECTRACK_E_[i].SP1;
					dg_TRACK_mod_[index].SP2 =dg_SECTRACK_E_[i].SP2;
					dg_TRACK_mod_[index].SP3 = dg_SECTRACK_E_[i].SP3;
					dg_TRACK_mod_[index].IPOL = dg_SECTRACK_E_[i].IPOL;
				}
				else
				{
					dg_TRACK_mod_[index].SP1 = 0.0e0;
					dg_TRACK_mod_[index].SP2 = 0.0e0;
					dg_TRACK_mod_[index].SP3 = 0.0e0;
					dg_TRACK_mod_[index].IPOL = 0;
				}
				dg_TRACK_mod_[index].PAGE = dg_SECTRACK_E_[i].PAGE;
			
		}

	}
}


__global__ void g_imprimir(int size){

	int index =  blockDim.x * blockIdx.x + threadIdx.x;

	if (index < size){
		printf("Index %d  Energia particula E: %.18lf\n", index, dg_SECTRACK_E_[index].E);
	}


}



