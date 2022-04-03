extern "C"{

	__device__ void d_gcone02_(double THETA, double PHI, double ALPHA);

	__device__ void d_gcone2_(double& UF, double& VF, double& WF);

	__device__ void d_step2_(double& DS, double& DSEF, int& NCROSS);

	__device__ void d_steplb2_(int& KB, int& IERR);

	__device__ void d_stepsi2_(int& KB, int& NSC);

	__device__ void d_fsurf2_(int& KS, double& A, double& B, double& C);

	__device__ void d_locate2_();

    __global__ void teste();


}




__device__ void d_step2_(double& DS, double& DSEF, int& NCROSS) {
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

	DSEF = 0.0e0;
	dg_PENGEOM_mod_.DSTOT = 0.0e0;
	NCROSS = 0;
	dg_PENGEOM_mod_.KSLAST = 0;
	double DSRES;
	int KB1;
	//double d_S[NS2M];
	//int d_ISNS2M];

	//double *S = (double *)malloc(NS2M*sizeof(double));
	//int *IS = (int *)malloc(NS2M*sizeof(int));


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

	int NSC = 0; //Número de cruzamentos da superfície à frente da partícula.
	int NSCT;
	int MAT0;

	//	printf("1\n");

	//	printf("LVERB %d\n", dg_PENGEOM_mod_.LVERB);
	//	printf(dg_PENGEOM_mod_.LVERB ? "true\n" : "false\n")
	for (int I = 1; I <= dg_QSURF_.NSURF; I++) {
		dg_QTREE_.KSP[I - 1] = 0; //Ponteiros laterais das superfícies avaliadas.
	}

	//	printf("2\n");


	MAT0 = dg_TRACK_mod_.MAT; //Material Inicial

	if (dg_TRACK_mod_.MAT == 0) {
		DSRES = 1.0e35; //No vácuo, as partículas voam livremente.
	}
	else {
		DSRES = DS; // comprimento do camimho residual
	}

	//		printf("3\n");
		//A partícula entra de fora do recinto.

	if (dg_TRACK_mod_.IBODY > dg_QTREE_.NBODYS) {
		KB1 = dg_QTREE_.NBODYS;
		//d_stepsi2_(KB1, S, IS, NSC);
		d_stepsi2_(KB1, NSC);
		if (NSC == 0)
			goto L300;
		NSCT = NSC;
		NST = dg_QTREE_.KSURF[NXG - 1][KB1 - 1];

		for (int KI = NSCT; KI >= 1; KI--) {
			// a particula atravessa uma superficie
			dg_PENGEOM_mod_.KSLAST = d_IS[KI - 1];
			if (dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] == 1)
				dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 2;
			else
				dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 1;

			DSP = d_S[KI - 1];
			DSEF = DSEF + DSP;
			dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSP;
			dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSP * dg_TRACK_mod_.U;
			dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSP * dg_TRACK_mod_.V;
			dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSP * dg_TRACK_mod_.W;
			NSC = NSC - 1;

			if (NSC > 0) {
				for (int I = 1; I <= NSC; I++) {
					d_S[I - 1] = d_S[I - 1] - DSP;
				}
			}

			for (int KSS = 1; KSS <= NST; KSS++) {
				KS1 = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
				KF = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
				if ((KF < 3) && (dg_QTREE_.KSP[KS1 - 1] != KF))
					goto L101;
			}
			// A partícula entra no invólucro.
		L100:;
			d_steplb2_(KB1, IERR);
			// A partícula entra em um submódulo.
			if (IERR == -1) {
				KB1 = dg_TRACK_mod_.IBODY;
				//d_stepsi2_(KB1, S, IS, NSC);
				d_stepsi2_(KB1, NSC);
				goto L100;
			}
			else {
				//A particula entrou em um corpo material
				if (dg_TRACK_mod_.MAT != 0) {
					NCROSS = 1;
					//				free(S);
					//				free(IS);
					return;
				}
				else {
					KB1 = dg_TRACK_mod_.IBODY;
					//d_stepsi2_(KB1, S, IS, NSC);
					d_stepsi2_(KB1, NSC);
					goto L200;
				}

			}

			//Neste ponto, o programa saiu do ciclo DO.
		L101:;
		}
		//	printf("6\n");
		goto L300;
	}

	//Cruzamentos de superfície.

	IBODYL = dg_TRACK_mod_.IBODY;
	if (dg_PENGEOM_mod_.LVERB)
		NERR = 0;
L102:;
	KB1 = dg_TRACK_mod_.IBODY;
	//d_stepsi2_(KB1, S, IS, NSC);
	d_stepsi2_(KB1, NSC);
	d_steplb2_(KB1, IERR);

	//Evidência de erros de arredondamento.
	if (IERR != 0) {
		if (NSC > 0) {
			//Quando uma superfície está muito próxima, movemos a partícula além dela.
			if (d_S[NSC - 1] < 1e-10) {
				dg_PENGEOM_mod_.KSLAST = d_IS[NSC - 1];
				if (dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] == 1) {
					dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 2;
				}
				else {
					dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 1;
				}

				DSP = d_S[NSC - 1];
				dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSP * dg_TRACK_mod_.U;
				dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSP * dg_TRACK_mod_.V;
				dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSP * dg_TRACK_mod_.W;
				if (dg_TRACK_mod_.MAT == MAT0) {
					DSEF = DSEF + DSP;
					DSRES = DSRES - DSP;
				}

				dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSP;
				NSC = NSC - 1;

				if (dg_TRACK_mod_.IBODY <= dg_QTREE_.NBODYS)
					goto L102;
			}
		}
		if (dg_PENGEOM_mod_.LVERB) {

			NERR = NERR + 1;
			if ((dg_QTREE_.NWARN < 100) && (dg_TRACK_mod_.MAT != 0)) {

				printf("WARNING, STEP: Accidental undershot or r");

				for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB1 - 1]; KSS++) {
					KS = dg_QTREE_.KSURF[KSS - 1][KB1 - 1];
					KFLO = dg_QTREE_.KFLAG[KSS - 1][KB1 - 1];
					if (KFLO < 3) {
						for (int KI = NSC; KI >= 1; KI--) {
							if (KS == d_IS[KI - 1]) {
								SW = d_S[KI - 1];
								goto L103;
							}
						}

						SW = 0.0e0;
					L103:;
						d_fsurf2_(KS, A, B, C);
						if (KFLO == dg_QTREE_.KSP[KS - 1]) {
							printf("KS, KFLO, KSP, SW %lf", SW);
						}
						else {
							printf("KS, KFLO, KSP, SW %lf", SW);
							dg_PENGEOM_mod_.KSLAST = KS;
						}
					}
				}
				dg_QTREE_.NWARN = dg_QTREE_.NWARN + 1;
			}

		}
		if (dg_TRACK_mod_.IBODY <= dg_QTREE_.NBODYS)
			goto L102;
	}

	if ((dg_PENGEOM_mod_.KDET[dg_TRACK_mod_.IBODY - 1] != dg_PENGEOM_mod_.KDET[IBODYL - 1]) || (dg_TRACK_mod_.MAT != MAT0)) {
		NCROSS = 1;
		DSEF = 0.0e0;
		//		free(S);
		//		free(IS);
		return;
	}

	//A particula permanece no mesmo material

	if ((dg_TRACK_mod_.MAT != 0) && (DSRES < d_S[NSC - 1])) {
		if (dg_TRACK_mod_.MAT == MAT0)
			DSEF = DSEF + DSRES;
		dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSRES;
		dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSRES * dg_TRACK_mod_.U;
		dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSRES * dg_TRACK_mod_.V;
		dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSRES * dg_TRACK_mod_.W;
		//		free(S);
		//		free(IS);
		return;
	}


	// Nova posição

L200:;
	if (NSC == 0) {
		if (dg_TRACK_mod_.MAT == MAT0)
			DSEF = DSEF + DSRES;
		dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSRES;
		dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSRES * dg_TRACK_mod_.U;
		dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSRES * dg_TRACK_mod_.V;
		dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSRES * dg_TRACK_mod_.W;
		return;
	}
	NSCT = NSC;
	MATL = dg_TRACK_mod_.MAT;
	IBODYL = dg_TRACK_mod_.IBODY;
	for (int KI = NSCT; KI >= 1; KI--) {
		//A etapa termina dentro do corpo
		if (DSRES < d_S[KI - 1]) {
			if (dg_TRACK_mod_.MAT == MAT0)
				DSEF = DSEF + DSRES;
			dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSRES;
			dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSRES * dg_TRACK_mod_.U;
			dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSRES * dg_TRACK_mod_.V;
			dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSRES * dg_TRACK_mod_.W;
			//			free(S);
			//			free(IS);
			return;
		}

		// A particula atravessa uma superfice
		dg_PENGEOM_mod_.KSLAST = d_IS[KI - 1];
		if (dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] == 1)
			dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 2;
		else
			dg_QTREE_.KSP[dg_PENGEOM_mod_.KSLAST - 1] = 1;

		DSP = d_S[KI - 1];
		dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSP * dg_TRACK_mod_.U;
		dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSP * dg_TRACK_mod_.V;
		dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSP * dg_TRACK_mod_.W;
		if (dg_TRACK_mod_.MAT == MAT0) {
			DSEF = DSEF + DSP;
			DSRES = DSRES - DSP;
		}
		dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSP;
		NSC = NSC - 1;
		if (NSC > 0) {
			for (int I = 1; I <= NSC; I++) {
				d_S[I - 1] = d_S[I - 1] - DSP;
			}
		}
		d_steplb2_(KB1, IERR);

	L201:;
		KB1 = dg_TRACK_mod_.IBODY;
		if (IERR == -1) {
			// A particula entrou em um submodulo
			//d_stepsi2_(KB1, S, IS, NSC);
			d_stepsi2_(KB1, NSC);
			d_steplb2_(KB1, IERR);
			goto L201;
		}
		else if (IERR == 1) {
			//A partícula deixa o corpo ou módulo.
			if (dg_TRACK_mod_.IBODY <= dg_QTREE_.NBODYS) {
				//d_stepsi2_(KB1, S, IS, NSC);
				d_stepsi2_(KB1, NSC);
				d_steplb2_(KB1, IERR);
				goto L201;
			}
			else {
				//A partícula sai do recinto.
				if (dg_TRACK_mod_.MAT != MATL)
					NCROSS = NCROSS + 1;
				goto L300;
			}
		}

		//A partícula continua voando quando entra em uma região vazia 
		if (dg_TRACK_mod_.MAT == 0) {
			if (MATL == MAT0)
				NCROSS = NCROSS + 1;
			MATL = 0;
			DSRES = 1.0e35;
			goto L202;
			//A partícula continua voando quando entra em um novo corpo do 
			//mesmo material que não faz parte de um detector diferente ...
		}
		else if (dg_TRACK_mod_.MAT == MATL) {
			if (dg_PENGEOM_mod_.KDET[dg_TRACK_mod_.IBODY - 1] == dg_PENGEOM_mod_.KDET[IBODYL - 1]) {
				goto L202;
			}
			else {
				NCROSS = NCROSS + 1;
				return;
			}
			//.. e para quando penetra um novo corpo material ou umDetector
		}
		else {
			NCROSS = NCROSS + 1;
			//			free(S);
			//			free(IS);
			return;
		}
	L202:;
		//d_stepsi2_(KB1, S, IS, NSC);
		d_stepsi2_(KB1, NSC);
		goto L200;
		//Neste ponto, o programa saiu do ciclo DO.	
	//L203:; indica o final do loop
	}
	//A particula sai do recinto.
L300:;

	DSP = 1.0e36;
	dg_TRACK_mod_.IBODY = dg_QTREE_.NBODYS + 1;
	dg_TRACK_mod_.MAT = 0;
	if (dg_TRACK_mod_.MAT == MAT0)
		DSEF = DSEF + DSP;
	dg_PENGEOM_mod_.DSTOT = dg_PENGEOM_mod_.DSTOT + DSP;
	dg_TRACK_mod_.X = dg_TRACK_mod_.X + DSP * dg_TRACK_mod_.U;
	dg_TRACK_mod_.Y = dg_TRACK_mod_.Y + DSP * dg_TRACK_mod_.V;
	dg_TRACK_mod_.Z = dg_TRACK_mod_.Z + DSP * dg_TRACK_mod_.W;
	//	free(S);
	//	free(IS);
}


__device__ void d_stepsi2_(int& KB, int& NSC) {
	/*Calcula as interseções da trajetória com o limite
		Superfícies do corpo KB. Os cruzamentos são adicionados à lista e
		classificados em ordem decrescente.
		Esta sub-rotina funciona apenas quando chamada de dentro da sub-rotina STEP.*/

		//stepsi_(KB, S, IS, NSC);
			//return;

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



	//Determine cruzamentos de superfície.


	for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++) {
		//printf("\nKSURF: %d\n",  dg_QTREE_.KSURF[NXG - 1][KB -1]);

		/*As interseções com uma determinada superfície são calculadas apenas uma vez.
		O ponteiro lateral de uma superfície deve ser alterado cada vez que o a superfície está cruzada.*/

		KFL = dg_QTREE_.KFLAG[KSS - 1][KB - 1];
		if (KFL > 4)
			goto L100;
		KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
		if (dg_QTREE_.KSP[KS - 1] != 0)
			goto L100;
		d_fsurf2_(KS, A, B, C);
		ABSA = fabs(A);
		ABSB = fabs(B);

		// Plano, unica raiz
		if (ABSA < 1.0e-36) {
			if (ABSB > 0.0e0) {
				if (C < -FUZZL) {
					dg_QTREE_.KSP[KS - 1] = 1;
				}
				else if (C > FUZZL) {
					dg_QTREE_.KSP[KS - 1] = 2;
				}
				else {
					if (B < 0.0e0) {
						dg_QTREE_.KSP[KS - 1] = 1;
					}
					else {
						dg_QTREE_.KSP[KS - 1] = 2;
					}
					goto L100;
				}
				T1 = -C / B;
				if (T1 > 0.0e0) {
					NSC = NSC + 1;
					d_IS[NSC - 1] = KS;
					d_S[NSC - 1] = T1;
				}
			}
			else {
				if (C < 0.0e0) {
					dg_QTREE_.KSP[KS - 1] = 1;
				}
				else {
					dg_QTREE_.KSP[KS - 1] = 2;
				}
			}

			// Superficie não plana, duas raizes			
		}
		else {
			DISCR = B * B - 4.0e0 * A * C;
			FUZZ = FUZZL * DISCR / ABSA;
			if (C < -FUZZ) {
				IAMBIG = 0;
				dg_QTREE_.KSP[KS - 1] = 1;
			}
			else if (C > FUZZ) {
				IAMBIG = 0;
				dg_QTREE_.KSP[KS - 1] = 2;
			}
			else {
				IAMBIG = 1;
				if (B < 0.0e0) {
					dg_QTREE_.KSP[KS - 1] = 1;
				}
				else {
					dg_QTREE_.KSP[KS - 1] = 2;
				}
			}

			if (DISCR < 1.0e-36)
				goto L100;

			if (IAMBIG == 0) {
				R2A = 0.5e0 / A;
				DELTA = sqrt(DISCR) * fabs(R2A);
				SH = -B * R2A;
				T1 = SH - DELTA;
				if (T1 > 0.0e0) {
					NSC = NSC + 1;
					d_IS[NSC - 1] = KS;
					d_S[NSC - 1] = T1;
				}
				T2 = SH + DELTA;
				if (T2 > 0.0e0) {
					NSC = NSC + 1;
					d_IS[NSC - 1] = KS;
					d_S[NSC - 1] = T2;
				}
			}
			else {
				if (B * A < 0.0e0) {
					R2A = 0.5e0 / A;
					DELTA = sqrt(DISCR) * fabs(R2A);
					SH = -B * R2A;
					T2 = SH + DELTA;
					NSC = NSC + 1;
					d_IS[NSC - 1] = KS;
					d_S[NSC - 1] = fmax(T2, 0.0e0);
				}
			}
		}
	L100:;
	}



	//Classifique as distâncias da superfície em ordem decrescente.
   // printf("\nNSC %d\n", NSC);
	if (NSC > 1) {
		for (int KI = 1; KI <= (NSC - 1); KI++) {
			SMAX = d_S[KI - 1];
			KMAX = KI;
			for (int KJ = (KI + 1); KJ <= NSC; KJ++) {
				if (d_S[KJ - 1] > SMAX) {
					SMAX = d_S[KJ - 1];
					KMAX = KJ;
				}
			}
			if (KMAX != KI) {
				SMAX = d_S[KI - 1];
				d_S[KI - 1] = d_S[KMAX - 1];
				d_S[KMAX - 1] = SMAX;
				KKMAX = d_IS[KI - 1];
				d_IS[KI - 1] = d_IS[KMAX - 1];
				d_IS[KMAX - 1] = KKMAX;
			}
		}
	}
}


__device__ void d_steplb2_(int& KB, int& IERR) {

	/*Ajuda a encontrar o corpo ou módulo que contém os ponteiros laterais fornecidos para as superfícies analisadas.
	 A subrotina STEPLB funciona apenas quando invocada de dentro da sub-rotina STEP.
	 Ele se move através da árvore de módulos a única etapa*/


	int NLBOD;
	int KBS;
	int KS;
	int KF;
	int KBD;

	//Analisa o corpo ou o modulo atual. 

	if (dg_QBODY_.KBOMO[KB - 1] == 0) {
		//Corpo
		NLBOD = dg_QBODY_.KBODY[NXG - 1][KB - 1];
		if (NLBOD > 0) {
			for (int KBB = 1; KBB <= NLBOD; KBB++) {
				KBS = dg_QBODY_.KBODY[KBB - 1][KB - 1];
				for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KBS - 1]; KSS++) {
					KS = dg_QTREE_.KSURF[KSS - 1][KBS - 1];
					KF = dg_QTREE_.KFLAG[KSS - 1][KBS - 1];
					if ((KF < 3) && (dg_QTREE_.KSP[KS - 1] != KF))
						goto L100;
				}
				dg_TRACK_mod_.IBODY = KBS;
				if (dg_QTREE_.KDGHT[NXG - 1][dg_TRACK_mod_.IBODY - 1] > 1) {
					IERR = -1; //A particula está dentro de um módulo irmã
				}
				else {
					IERR = 0; //A particula está dentro de um corpo irma	
					dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_.IBODY - 1];
				}
				return;
			L100:;
			}
		}

		for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++) {
			KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
			KF = dg_QTREE_.KFLAG[KSS - 1][KB - 1];
			if ((KF < 3) && (dg_QTREE_.KSP[KS - 1] != KF))
				goto L300;
		}
		dg_TRACK_mod_.IBODY = KB;
		IERR = 0;
		dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_.IBODY - 1];
		return;
	}
	else {
		// Modulo
		for (int KBB = 1; KBB <= dg_QTREE_.KDGHT[NXG - 1][KB - 1]; KBB++) {
			KBD = dg_QTREE_.KDGHT[KBB - 1][KB - 1];
			for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KBD - 1]; KSS++) {
				KS = dg_QTREE_.KSURF[KSS - 1][KBD - 1];
				KF = dg_QTREE_.KFLAG[KSS - 1][KBD - 1];
				if ((KF < 3) && (dg_QTREE_.KSP[KS - 1] != KF))
					goto L200;
			}
			dg_TRACK_mod_.IBODY = KBD;
			if (KBD == KB) {
				IERR = 0; //A partícula permanece dentro do módulo atual.
				dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_.IBODY - 1];
			}
			else {
				if (dg_QTREE_.KDGHT[NXG - 1][KBD - 1] > 1) {
					IERR = -1; // A particula esta dentro de um submodulo
				}
				else {
					IERR = 0; //A particula esta dentro de um corpo simples;
					dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[dg_TRACK_mod_.IBODY - 1];
				}
			}
			return;
		L200:;
		}
	}

	//A partícula está fora do corpo ou módulo atual.
L300:;
	IERR = 1;
	dg_TRACK_mod_.IBODY = dg_QTREE_.KMOTH[KB - 1];
	if (dg_TRACK_mod_.IBODY == 0) {
		dg_TRACK_mod_.IBODY = dg_QTREE_.NBODYS + 1;
		dg_TRACK_mod_.MAT = 0;
	}
}

__device__ void d_fsurf2_(int& KS, double& A, double& B, double& C) {


	/*
	Calcula os par�metros da fun��o mestre da superf�cie KS e o raio (X, Y, Z) + S * (U, V, W).
	*/
	double XXX, YYY, ZZZ;



	if (dg_QSURF_.KPLANE[KS - 1] == 0) {
		A = dg_TRACK_mod_.U * (dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_.U + dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_.V + dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_.W) +
			dg_TRACK_mod_.V * (dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_.V + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_.W) + dg_TRACK_mod_.W * dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_.W;
		XXX = dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_.X + dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_.Y + dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_.Z + dg_QSURF_.AX[KS - 1];
		YYY = dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_.Y + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_.Z + dg_QSURF_.AY[KS - 1];
		ZZZ = dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_.Z + dg_QSURF_.AZ[KS - 1];

		B = dg_TRACK_mod_.U * (dg_QSURF_.AXX[KS - 1] * dg_TRACK_mod_.X + XXX) + dg_TRACK_mod_.V * (dg_QSURF_.AXY[KS - 1] * dg_TRACK_mod_.X + dg_QSURF_.AYY[KS - 1] * dg_TRACK_mod_.Y + YYY) +
			dg_TRACK_mod_.W * (dg_QSURF_.AXZ[KS - 1] * dg_TRACK_mod_.X + dg_QSURF_.AYZ[KS - 1] * dg_TRACK_mod_.Y + dg_QSURF_.AZZ[KS - 1] * dg_TRACK_mod_.Z + ZZZ);

		C = dg_TRACK_mod_.X * XXX + dg_TRACK_mod_.Y * YYY + dg_TRACK_mod_.Z * ZZZ + dg_QSURF_.A0[KS - 1];
	}
	else {
		A = 0.0e0;
		B = dg_TRACK_mod_.U * dg_QSURF_.AX[KS - 1] + dg_TRACK_mod_.V * dg_QSURF_.AY[KS - 1] + dg_TRACK_mod_.W * dg_QSURF_.AZ[KS - 1];
		C = dg_TRACK_mod_.X * dg_QSURF_.AX[KS - 1] + dg_TRACK_mod_.Y * dg_QSURF_.AY[KS - 1] + dg_TRACK_mod_.Z * dg_QSURF_.AZ[KS - 1] + dg_QSURF_.A0[KS - 1];
	}

}


__device__ void d_locate2_() {
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


	double A = 0.0;
	double B = 0.0;
	double C = 0.0;
	double FUZZ = 0.0;
	const double FUZZL = 1.0e-12;;
	int KS;
	int KB;
	int KF;
	double ABSA = 0.0;

	for (int I = 1; I <= dg_QSURF_.NSURF; I++) {
		dg_QTREE_.KSP[I - 1] = 0;
	}
	int KB0 = dg_QTREE_.NBODYS;

d100:
	for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB0 - 1]; KSS++) {
		KS = dg_QTREE_.KSURF[KSS - 1][KB0 - 1];
		if ((dg_QTREE_.KSP[KS - 1] != 0) || (dg_QTREE_.KFLAG[KSS - 1][KB0 - 1] > 4)) {
			goto d101;
		}

		d_fsurf2_(KS, A, B, C);

		ABSA = fabs(A);

		if (ABSA > 1.0e-36)
			FUZZ = FUZZL * (B * B - 4.0e0 * A * C) / ABSA;
		else
			FUZZ = FUZZL * fabs(B);

		if (C < (-FUZZ)) {
			dg_QTREE_.KSP[KS - 1] = 1;
		}
		else if (C > FUZZ) {
			dg_QTREE_.KSP[KS - 1] = 2;
		}
		else {
			if (B < 0.0e0)
				dg_QTREE_.KSP[KS - 1] = 1; //particula movendo-se para dentro
			else
				dg_QTREE_.KSP[KS - 1] = 2; // particula movendo-se para fora
		}
	d101:;
	}
	for (int KBB = 1; KBB <= dg_QTREE_.KDGHT[NXG - 1][KB0 - 1]; KBB++) {
		KB = dg_QTREE_.KDGHT[KBB - 1][KB0 - 1];
		for (int KSS = 1; KSS <= dg_QTREE_.KSURF[NXG - 1][KB - 1]; KSS++) {
			KS = dg_QTREE_.KSURF[KSS - 1][KB - 1];
			KF = dg_QTREE_.KFLAG[KSS - 1][KB - 1];

			if ((KF < 3) && (dg_QTREE_.KSP[KS - 1] != KF)) {
				goto d102;
			}
		}
		if (KB == KB0) {
			dg_TRACK_mod_.IBODY = KB; // a particula está dentro do corpo ou modulo KB
			dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[KB - 1];// a particula está dentro do MATERial KB    
			return;
		}
		else if (dg_QTREE_.KDGHT[NXG - 1][KB - 1] > 1) {
			KB0 = KB; //o ponto está dentro de um submodulo
			goto d100;
		}
		else {
			dg_TRACK_mod_.IBODY = KB; //a particula esta dentro de um corpo ou modulo irmão
			dg_TRACK_mod_.MAT = dg_PENGEOM_mod_.MATER[KB - 1];
			return;
		}
	d102:;
	}
	dg_TRACK_mod_.IBODY = dg_QTREE_.NBODYS + 1;
	dg_TRACK_mod_.MAT = 0;

}


