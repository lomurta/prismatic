extern "C"{

    void transfCPU_to_GPU();

    void transfGPU_to_CPU();

    void memoryFreeGPU();

    __global__ void teste();


}

void transfCPU_to_GPU(){

    /* é possivel transferir struct inteira desde que seja de tamanho estatico
     gpuErrchk(cudaMalloc((void **)&d2_CEELDB, sizeof(hd_CEELDB)));
     gpuErrchk(cudaMemcpy(d2_CEELDB, &h_CEELDB, sizeof(hd_CEELDB), cudaMemcpyHostToDevice));
     gpuErrchk(cudaMemcpyToSymbol(dg2_CEELDB, d2_CEELDB, sizeof(hd_CEELDB)));
     */

    //CEELDB
    gpuErrchk(cudaMalloc((void **)&d_CEELDB, sizeof(hd_CEELDB)));
    gpuErrchk(cudaMemcpy(d_CEELDB->XSE, CEELDB_.XSE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEELDB->PSE, CEELDB_.PSE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEELDB->ASE, CEELDB_.ASE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEELDB->BSE, CEELDB_.BSE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEELDB->ITLE, CEELDB_.ITLE, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEELDB->ITUE, CEELDB_.ITUE, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_CEELDB_, d_CEELDB, sizeof(hd_CEELDB)));

    //SECST
    gpuErrchk(cudaMalloc((void **)&d_SECST, sizeof(hd_SECST)));
    gpuErrchk(cudaMemcpy(d_SECST->ES, SECST_.ES, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->XS, SECST_.XS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->YS, SECST_.YS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->ZS, SECST_.ZS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->US, SECST_.US, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->VS, SECST_.VS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->WS, SECST_.WS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->WGHTS, SECST_.WGHTS, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->SP1S, SECST_.SP1S, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->SP2S, SECST_.SP2S, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->SP3S, SECST_.SP3S, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->PAGES, SECST_.PAGES, sizeof(double)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->KS, SECST_.KS, sizeof(int)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->IBODYS, SECST_.IBODYS, sizeof(int)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->MS, SECST_.MS, sizeof(int)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->ILBS, SECST_.ILBS, sizeof(int)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_SECST->IPOLS, SECST_.IPOLS, sizeof(int)*NMS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_SECST->NSEC, SECST_.NSEC, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_SECST_, d_SECST, sizeof(hd_SECST)));

    //PENGEOM_MOD
    gpuErrchk(cudaMalloc((void **)&d_PENGEOM_MOD, sizeof(hd_PENGEOM_MOD)));
    gpuErrchk(cudaMemcpy(d_PENGEOM_MOD->BALIAS, PENGEOM_mod_.BALIAS, sizeof(char)*NB*5, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_PENGEOM_MOD->DSTOT, PENGEOM_mod_.DSTOT, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_PENGEOM_MOD->MATER, PENGEOM_mod_.MATER, sizeof(int)*NB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_PENGEOM_MOD->KDET, PENGEOM_mod_.KDET, sizeof(int)*NB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_PENGEOM_MOD->KSLAST, PENGEOM_mod_.KSLAST, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_PENGEOM_MOD->NBODY, PENGEOM_mod_.NBODY, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_PENGEOM_MOD->LVERB, PENGEOM_mod_.LVERB, sizeof(bool), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_PENGEOM_MOD_, d_PENGEOM_MOD, sizeof(hd_PENGEOM_MOD)));

    //GCONE
    gpuErrchk(cudaMalloc((void **)&d_CGCONE, sizeof(hd_CGCONE)));
    gpuErrchk(cudaMemcpy(&d_CGCONE->CPCT, CGCONE_.CPCT, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->CPST, CGCONE_.CPST, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->SPCT, CGCONE_.SPCT, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->SPST, CGCONE_.SPST, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->SPHI, CGCONE_.SPHI, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->CPHI, CGCONE_.CPHI, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->STHE, CGCONE_.STHE, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CGCONE->CAPER, CGCONE_.CAPER, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_CGCONE_, d_CGCONE, sizeof(hd_CGCONE)));

    //QSURF
    gpuErrchk(cudaMalloc((void **)&d_QSURF, sizeof(hd_QSURF)));
    gpuErrchk(cudaMemcpy(d_QSURF->AXX, QSURF_.AXX, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AXY, QSURF_.AXY, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AXZ, QSURF_.AXZ, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AYY, QSURF_.AYY, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AYZ, QSURF_.AYZ, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AZZ, QSURF_.AZZ, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AX, QSURF_.AX, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AY, QSURF_.AY, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->AZ, QSURF_.AZ, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->A0, QSURF_.A0, sizeof(double)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_QSURF->NSURF, QSURF_.NSURF, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QSURF->KPLANE, QSURF_.KPLANE, sizeof(int)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_QSURF_, d_QSURF, sizeof(hd_QSURF)));

    //QTREE
    gpuErrchk(cudaMalloc((void **)&d_QTREE, sizeof(hd_QTREE)));
    gpuErrchk(cudaMemcpy(&d_QTREE->NBODYS,QTREE_.NBODYS, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QTREE->KMOTH, QTREE_.KMOTH, sizeof(int)*NB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QTREE->KDGHT, QTREE_.KDGHT, sizeof(int)*NB*NXG, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QTREE->KSURF, QTREE_.KSURF, sizeof(int)*NB*NXG, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QTREE->KFLAG, QTREE_.KFLAG, sizeof(int)*NB*NXG, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_QTREE->KSP, QTREE_.KSP, sizeof(int)*NS, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_QTREE->NWARN, QTREE_.NWARN, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_QTREE_, d_QTREE, sizeof(hd_QTREE)));

    //RSEED
    gpuErrchk(cudaMalloc((void **)&d_RSEED, sizeof(hd_RSEED)));
    gpuErrchk(cudaMemcpy(&d_RSEED->ISEED1,RSEED_.ISEED1, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_RSEED->ISEED2,RSEED_.ISEED2, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_RSEED_, d_RSEED, sizeof(hd_RSEED)));

    //CIMDET
    gpuErrchk(cudaMalloc((void **)&d_CIMDET, sizeof(hd_CIMDET)));
    gpuErrchk(cudaMemcpy(d_CIMDET->EL, CIMDET_.EL, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->EU, CIMDET_.EU, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->BSE, CIMDET_.BSE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->RBSE, CIMDET_.RBSE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->ET, CIMDET_.ET, sizeof(double)*NIDM*(NBEM2+1), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->EDEP, CIMDET_.EDEP, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->EDEP2, CIMDET_.EDEP2, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->EDEPP, CIMDET_.EDEPP, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DIT, CIMDET_.DIT, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DIT2, CIMDET_.DIT2, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DITP, CIMDET_.DITP, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DIP, CIMDET_.DIP, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DIP2, CIMDET_.DIP2, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->DIPP, CIMDET_.DIPP, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLT, CIMDET_.FLT, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLT2, CIMDET_.FLT2, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLTP, CIMDET_.FLTP, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLP, CIMDET_.FLP, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLP2, CIMDET_.FLP2, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->FLPP, CIMDET_.FLPP, sizeof(double)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->AGEL, CIMDET_.AGEL, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->AGEU, CIMDET_.AGEU, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->BAGE, CIMDET_.BAGE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->RBAGE, CIMDET_.RBAGE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->AGE, CIMDET_.AGE, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->AGE2, CIMDET_.AGE2, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->AGEP, CIMDET_.AGEP, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LEDEP, CIMDET_.LEDEP, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LDIT, CIMDET_.LDIT, sizeof(int)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LDIP, CIMDET_.LDIP, sizeof(int)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LFLT, CIMDET_.LFLT, sizeof(int)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LFLP, CIMDET_.LFLP, sizeof(int)*NIDM*NBEM2*3, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LAGEA, CIMDET_.LAGEA, sizeof(int)*NIDM*NBEM2, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->IDCUT, CIMDET_.IDCUT, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->NE, CIMDET_.NE, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LLE, CIMDET_.LLE, sizeof(bool)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->LLAGE, CIMDET_.LLAGE, sizeof(bool)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CIMDET->NAGE, CIMDET_.NAGE, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CIMDET->NID, CIMDET_.NID, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_CIMDET_, d_CIMDET, sizeof(hd_CIMDET)));

    //CERSEC
    gpuErrchk(cudaMalloc((void **)&d_CERSEC, sizeof(hd_CERSEC)));
    gpuErrchk(cudaMemcpy(&d_CERSEC->IERSEC,CERSEC_.IERSEC, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_CERSEC_, d_CERSEC, sizeof(hd_CERSEC)));

    //CEGRID
    gpuErrchk(cudaMalloc((void **)&d_CEGRID, sizeof(hd_CEGRID)));
    gpuErrchk(cudaMemcpy(&d_CEGRID->EMIN,CEGRID_.EMIN, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->EL,CEGRID_.EL, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->EU,CEGRID_.EU, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEGRID->ET,CEGRID_.ET, sizeof(double)*NEGP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_CEGRID->DLEMP,CEGRID_.DLEMP, sizeof(double)*NEGP, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->DLEMP1,CEGRID_.DLEMP1, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->DLFC,CEGRID_.DLFC, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->XEL,CEGRID_.XEL, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->XE,CEGRID_.XE, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->XEK,CEGRID_.XEK, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_CEGRID->KE,CEGRID_.KE, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpyToSymbol(dg_CEGRID_, d_CEGRID, sizeof(hd_CEGRID)));

	//CJUMP1
    gpuErrchk(cudaMalloc((void **)&d_CJUMP1, sizeof(hd_CJUMP1)));
    gpuErrchk(cudaMemcpy(&d_CJUMP1->ELAST1,CJUMP1_.ELAST1, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP1->ELAST2,CJUMP1_.ELAST2, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP1->MHINGE,CJUMP1_.MHINGE, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP1->KSOFTE,CJUMP1_.KSOFTE, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP1->KSOFTI,CJUMP1_.KSOFTI, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP1->KDELTA,CJUMP1_.KDELTA, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CJUMP1_, d_CJUMP1, sizeof(hd_CJUMP1)));

	//CEIMFP
	gpuErrchk(cudaMalloc((void **)&d_CEIMFP, sizeof(hd_CEIMFP)));
    gpuErrchk(cudaMemcpy(d_CEIMFP->SEHEL,CEIMFP_.SEHEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->SEHIN,CEIMFP_.SEHIN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->SEISI,CEIMFP_.SEISI, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->SEHBR,CEIMFP_.SEHBR, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->SEAUX,CEIMFP_.SEAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->SETOT,CEIMFP_.SETOT, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->CSTPE,CEIMFP_.CSTPE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->RSTPE,CEIMFP_.RSTPE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->DEL,CEIMFP_.DEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->W1E,CEIMFP_.W1E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->W2E,CEIMFP_.W2E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->DW1EL,CEIMFP_.DW1EL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->DW2EL,CEIMFP_.DW2EL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->RNDCE,CEIMFP_.RNDCE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->AE,CEIMFP_.AE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->BE,CEIMFP_.BE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->T1E,CEIMFP_.T1E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIMFP->T2E,CEIMFP_.T2E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CEIMFP_, d_CEIMFP, sizeof(hd_CEIMFP)));

	//CPIMFP
	gpuErrchk(cudaMalloc((void **)&d_CPIMFP, sizeof(hd_CPIMFP)));
    gpuErrchk(cudaMemcpy(d_CPIMFP->SPHEL,CPIMFP_.SPHEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPHIN,CPIMFP_.SPHIN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPISI,CPIMFP_.SPISI, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPHBR,CPIMFP_.SPHBR, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPAN,CPIMFP_.SPAN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPAUX,CPIMFP_.SPAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->SPTOT,CPIMFP_.SPTOT, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->CSTPP,CPIMFP_.CSTPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->RSTPP,CPIMFP_.RSTPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->W1P,CPIMFP_.W1P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->W2P,CPIMFP_.W2P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->DW1PL,CPIMFP_.DW1PL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->DW2PL,CPIMFP_.DW2PL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->RNDCP,CPIMFP_.RNDCP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->AP,CPIMFP_.AP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->BP,CPIMFP_.BP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->T1P,CPIMFP_.T1P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPIMFP->T2P,CPIMFP_.T2P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CPIMFP_, d_CPIMFP, sizeof(hd_CPIMFP)));

	//CJUMP0
	gpuErrchk(cudaMalloc((void **)&d_CJUMP0, sizeof(hd_CJUMP0)));
    gpuErrchk(cudaMemcpy(d_CJUMP0->P,CJUMP0_.P, sizeof(double)*8, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->ST,CJUMP0_.ST, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->DST,CJUMP0_.DST, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->DSR,CJUMP0_.DSR, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->W1,CJUMP0_.W1, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->W2,CJUMP0_.W2, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->T1,CJUMP0_.T1, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CJUMP0->T2,CJUMP0_.T2, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CJUMP0_, d_CJUMP0, sizeof(hd_CJUMP0)));

	//CGIMFP
	gpuErrchk(cudaMalloc((void **)&d_CGIMFP, sizeof(hd_CGIMFP)));
    gpuErrchk(cudaMemcpy(d_CGIMFP->SGRA,CGIMFP_.SGRA, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGIMFP->SGCO,CGIMFP_.SGCO, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGIMFP->SGPH,CGIMFP_.SGPH, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGIMFP->SGPP,CGIMFP_.SGPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGIMFP->SGAUX,CGIMFP_.SGAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGIMFP_, d_CGIMFP, sizeof(hd_CGIMFP)));

	//CHIST
	gpuErrchk(cudaMalloc((void **)&d_CHIST, sizeof(hd_CHIST)));
	gpuErrchk(cudaMemcpy(d_CHIST->ILBA,CHIST_.ILBA, sizeof(int)*5, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CHIST_, d_CHIST, sizeof(hd_CHIST)));

	//COMPOS
	gpuErrchk(cudaMalloc((void **)&d_COMPOS, sizeof(hd_COMPOS)));
	gpuErrchk(cudaMemcpy(d_COMPOS->STF,COMPOS_.STF, sizeof(double)*30*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->ZT,COMPOS_.ZT, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->AT,COMPOS_.AT, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->RHO,COMPOS_.RHO, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->VMOL,COMPOS_.VMOL, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->IZ,COMPOS_.IZ, sizeof(int)*30*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_COMPOS->NELEM,COMPOS_.NELEM, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_COMPOS_, d_COMPOS, sizeof(hd_COMPOS)));

	//CADATA
	gpuErrchk(cudaMalloc((void **)&d_CADATA, sizeof(hd_CADATA)));
	gpuErrchk(cudaMemcpy(d_CADATA->ATW,CADATA_.ATW, sizeof(double)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->EPX,CADATA_.EPX, sizeof(double)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->RSCR,CADATA_.RSCR, sizeof(double)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->ETA,CADATA_.ETA, sizeof(double)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->EB,CADATA_.EB, sizeof(double)*99*30, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->ALW,CADATA_.ALW, sizeof(double)*99*30, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->CP0,CADATA_.CP0, sizeof(double)*99*30, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->IFI,CADATA_.IFI, sizeof(int)*99*30, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->IKS,CADATA_.IKS, sizeof(int)*99*30, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->NSHT,CADATA_.NSHT, sizeof(int)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CADATA->LASYMB,CADATA_.LASYMB, sizeof(char)*99*3, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CADATA_, d_CADATA, sizeof(hd_CADATA)));

	//CEIN
	gpuErrchk(cudaMalloc((void **)&d_CEIN, sizeof(hd_CEIN)));
	gpuErrchk(cudaMemcpy(d_CEIN->EXPOT,CEIN_.EXPOT, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->OP2,CEIN_.OP2, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->F,CEIN_.F, sizeof(double)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->UI,CEIN_.UI, sizeof(double)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->WRI,CEIN_.WRI, sizeof(double)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->KZ,CEIN_.KZ, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->KS,CEIN_.KS, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEIN->NOSC,CEIN_.NOSC, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CEIN_, d_CEIN, sizeof(hd_CEIN)));

	//CGCO
	gpuErrchk(cudaMalloc((void **)&d_CGCO, sizeof(hd_CGCO)));
	gpuErrchk(cudaMemcpy(d_CGCO->FCO,CGCO_.FCO, sizeof(double)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->UICO,CGCO_.UICO, sizeof(double)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->FJ0,CGCO_.FJ0, sizeof(double)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->PTRSH,CGCO_.PTRSH, sizeof(double)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->KZCO,CGCO_.KZCO, sizeof(int)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->KSCO,CGCO_.KSCO, sizeof(int)*MAXMAT*NOCO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGCO->NOSCCO,CGCO_.NOSCCO, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGCO_, d_CGCO, sizeof(hd_CGCO)));

	//CEBR
	gpuErrchk(cudaMalloc((void **)&d_CEBR, sizeof(hd_CEBR)));
	gpuErrchk(cudaMemcpy(d_CEBR->WB,CEBR_.WB, sizeof(double)*NBW, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->PBCUT,CEBR_.PBCUT, sizeof(double)*MAXMAT*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->WBCUT,CEBR_.WBCUT, sizeof(double)*MAXMAT*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->PDFB,CEBR_.PDFB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->DPDFB,CEBR_.DPDFB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->PACB,CEBR_.PACB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEBR->ZBR2,CEBR_.ZBR2, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CEBR_, d_CEBR, sizeof(hd_CEBR)));

	//CELSEP
	gpuErrchk(cudaMalloc((void **)&d_CELSEP, sizeof(hd_CELSEP)));
	gpuErrchk(cudaMemcpy(d_CELSEP->EELMAX,CELSEP_.EELMAX, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CELSEP->PELMAX,CELSEP_.PELMAX, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CELSEP->RNDCED,CELSEP_.RNDCED, sizeof(double)*MAXMAT*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CELSEP->RNDCPD,CELSEP_.RNDCPD, sizeof(double)*MAXMAT*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CELSEP_, d_CELSEP, sizeof(hd_CELSEP)));

	//CECUTR
	gpuErrchk(cudaMalloc((void **)&d_CECUTR, sizeof(hd_CECUTR)));
	gpuErrchk(cudaMemcpy(d_CECUTR->ECUTR,CECUTR_.ECUTR, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CECUTR_, d_CECUTR, sizeof(hd_CECUTR)));

	//CESIAC
	gpuErrchk(cudaMalloc((void **)&d_CESIAC, sizeof(hd_CESIAC)));
	gpuErrchk(cudaMemcpy(d_CESIAC->ESIAC,CESIAC_.ESIAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CESIAC->IESI,CESIAC_.IESI, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CESIAC->NESI,CESIAC_.NESI, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CESIAC_, d_CESIAC, sizeof(hd_CESIAC)));

	//CEINAC
	gpuErrchk(cudaMalloc((void **)&d_CEINAC, sizeof(hd_CEINAC)));
	gpuErrchk(cudaMemcpy(d_CEINAC->EINAC,CEINAC_.EINAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEINAC->IEIN,CEINAC_.IEIN, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CEINAC->NEIN,CEINAC_.NEIN, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CEINAC_, d_CEINAC, sizeof(hd_CEINAC)));

	//CBRANG
	gpuErrchk(cudaMalloc((void **)&d_CBRANG, sizeof(hd_CBRANG)));
	gpuErrchk(cudaMemcpy(d_CBRANG->BET,CBRANG_.BET, sizeof(double)*6, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CBRANG->BK,CBRANG_.BK, sizeof(double)*21, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CBRANG->BP1,CBRANG_.BP1, sizeof(double)*4*21*6*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CBRANG->BP2,CBRANG_.BP2, sizeof(double)*4*21*6*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CBRANG->ZBEQ,CBRANG_.ZBEQ, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CBRANG_, d_CBRANG, sizeof(hd_CBRANG)));

	//CRELAX
	gpuErrchk(cudaMalloc((void **)&d_CRELAX, sizeof(hd_CRELAX)));
	gpuErrchk(cudaMemcpy(d_CRELAX->P,CRELAX_.P, sizeof(double)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->ET,CRELAX_.ET, sizeof(double)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->F,CRELAX_.F, sizeof(double)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->IS0,CRELAX_.IS0, sizeof(int)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->IS1,CRELAX_.IS1, sizeof(int)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->IS2,CRELAX_.IS2, sizeof(int)*NRX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->IFIRST,CRELAX_.IFIRST, sizeof(int)*16*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CRELAX->ILAST,CRELAX_.ILAST, sizeof(int)*16*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CRELAX->NCUR,CRELAX_.NCUR, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CRELAX->KS,CRELAX_.KS, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CRELAX->MODER,CRELAX_.MODER, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CRELAX_, d_CRELAX, sizeof(hd_CRELAX)));

	//CGRA01
	gpuErrchk(cudaMalloc((void **)&d_CGRA01, sizeof(hd_CGRA01)));
	gpuErrchk(cudaMemcpy(d_CGRA01->FF,CGRA01_.FF, sizeof(double)*NQ*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA01->ERA,CGRA01_.ERA, sizeof(double)*NEX, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA01->XSRA,CGRA01_.XSRA, sizeof(double)*NEX*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA01->IED,CGRA01_.IED, sizeof(int)*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA01->IEU,CGRA01_.IEU, sizeof(int)*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CGRA01->NE,CGRA01_.NE, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGRA01_, d_CGRA01, sizeof(hd_CGRA01)));

	//CGRA03
	gpuErrchk(cudaMalloc((void **)&d_CGRA03, sizeof(hd_CGRA03)));
	gpuErrchk(cudaMemcpy(d_CGRA03->QRA,CGRA03_.QRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->PRA,CGRA03_.PRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->DPRA,CGRA03_.DPRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->ARA,CGRA03_.ARA, sizeof(double)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->BRA,CGRA03_.BRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->PMAX,CGRA03_.PMAX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->ITLRA,CGRA03_.ITLRA, sizeof(int)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGRA03->ITURA,CGRA03_.ITURA, sizeof(int)*NP2*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGRA03_, d_CGRA03, sizeof(hd_CGRA03)));

	//CGPH00
	gpuErrchk(cudaMalloc((void **)&d_CGPH00, sizeof(hd_CGPH00)));
	gpuErrchk(cudaMemcpy(d_CGPH00->EPH,CGPH00_.EPH, sizeof(double)*NTP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPH00->XPH,CGPH00_.XPH, sizeof(double)*NTP*17, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPH00->IPHF,CGPH00_.IPHF, sizeof(int)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPH00->IPHL,CGPH00_.IPHL, sizeof(int)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPH00->NPHS,CGPH00_.NPHS, sizeof(int)*99, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CGPH00->NCUR,CGPH00_.NCUR, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGPH00_, d_CGPH00, sizeof(hd_CGPH00)));

	//CGPP00
	gpuErrchk(cudaMalloc((void **)&d_CGPP00, sizeof(hd_CGPP00)));
	gpuErrchk(cudaMemcpy(d_CGPP00->ZEQPP,CGPP00_.ZEQPP, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPP00->F0,CGPP00_.F0, sizeof(double)*MAXMAT*2, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CGPP00->BCB,CGPP00_.BCB, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGPP00_, d_CGPP00, sizeof(hd_CGPP00)));

	//CGPP01
	gpuErrchk(cudaMalloc((void **)&d_CGPP01, sizeof(hd_CGPP01)));
	gpuErrchk(cudaMemcpy(d_CGPP01->TRIP,CGPP01_.TRIP, sizeof(double)*MAXMAT*NEGP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CGPP01_, d_CGPP01, sizeof(hd_CGPP01)));

	//CPELDB
	gpuErrchk(cudaMalloc((void **)&d_CPELDB, sizeof(hd_CPELDB)));
	gpuErrchk(cudaMemcpy(d_CPELDB->XSP,CPELDB_.XSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPELDB->PSP,CPELDB_.PSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPELDB->ASP,CPELDB_.ASP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPELDB->BSP,CPELDB_.BSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPELDB->ITLP,CPELDB_.ITLP, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPELDB->ITUP,CPELDB_.ITUP, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CPELDB_, d_CPELDB, sizeof(hd_CPELDB)));

	//CPINAC
	gpuErrchk(cudaMalloc((void **)&d_CPINAC, sizeof(hd_CPINAC)));
	gpuErrchk(cudaMemcpy(d_CPINAC->PINAC,CPINAC_.PINAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPINAC->IPIN,CPINAC_.IPIN, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPINAC->NPIN,CPINAC_.NPIN, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CPINAC_, d_CPINAC, sizeof(hd_CPINAC)));

	//CPSIAC
	gpuErrchk(cudaMalloc((void **)&d_CPSIAC, sizeof(hd_CPSIAC)));
	gpuErrchk(cudaMemcpy(d_CPSIAC->PSIAC,CPSIAC_.PSIAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPSIAC->IPSI,CPSIAC_.IPSI, sizeof(int)*MAXMAT*NO, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CPSIAC->NPSI,CPSIAC_.NPSI, sizeof(int)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CPSIAC_, d_CPSIAC, sizeof(hd_CPSIAC)));


	//CENANG
	gpuErrchk(cudaMalloc((void **)&d_CENANG, sizeof(hd_CENANG)));
	gpuErrchk(cudaMemcpy(&d_CENANG->EL,CENANG_.EL, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->EU,CENANG_.EU, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->THL,CENANG_.THL, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->THU,CENANG_.THU, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->BSE,CENANG_.BSE, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->RBSE,CENANG_.RBSE, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->BSTH,CENANG_.BSTH, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->RBSTH,CENANG_.RBSTH, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->BSPH,CENANG_.BSPH, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->RBSPH,CENANG_.RBSPH, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDE,CENANG_.PDE, sizeof(double)*2*3*NBEM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDE2,CENANG_.PDE2, sizeof(double)*2*3*NBEM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDEP,CENANG_.PDEP, sizeof(double)*2*3*NBEM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDA,CENANG_.PDA, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDA2,CENANG_.PDA2, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->PDAP,CENANG_.PDAP, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->LPDE,CENANG_.LPDE, sizeof(int)*3*2*NBEM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->LPDA,CENANG_.LPDA, sizeof(int)*3*NBPHM*NBTHM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->NE,CENANG_.NE, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->NTH,CENANG_.NTH, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->NPH,CENANG_.NPH, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->LLE,CENANG_.LLE, sizeof(bool), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENANG->LLTH,CENANG_.LLTH, sizeof(bool), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CENANG_, d_CENANG, sizeof(hd_CENANG)));

	//CENDET
	gpuErrchk(cudaMalloc((void **)&d_CENDET, sizeof(hd_CENDET)));
	gpuErrchk(cudaMemcpy(d_CENDET->EL,CENDET_.EL, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->EU,CENDET_.EU, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->BSE,CENDET_.BSE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->RBSE,CENDET_.RBSE, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->EDEP,CENDET_.EDEP, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->EDEP2,CENDET_.EDEP2, sizeof(double)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->DET,CENDET_.DET, sizeof(double)*NIDM*NBEM2, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->NE,CENDET_.NE, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_CENDET->NID,CENDET_.NID, sizeof(int)*NIDM, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_CENDET->LLE,CENDET_.LLE, sizeof(bool), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_CENDET_, d_CENDET, sizeof(hd_CENDET)));

	//PENELOPE_MOD
	gpuErrchk(cudaMalloc((void **)&d_PENELOPE_mod, sizeof(hd_PENELOPE_MOD)));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->EABS,PENELOPE_mod_.EABS, sizeof(double)*3*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->C1,PENELOPE_mod_.C1, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->C2,PENELOPE_mod_.C2, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->WCC,PENELOPE_mod_.WCC, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->WCR,PENELOPE_mod_.WCR, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->DEN,PENELOPE_mod_.DEN, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(d_PENELOPE_mod->RDEN,PENELOPE_mod_.RDEN, sizeof(double)*MAXMAT, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->E0STEP,PENELOPE_mod_.E0STEP, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->DESOFT, PENELOPE_mod_.DESOFT, sizeof(double), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->SSOFT, PENELOPE_mod_.SSOFT, sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->NMS, PENELOPE_mod_.NMS, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->NEGP, PENELOPE_mod_.NEGP, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(&d_PENELOPE_mod->NMAT, PENELOPE_mod_.NMAT, sizeof(int), cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(dg_PENELOPE_mod_, d_PENELOPE_mod, sizeof(hd_PENELOPE_MOD)));





}


void transfGPU_to_CPU(){

    //CEELDB
    gpuErrchk(cudaMemcpyFromSymbol(d_CEELDB, dg_CEELDB_, sizeof(hd_CEELDB)));
    gpuErrchk(cudaMemcpy(CEELDB_.XSE, d_CEELDB->XSE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEELDB_.PSE, d_CEELDB->PSE, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEELDB_.ASE, d_CEELDB->ASE,  sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEELDB_.BSE, d_CEELDB->BSE,  sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEELDB_.ITLE, d_CEELDB->ITLE,  sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEELDB_.ITUE, d_CEELDB->ITUE, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));

    //SECST
    gpuErrchk(cudaMemcpyFromSymbol(d_SECST, dg_SECST_, sizeof(hd_SECST)));
    gpuErrchk(cudaMemcpy(SECST_.ES,d_SECST->ES,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.XS,d_SECST->XS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.YS,d_SECST->YS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.ZS,d_SECST->ZS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.US,d_SECST->US,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.VS,d_SECST->VS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.WS,d_SECST->WS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.WGHTS,d_SECST->WGHTS,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.SP1S,d_SECST->SP1S,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.SP2S,d_SECST->SP2S,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.SP3S,d_SECST->SP3S,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.PAGES,d_SECST->PAGES,  sizeof(double)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.KS,d_SECST->KS,  sizeof(int)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.IBODYS,d_SECST->IBODYS,  sizeof(int)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.MS,d_SECST->MS,  sizeof(int)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.ILBS,d_SECST->ILBS, sizeof(int)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.IPOLS, d_SECST->IPOLS, sizeof(int)*NMS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(SECST_.NSEC, &d_SECST->NSEC,  sizeof(int), cudaMemcpyDeviceToHost));

    //PENGEOM_MOD
    gpuErrchk(cudaMemcpyFromSymbol(d_PENGEOM_MOD, dg_PENGEOM_MOD_, sizeof(hd_PENGEOM_MOD)));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.BALIAS,d_PENGEOM_MOD->BALIAS,  sizeof(char)*NB*5, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.DSTOT, &d_PENGEOM_MOD->DSTOT,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.MATER,d_PENGEOM_MOD->MATER,  sizeof(int)*NB, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.KDET,d_PENGEOM_MOD->KDET,  sizeof(int)*NB, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.KSLAST,&d_PENGEOM_MOD->KSLAST,  sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.NBODY,&d_PENGEOM_MOD->NBODY,  sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENGEOM_mod_.LVERB, &d_PENGEOM_MOD->LVERB,  sizeof(bool), cudaMemcpyDeviceToHost));

    //GCONE
    gpuErrchk(cudaMemcpyFromSymbol(d_CGCONE, dg_CGCONE_, sizeof(hd_CGCONE)));
    gpuErrchk(cudaMemcpy(CGCONE_.CPCT,&d_CGCONE->CPCT,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CGCONE_.CPST,&d_CGCONE->CPST,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CGCONE_.SPCT,&d_CGCONE->SPCT,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CGCONE_.SPST,&d_CGCONE->SPST, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CGCONE_.SPHI,&d_CGCONE->SPHI,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( CGCONE_.CPHI,&d_CGCONE->CPHI, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( CGCONE_.STHE,&d_CGCONE->STHE, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CGCONE_.CAPER,&d_CGCONE->CAPER,  sizeof(double), cudaMemcpyDeviceToHost));

    //QSURF
    gpuErrchk(cudaMemcpyFromSymbol(d_QSURF, dg_QSURF_, sizeof(hd_QSURF)));
    gpuErrchk(cudaMemcpy(QSURF_.AXX,d_QSURF->AXX,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AXY,d_QSURF->AXY,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AXZ,d_QSURF->AXZ,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AYY, d_QSURF->AYY, sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AYZ,d_QSURF->AYZ,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AZZ,d_QSURF->AZZ,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AX,d_QSURF->AX,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AY,d_QSURF->AY,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.AZ,d_QSURF->AZ,  sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( QSURF_.A0,d_QSURF->A0, sizeof(double)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QSURF_.NSURF,&d_QSURF->NSURF,  sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( QSURF_.KPLANE,d_QSURF->KPLANE, sizeof(int)*NS, cudaMemcpyDeviceToHost));

    //QTREE
    gpuErrchk(cudaMemcpyFromSymbol(d_QTREE, dg_QTREE_, sizeof(hd_QTREE)));
    gpuErrchk(cudaMemcpy(QTREE_.NBODYS, &d_QTREE->NBODYS,sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.KMOTH,d_QTREE->KMOTH,  sizeof(int)*NB, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.KDGHT,d_QTREE->KDGHT,  sizeof(int)*NB*NXG, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.KSURF,d_QTREE->KSURF,  sizeof(int)*NB*NXG, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.KFLAG,d_QTREE->KFLAG,  sizeof(int)*NB*NXG, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.KSP,d_QTREE->KSP,  sizeof(int)*NS, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(QTREE_.NWARN, &d_QTREE->NWARN,  sizeof(int), cudaMemcpyDeviceToHost));

    //RSEED
    gpuErrchk(cudaMemcpyFromSymbol(d_RSEED, dg_RSEED_, sizeof(hd_RSEED)));
    gpuErrchk(cudaMemcpy(RSEED_.ISEED1,&d_RSEED->ISEED1, sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(RSEED_.ISEED2,&d_RSEED->ISEED2, sizeof(int), cudaMemcpyDeviceToHost));


     //CIMDET
    gpuErrchk(cudaMemcpyFromSymbol(d_CIMDET, dg_CIMDET_, sizeof(hd_CIMDET)));
    gpuErrchk(cudaMemcpy( CIMDET_.EL,d_CIMDET->EL, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.EU,d_CIMDET->EU,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.BSE,d_CIMDET->BSE,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.RBSE,d_CIMDET->RBSE,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.ET,d_CIMDET->ET,  sizeof(double)*NIDM*(NBEM2+1), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.EDEP,d_CIMDET->EDEP,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.EDEP2,d_CIMDET->EDEP2,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.EDEPP,d_CIMDET->EDEPP,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DIT,d_CIMDET->DIT,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DIT2,d_CIMDET->DIT2,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DITP,d_CIMDET->DITP,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DIP,d_CIMDET->DIP,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DIP2,d_CIMDET->DIP2,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.DIPP,d_CIMDET->DIPP,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLT,d_CIMDET->FLT,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLT2,d_CIMDET->FLT2,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLTP,d_CIMDET->FLTP,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLP,d_CIMDET->FLP,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLP2,d_CIMDET->FLP2,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.FLPP,d_CIMDET->FLPP,  sizeof(double)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.AGEL,d_CIMDET->AGEL,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.AGEU,d_CIMDET->AGEU,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.BAGE,d_CIMDET->BAGE,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.RBAGE,d_CIMDET->RBAGE,  sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.AGE,d_CIMDET->AGE,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.AGE2,d_CIMDET->AGE2,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.AGEP,d_CIMDET->AGEP,  sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LEDEP,d_CIMDET->LEDEP,  sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LDIT,d_CIMDET->LDIT,  sizeof(int)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LDIP,d_CIMDET->LDIP,  sizeof(int)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LFLT,d_CIMDET->LFLT,  sizeof(int)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LFLP,d_CIMDET->LFLP,  sizeof(int)*NIDM*NBEM2*3, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.LAGEA,d_CIMDET->LAGEA,  sizeof(int)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.IDCUT,d_CIMDET->IDCUT,  sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.NE,d_CIMDET->NE,  sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( CIMDET_.LLE,d_CIMDET->LLE, sizeof(bool)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy( CIMDET_.LLAGE, d_CIMDET->LLAGE,sizeof(bool)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.NAGE,d_CIMDET->NAGE,  sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CIMDET_.NID,&d_CIMDET->NID,  sizeof(int), cudaMemcpyDeviceToHost));

    //CERSEC
    gpuErrchk(cudaMemcpyFromSymbol(d_CERSEC, dg_CERSEC_, sizeof(hd_CERSEC)));
    gpuErrchk(cudaMemcpy(CERSEC_.IERSEC,&d_CERSEC->IERSEC, sizeof(int), cudaMemcpyDeviceToHost));

    //CEGRID
    gpuErrchk(cudaMemcpyFromSymbol(d_CEGRID, dg_CEGRID_, sizeof(hd_CEGRID)));
    gpuErrchk(cudaMemcpy(CEGRID_.EMIN,&d_CEGRID->EMIN, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.EL,&d_CEGRID->EL, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.EU,&d_CEGRID->EU, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.ET,d_CEGRID->ET, sizeof(double)*NEGP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.DLEMP,d_CEGRID->DLEMP, sizeof(double)*NEGP, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.DLEMP1,&d_CEGRID->DLEMP1, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.DLFC,&d_CEGRID->DLFC, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.XEL,&d_CEGRID->XEL, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.XE,&d_CEGRID->XE, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.XEK,&d_CEGRID->XEK, sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(CEGRID_.KE,&d_CEGRID->KE, sizeof(int), cudaMemcpyDeviceToHost));

	//CJUMP1
	gpuErrchk(cudaMemcpyFromSymbol(d_CJUMP1, dg_CJUMP1_, sizeof(hd_CJUMP1)));
	gpuErrchk(cudaMemcpy(CJUMP1_.ELAST1,&d_CJUMP1->ELAST1, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP1_.ELAST2,&d_CJUMP1->ELAST2, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP1_.MHINGE,&d_CJUMP1->MHINGE, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP1_.KSOFTE,&d_CJUMP1->KSOFTE, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP1_.KSOFTI,&d_CJUMP1->KSOFTI, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP1_.KDELTA,&d_CJUMP1->KDELTA, sizeof(int), cudaMemcpyDeviceToHost));

	//CEIMFP
	gpuErrchk(cudaMemcpyFromSymbol(d_CEIMFP, dg_CEIMFP_, sizeof(hd_CEIMFP)));
	gpuErrchk(cudaMemcpy(CEIMFP_.SEHEL,d_CEIMFP->SEHEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.SEHIN,d_CEIMFP->SEHIN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.SEISI,d_CEIMFP->SEISI, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.SEHBR,d_CEIMFP->SEHBR, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.SEAUX,d_CEIMFP->SEAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.SETOT,d_CEIMFP->SETOT, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.CSTPE,d_CEIMFP->CSTPE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.RSTPE,d_CEIMFP->RSTPE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.DEL,d_CEIMFP->DEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.W1E,d_CEIMFP->W1E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.W2E,d_CEIMFP->W2E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.DW1EL,d_CEIMFP->DW1EL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.DW2EL,d_CEIMFP->DW2EL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.RNDCE,d_CEIMFP->RNDCE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.AE,d_CEIMFP->AE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.BE,d_CEIMFP->BE, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.T1E,d_CEIMFP->T1E, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIMFP_.T2E, d_CEIMFP->T2E,sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));

	//CPIMFP
	gpuErrchk(cudaMemcpyFromSymbol(d_CEIMFP, dg_CPIMFP_, sizeof(hd_CPIMFP)));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPHEL,d_CPIMFP->SPHEL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPHIN,d_CPIMFP->SPHIN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPISI,d_CPIMFP->SPISI, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPHBR,d_CPIMFP->SPHBR, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPAN,d_CPIMFP->SPAN, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPAUX,d_CPIMFP->SPAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.SPTOT,d_CPIMFP->SPTOT, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.CSTPP,d_CPIMFP->CSTPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.RSTPP,d_CPIMFP->RSTPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.W1P,d_CPIMFP->W1P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.W2P,d_CPIMFP->W2P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.DW1PL,d_CPIMFP->DW1PL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.DW2PL,d_CPIMFP->DW2PL, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.RNDCP,d_CPIMFP->RNDCP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.AP,d_CPIMFP->AP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.BP,d_CPIMFP->BP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.T1P,d_CPIMFP->T1P, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPIMFP_.T2P, d_CPIMFP->T2P,sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));

	//CJUMP0
	gpuErrchk(cudaMemcpyFromSymbol(d_CJUMP0, dg_CJUMP0_, sizeof(hd_CJUMP0)));
	gpuErrchk(cudaMemcpy(CJUMP0_.P,d_CJUMP0->P, sizeof(double)*8, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.ST,&d_CJUMP0->ST, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.DST,&d_CJUMP0->DST, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.DSR,&d_CJUMP0->DSR, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.W1,&d_CJUMP0->W1, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.W2,&d_CJUMP0->W2, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.T1,&d_CJUMP0->T1, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CJUMP0_.T2,&d_CJUMP0->T2, sizeof(double), cudaMemcpyDeviceToHost));

	//CGIMFP
	gpuErrchk(cudaMemcpyFromSymbol(d_CGIMFP, dg_CGIMFP_, sizeof(hd_CGIMFP)));
	gpuErrchk(cudaMemcpy(CGIMFP_.SGRA,d_CGIMFP->SGRA, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGIMFP_.SGCO,d_CGIMFP->SGCO, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGIMFP_.SGPH,d_CGIMFP->SGPH, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGIMFP_.SGPP,d_CGIMFP->SGPP, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGIMFP_.SGAUX,d_CGIMFP->SGAUX, sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));

	//CHIST
	gpuErrchk(cudaMemcpyFromSymbol(d_CHIST, dg_CHIST_, sizeof(hd_CHIST)));
	gpuErrchk(cudaMemcpy(CHIST_.ILBA,d_CHIST->ILBA, sizeof(int)*5, cudaMemcpyDeviceToHost));

	//COMPOS
	gpuErrchk(cudaMemcpyFromSymbol(d_COMPOS, dg_COMPOS_, sizeof(hd_COMPOS)));
	gpuErrchk(cudaMemcpy(COMPOS_.STF,d_COMPOS->STF, sizeof(double)*30*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.ZT,d_COMPOS->ZT, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.AT,d_COMPOS->AT, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.RHO,d_COMPOS->RHO, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.VMOL,d_COMPOS->VMOL, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.IZ, d_COMPOS->IZ,sizeof(int)*30*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(COMPOS_.NELEM,d_COMPOS->NELEM, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CADATA
	gpuErrchk(cudaMemcpyFromSymbol(d_CADATA, dg_CADATA_, sizeof(hd_CADATA)));
	gpuErrchk(cudaMemcpy(CADATA_.ATW,d_CADATA->ATW, sizeof(double)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.EPX,d_CADATA->EPX, sizeof(double)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.RSCR,d_CADATA->RSCR, sizeof(double)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.ETA,d_CADATA->ETA, sizeof(double)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.EB,d_CADATA->EB, sizeof(double)*99*30, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.ALW,d_CADATA->ALW, sizeof(double)*99*30, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.CP0,d_CADATA->CP0, sizeof(double)*99*30, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.IFI,d_CADATA->IFI, sizeof(int)*99*30, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.IKS,d_CADATA->IKS, sizeof(int)*99*30, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.NSHT,d_CADATA->NSHT, sizeof(int)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CADATA_.LASYMB,d_CADATA->LASYMB,sizeof(char)*99*3, cudaMemcpyDeviceToHost));

	//CEIN
	gpuErrchk(cudaMemcpyFromSymbol(d_CEIN, dg_CEIN_, sizeof(hd_CEIN)));
	gpuErrchk(cudaMemcpy(CEIN_.EXPOT,d_CEIN->EXPOT, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.OP2,d_CEIN->OP2, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.F,d_CEIN->F, sizeof(double)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.UI,d_CEIN->UI, sizeof(double)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.WRI,d_CEIN->WRI, sizeof(double)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.KZ,d_CEIN->KZ, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.KS,d_CEIN->KS, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEIN_.NOSC,d_CEIN->NOSC, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CGCO
	gpuErrchk(cudaMemcpyFromSymbol(d_CGCO, dg_CGCO_, sizeof(hd_CGCO)));
	gpuErrchk(cudaMemcpy(CGCO_.FCO,d_CGCO->FCO, sizeof(double)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.UICO,d_CGCO->UICO, sizeof(double)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.FJ0,d_CGCO->FJ0, sizeof(double)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.PTRSH,d_CGCO->PTRSH, sizeof(double)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.KZCO,d_CGCO->KZCO, sizeof(int)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.KSCO,d_CGCO->KSCO, sizeof(int)*MAXMAT*NOCO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGCO_.NOSCCO,d_CGCO->NOSCCO, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CEBR
	gpuErrchk(cudaMemcpyFromSymbol(d_CEBR, dg_CEBR_, sizeof(hd_CEBR)));
	gpuErrchk(cudaMemcpy(CEBR_.WB,d_CEBR->WB, sizeof(double)*NBW, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.PBCUT,d_CEBR->PBCUT, sizeof(double)*MAXMAT*NEGP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.WBCUT,d_CEBR->WBCUT, sizeof(double)*MAXMAT*NEGP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.PDFB,d_CEBR->PDFB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.DPDFB,d_CEBR->DPDFB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.PACB,d_CEBR->PACB, sizeof(double)*NBW*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEBR_.ZBR2,d_CEBR->ZBR2, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));

	//CELSEP
	gpuErrchk(cudaMemcpyFromSymbol(d_CELSEP, dg_CELSEP_, sizeof(hd_CELSEP)));
	gpuErrchk(cudaMemcpy(CELSEP_.EELMAX,d_CELSEP->EELMAX, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CELSEP_.PELMAX,d_CELSEP->PELMAX, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CELSEP_.RNDCED,d_CELSEP->RNDCED, sizeof(double)*MAXMAT*NEGP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CELSEP_.RNDCPD,d_CELSEP->RNDCPD, sizeof(double)*MAXMAT*NEGP, cudaMemcpyDeviceToHost));

	//CECUTR
	gpuErrchk(cudaMemcpyFromSymbol(d_CECUTR, dg_CECUTR_, sizeof(hd_CECUTR)));
	gpuErrchk(cudaMemcpy(CECUTR_.ECUTR,d_CECUTR->ECUTR, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));

	//CESIAC
	gpuErrchk(cudaMemcpyFromSymbol(d_CESIAC, dg_CESIAC_, sizeof(hd_CESIAC)));
	gpuErrchk(cudaMemcpy(CESIAC_.ESIAC,d_CESIAC->ESIAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CESIAC_.IESI,d_CESIAC->IESI, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CESIAC_.NESI,d_CESIAC->NESI, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CEINAC
	gpuErrchk(cudaMemcpyFromSymbol(d_CEINAC, dg_CEINAC_, sizeof(hd_CEINAC)));
	gpuErrchk(cudaMemcpy(CEINAC_.EINAC,d_CEINAC->EINAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEINAC_.IEIN,d_CEINAC->IEIN, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CEINAC_.NEIN,d_CEINAC->NEIN, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CBRANG
	gpuErrchk(cudaMemcpyFromSymbol(d_CBRANG, dg_CBRANG_, sizeof(hd_CBRANG)));
	gpuErrchk(cudaMemcpy(CBRANG_.BET,d_CBRANG->BET, sizeof(double)*6, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CBRANG_.BK,d_CBRANG->BK, sizeof(double)*21, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CBRANG_.BP1,d_CBRANG->BP1, sizeof(double)*4*21*6*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CBRANG_.BP2,d_CBRANG->BP2, sizeof(double)*4*21*6*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CBRANG_.ZBEQ,d_CBRANG->ZBEQ, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));

	//CRELAX
	gpuErrchk(cudaMemcpyFromSymbol(d_CRELAX, dg_CRELAX_, sizeof(hd_CRELAX)));
	gpuErrchk(cudaMemcpy(CRELAX_.P,d_CRELAX->P, sizeof(double)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.ET,d_CRELAX->ET, sizeof(double)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.F,d_CRELAX->F, sizeof(double)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.IS0, d_CRELAX->IS0,sizeof(int)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.IS1,d_CRELAX->IS1, sizeof(int)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.IS2,d_CRELAX->IS2, sizeof(int)*NRX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.IFIRST,d_CRELAX->IFIRST, sizeof(int)*16*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.ILAST,d_CRELAX->ILAST, sizeof(int)*16*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.NCUR,&d_CRELAX->NCUR, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.KS,&d_CRELAX->KS, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CRELAX_.MODER,&d_CRELAX->MODER, sizeof(int), cudaMemcpyDeviceToHost));

	//CGRA01
	gpuErrchk(cudaMemcpyFromSymbol(d_CGRA01, dg_CGRA01_, sizeof(hd_CGRA01)));
	gpuErrchk(cudaMemcpy(CGRA01_.FF,d_CGRA01->FF, sizeof(double)*NQ*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA01_.ERA,d_CGRA01->ERA, sizeof(double)*NEX, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA01_.XSRA,d_CGRA01->XSRA, sizeof(double)*NEX*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA01_.IED,d_CGRA01->IED, sizeof(int)*NEGP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA01_.IEU,d_CGRA01->IEU, sizeof(int)*NEGP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA01_.NE,&d_CGRA01->NE, sizeof(int), cudaMemcpyDeviceToHost));

	//CGRA02
	gpuErrchk(cudaMemcpyFromSymbol(d_CGRA03, dg_CGRA03_, sizeof(hd_CGRA03)));
	gpuErrchk(cudaMemcpy(CGRA03_.QRA,d_CGRA03->QRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.PRA,d_CGRA03->PRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.DPRA,d_CGRA03->DPRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.ARA,d_CGRA03->ARA, sizeof(double)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.BRA,d_CGRA03->BRA, sizeof(double)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.PMAX, d_CGRA03->PMAX,sizeof(double)*NEGP*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.ITLRA,d_CGRA03->ITLRA, sizeof(int)*NP2*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGRA03_.ITURA,d_CGRA03->ITURA, sizeof(int)*NP2*MAXMAT, cudaMemcpyDeviceToHost));

	//CGPH00
	gpuErrchk(cudaMemcpyFromSymbol(d_CGPH00, dg_CGPH00_, sizeof(hd_CGPH00)));
	gpuErrchk(cudaMemcpy(CGPH00_.EPH,d_CGPH00->EPH, sizeof(double)*NTP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPH00_.XPH,d_CGPH00->XPH, sizeof(double)*NTP*17, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPH00_.IPHF,d_CGPH00->IPHF, sizeof(int)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPH00_.IPHL,d_CGPH00->IPHL, sizeof(int)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPH00_.NPHS,d_CGPH00->NPHS, sizeof(int)*99, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPH00_.NCUR,&d_CGPH00->NCUR, sizeof(int), cudaMemcpyDeviceToHost));

	//CGPP00
	gpuErrchk(cudaMemcpyFromSymbol(d_CGPP00, dg_CGPP00_, sizeof(hd_CGPP00)));
	gpuErrchk(cudaMemcpy(CGPP00_.ZEQPP,d_CGPP00->ZEQPP, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPP00_.F0,d_CGPP00->F0, sizeof(double)*MAXMAT*2, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CGPP00_.BCB,d_CGPP00->BCB, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));

	//CGPP01
	gpuErrchk(cudaMemcpyFromSymbol(d_CGPP01, dg_CGPP01_, sizeof(hd_CGPP01)));
	gpuErrchk(cudaMemcpy(CGPP01_.TRIP, d_CGPP01->TRIP,sizeof(double)*MAXMAT*NEGP, cudaMemcpyDeviceToHost));

	//CPELDB
	gpuErrchk(cudaMemcpyFromSymbol(d_CPELDB, dg_CPELDB_, sizeof(hd_CPELDB)));
	gpuErrchk(cudaMemcpy(CPELDB_.XSP,d_CPELDB->XSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPELDB_.PSP,d_CPELDB->PSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPELDB_.ASP,d_CPELDB->ASP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPELDB_.BSP,d_CPELDB->BSP, sizeof(double)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPELDB_.ITLP,d_CPELDB->ITLP, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPELDB_.ITUP,d_CPELDB->ITUP, sizeof(int)*MAXMAT*NEGP*NP, cudaMemcpyDeviceToHost));

	//CPINAC
	gpuErrchk(cudaMemcpyFromSymbol(d_CPINAC, dg_CPINAC_, sizeof(hd_CPINAC)));
	gpuErrchk(cudaMemcpy(CPINAC_.PINAC,d_CPINAC->PINAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPINAC_.IPIN,d_CPINAC->IPIN, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPINAC_.NPIN,d_CPINAC->NPIN, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));
	
	//CPSIAC
	gpuErrchk(cudaMemcpyFromSymbol(d_CPSIAC, dg_CPSIAC_, sizeof(hd_CPSIAC)));
	gpuErrchk(cudaMemcpy(CPSIAC_.PSIAC,d_CPSIAC->PSIAC, sizeof(double)*MAXMAT*NEGP*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPSIAC_.IPSI,d_CPSIAC->IPSI, sizeof(int)*MAXMAT*NO, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CPSIAC_.NPSI,d_CPSIAC->NPSI, sizeof(int)*MAXMAT, cudaMemcpyDeviceToHost));

	//CENANG
	gpuErrchk(cudaMemcpyFromSymbol(d_CENANG, dg_CENANG_, sizeof(hd_CENANG)));
	gpuErrchk(cudaMemcpy(CENANG_.EL,&d_CENANG->EL, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.EU, &d_CENANG->EU,sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.THL,&d_CENANG->THL, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.THU,&d_CENANG->THU, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.BSE,&d_CENANG->BSE, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.RBSE,&d_CENANG->RBSE, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.BSTH,&d_CENANG->BSTH, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.RBSTH,&d_CENANG->RBSTH, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.BSPH,&d_CENANG->BSPH, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.RBSPH,&d_CENANG->RBSPH, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDE,d_CENANG->PDE, sizeof(double)*2*3*NBEM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDE2,d_CENANG->PDE2, sizeof(double)*2*3*NBEM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDEP,d_CENANG->PDEP, sizeof(double)*2*3*NBEM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDA,d_CENANG->PDA, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDA2,d_CENANG->PDA2, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.PDAP,d_CENANG->PDAP, sizeof(double)*3*NBPHM*NBTHM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.LPDE,d_CENANG->LPDE, sizeof(int)*3*2*NBEM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.LPDA,d_CENANG->LPDA, sizeof(int)*3*NBPHM*NBTHM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.NE,&d_CENANG->NE, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.NTH,&d_CENANG->NTH, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.NPH,&d_CENANG->NPH, sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.LLE,&d_CENANG->LLE, sizeof(bool), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENANG_.LLTH,&d_CENANG->LLTH, sizeof(bool), cudaMemcpyDeviceToHost));

	//CENDET
	gpuErrchk(cudaMemcpyFromSymbol(d_CENDET, dg_CENDET_, sizeof(hd_CENDET)));
	gpuErrchk(cudaMemcpy(CENDET_.EL,d_CENDET->EL, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.EU,d_CENDET->EU, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.BSE,d_CENDET->BSE, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.RBSE,d_CENDET->RBSE, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.EDEP,d_CENDET->EDEP, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.EDEP2,d_CENDET->EDEP2, sizeof(double)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.DET,d_CENDET->DET, sizeof(double)*NIDM*NBEM2, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.NE,d_CENDET->NE, sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.NID,d_CENDET->NID, sizeof(int)*NIDM, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(CENDET_.LLE,&d_CENDET->LLE, sizeof(bool), cudaMemcpyDeviceToHost));

	//PENELOPE_MOD
	gpuErrchk(cudaMemcpyFromSymbol(d_PENELOPE_mod, dg_PENELOPE_mod_, sizeof(hd_PENELOPE_MOD)));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.EABS,d_PENELOPE_mod->EABS, sizeof(double)*3*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.C1,d_PENELOPE_mod->C1, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.C2,d_PENELOPE_mod->C2, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.WCC,d_PENELOPE_mod->WCC, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.WCR,d_PENELOPE_mod->WCR, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.DEN,d_PENELOPE_mod->DEN, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.RDEN,d_PENELOPE_mod->RDEN, sizeof(double)*MAXMAT, cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.E0STEP, &d_PENELOPE_mod->E0STEP,sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.DESOFT,&d_PENELOPE_mod->DESOFT, sizeof(double), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy( PENELOPE_mod_.SSOFT,&d_PENELOPE_mod->SSOFT,  sizeof(double), cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(PENELOPE_mod_.NMS,&d_PENELOPE_mod->NMS,   sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy(PENELOPE_mod_.NEGP,&d_PENELOPE_mod->NEGP,   sizeof(int), cudaMemcpyDeviceToHost));
	gpuErrchk(cudaMemcpy( PENELOPE_mod_.NMAT,&d_PENELOPE_mod->NMAT,  sizeof(int), cudaMemcpyDeviceToHost));

}

void memoryFreeGPU(){

    gpuErrchk(cudaFree(d_CEELDB));
    gpuErrchk(cudaFree(d_SECST));
    gpuErrchk(cudaFree(d_PENGEOM_MOD));
    gpuErrchk(cudaFree(d_CGCONE));
    gpuErrchk(cudaFree(d_QSURF));
    gpuErrchk(cudaFree(d_QTREE));
    gpuErrchk(cudaFree(d_RSEED));
    gpuErrchk(cudaFree(d_CIMDET));
    gpuErrchk(cudaFree(d_CERSEC));
    gpuErrchk(cudaFree(d_CEGRID));
	gpuErrchk(cudaFree(d_CJUMP1));
	gpuErrchk(cudaFree(d_CEIMFP));
	gpuErrchk(cudaFree(d_CPIMFP));
	gpuErrchk(cudaFree(d_CJUMP0));
	gpuErrchk(cudaFree(d_CGIMFP));
	gpuErrchk(cudaFree(d_CHIST));
	gpuErrchk(cudaFree(d_COMPOS));
	gpuErrchk(cudaFree(d_CADATA));
	gpuErrchk(cudaFree(d_CEIN));
	gpuErrchk(cudaFree(d_CGCO));
	gpuErrchk(cudaFree(d_CEBR));
	gpuErrchk(cudaFree(d_CELSEP));
	gpuErrchk(cudaFree(d_CECUTR));
	gpuErrchk(cudaFree(d_CESIAC));
	gpuErrchk(cudaFree(d_CEINAC));
	gpuErrchk(cudaFree(d_CBRANG));
	gpuErrchk(cudaFree(d_CRELAX));
	gpuErrchk(cudaFree(d_CGRA01));
	gpuErrchk(cudaFree(d_CGRA03));
	gpuErrchk(cudaFree(d_CGPH00));
	gpuErrchk(cudaFree(d_CGPP00));
	gpuErrchk(cudaFree(d_CGPP01));
	gpuErrchk(cudaFree(d_CPELDB));
	gpuErrchk(cudaFree(d_CPSIAC));
	gpuErrchk(cudaFree(d_CPINAC));
	gpuErrchk(cudaFree(d_CENANG));
	gpuErrchk(cudaFree(d_CENDET));
	gpuErrchk(cudaFree(d_PENELOPE_mod));

}





