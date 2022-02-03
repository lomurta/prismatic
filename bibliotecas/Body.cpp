#include "Body.h"





// Constructors/Destructors
//  



Body::Body()
{
  initAttributes();
}

Body::~Body()
{
}

//  
// Methods
//  


// Accessor methods
//  


// Other methods
//  

void Body::initAttributes()
{


}


/*void Body::locate(Particle& p)
{
//void Body::locate()
  //  {
    double fuzz = 0.0;
    double fuzzl = 0.0;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    
    for (int i = 0; i <=  p.qsurf.nsurf; i++){
        p.qtree.ksp[i] = 0;
    }

    int kb0 = p.qtree.nbodys;

    for (int kss = 1; kss <= p.qtree.ksurf[kb0][NXG]; kss++) {
        int ks = p.qtree.ksurf[kb0][kss];
        if ((p.qtree.ksp[kss] != 0) || (p.qtree.kflag[kb0][kss] > 4)) {

            //goto 101 
            for (int kbb = 1; kbb <= p.qtree.kdght[kb0][NXG]; kbb++) {
                int kb = p.qtree.kdght[kb0][kbb];

                for (int kss = 1; kss <= p.qtree.ksurf[kb][NXG]; kss++) {
                    int ks = p.qtree.ksurf[kb][kss];
                    int kf = p.qtree.kflag[kb][kss];
                    if ((kf < 3) && (p.qtree.ksp[ks] != kf)) {
                        //goto 102
                        p.track_mod.ibody = p.qtree.nbodys + 1;
                        p.track_mod.mat = 0;
                    }
                }

                if (kb == kb0) {
                    p.track_mod.ibody = kb; // a particula está dentro do corpo ou modulo kb
                    p.track_mod.mat = p.pengeom_mod.mater[kb]; // a particula está dentro do material kb
                }
                else {
                    if (p.qtree.kdght[kb][NXG] > 1) {
                        kb0 = kb; //o ponto está dentro de um submodulo
                        //goto100
                    }
                    else {
                        p.track_mod.ibody = kb; //a particula esta dentro de um corpo ou modulo irmão
                        p.track_mod.mat = p.pengeom_mod.mater[kb];
                    }
                }
            }
        }
        else {
            fsurf(ks, a, b, c, p);
            int absa = abs(a);
            if (absa > 1.0e-36) 
                fuzz = fuzzl * (b * b - 4.0 * a * c) / absa;
            else
                fuzz = fuzzl * abs(b);

            if (c < -fuzz)
                p.qtree.ksp[ks] = 1;
            else
                if (c > fuzz)
                    p.qtree.ksp[ks] = 2; 
                else
                    if (b < 0.0)
                        p.qtree.ksp[ks] = 1; //particula movendo-se para dentro
                    else
                        p.qtree.ksp[ks] = 2; // particula movendo-se para fora
           }
    }

}*/

void Body::locate()
{
    //void Body::locate()
      //  {
    double fuzz = 0.0;
    double fuzzl = 0.0;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    for (int i = 0; i <= qsurf.nsurf; i++) {
        qtree.ksp[i] = 0;
    }

    int kb0 = qtree.nbodys;

    for (int kss = 1; kss <= qtree.ksurf[kb0][NXG]; kss++) {
        int ks = qtree.ksurf[kb0][kss];
        if ((qtree.ksp[kss] != 0) || (qtree.kflag[kb0][kss] > 4)) {

            //goto 101 
            for (int kbb = 1; kbb <= qtree.kdght[kb0][NXG]; kbb++) {
                int kb = qtree.kdght[kb0][kbb];

                for (int kss = 1; kss <= qtree.ksurf[kb][NXG]; kss++) {
                    int ks = qtree.ksurf[kb][kss];
                    int kf = qtree.kflag[kb][kss];
                    if ((kf < 3) && (qtree.ksp[ks] != kf)) {
                        //goto 102
                        track_mod.ibody = qtree.nbodys + 1;
                        track_mod.mat = 0;
                    }
                }

                if (kb == kb0) {
                    track_mod.ibody = kb; // a particula está dentro do corpo ou modulo kb
                    track_mod.mat = pengeom_mod.mater[kb]; // a particula está dentro do material kb
                }
                else {
                    if (qtree.kdght[kb][NXG] > 1) {
                        kb0 = kb; //o ponto está dentro de um submodulo
                        //goto100
                    }
                    else {
                        track_mod.ibody = kb; //a particula esta dentro de um corpo ou modulo irmão
                        track_mod.mat = pengeom_mod.mater[kb];
                    }
                }
            }
        }
        else {
            fsurf(ks, a, b, c);
            int absa = abs(a);
            if (absa > 1.0e-36)
                fuzz = fuzzl * (b * b - 4.0 * a * c) / absa;
            else
                fuzz = fuzzl * abs(b);

            if (c < -fuzz)
                qtree.ksp[ks] = 1;
            else
                if (c > fuzz)
                    qtree.ksp[ks] = 2;
                else
                    if (b < 0.0)
                        qtree.ksp[ks] = 1; //particula movendo-se para dentro
                    else
                        qtree.ksp[ks] = 2; // particula movendo-se para fora
        }
    }

}

void Body::fsurf(int ks, double& a, double& b, double& c)
{

    if (qsurf.kplane[ks] == 0) {
        a = track_mod.u * (qsurf.axx[ks] * track_mod.u + qsurf.axy[ks] * track_mod.v + qsurf.axz[ks] * track_mod.w) +
            track_mod.v * (qsurf.ayy[ks] * track_mod.v + qsurf.ayz[ks] * track_mod.w) + track_mod.w * qsurf.azz[ks] * track_mod.w;
        double XXX = qsurf.axx[ks] * track_mod.x + qsurf.axy[ks] * track_mod.y + qsurf.axz[ks] * track_mod.z + qsurf.ax[ks];
        double YYY = qsurf.ayy[ks] * track_mod.y + qsurf.ayz[ks] * track_mod.z + qsurf.ay[ks];
        double ZZZ = qsurf.azz[ks] * track_mod.z + qsurf.az[ks];

        b = track_mod.u * (qsurf.axx[ks] * track_mod.x + track_mod.v * qsurf.axy[ks] * track_mod.x + qsurf.ayy[ks] * track_mod.y + YYY) +
            track_mod.w * (qsurf.axz[ks] * track_mod.x + track_mod.x * qsurf.ayz[ks] * track_mod.y) + track_mod.y * qsurf.azz[ks] * ZZZ;

        c = track_mod.x * XXX + track_mod.y * YYY + track_mod.z * ZZZ + qsurf.a0[ks];
    }
    else
        a = 0.0;
    b = track_mod.u * qsurf.ax[ks] + track_mod.v * qsurf.ay[ks] + track_mod.w * qsurf.az[ks];
    c = track_mod.x * qsurf.ax[ks] + track_mod.y * qsurf.ay[ks] + track_mod.z * qsurf.az[ks] + qsurf.a0[ks];
}



/*void Body::fsurf(int ks, double &a, double &b, double &c, PENGEOM_MOD* pengeom_mod, TRACK_MOD* track_mod, QSURF* qsurf, QTREE* qtree)
{

    if (qsurf->kplane[ks] == 0) {
        a = track_mod->u * (qsurf->axx[ks] * track_mod->u + qsurf->axy[ks] * track_mod->v + qsurf->axz[ks] * track_mod->w ) +
            track_mod->v * (qsurf->ayy[ks] * track_mod->v + qsurf->ayz[ks] * track_mod->w) + track_mod->w* qsurf->azz[ks] * track_mod->w;
        double XXX = qsurf->axx[ks] * track_mod->x + qsurf->axy[ks] * track_mod->y + qsurf->axz[ks] * track_mod->z + qsurf->ax[ks];
        double YYY = qsurf->ayy[ks] * track_mod->y + qsurf->ayz[ks] * track_mod->z + qsurf->ay[ks];
        double ZZZ = qsurf->azz[ks] * track_mod->z + qsurf->az[ks];

        b = track_mod->u * (qsurf->axx[ks] * track_mod->x + track_mod->v*qsurf->axy[ks] * track_mod->x + qsurf->ayy[ks] * track_mod->y + YYY) +
            track_mod->w * (qsurf->axz[ks] * track_mod->x + track_mod->x*qsurf->ayz[ks] * track_mod->y) + track_mod->y * qsurf->azz[ks] * ZZZ;

        c = track_mod->x * XXX + track_mod->y * YYY + track_mod->z * ZZZ + qsurf->a0[ks];
    }
    else
        a = 0.0;
        b = track_mod->u * qsurf->ax[ks] + track_mod->v * qsurf->ay[ks] + track_mod->w * qsurf->az[ks];
        c = track_mod->x * qsurf->ax[ks] + track_mod->y * qsurf->ay[ks] + track_mod->z * qsurf->az[ks] + qsurf->a0[ks];
}

void step(PENGEOM_MOD &pengeom_mod, PENELOPE_MOD &penelope_mod, QSURF &qsurf, QTREE &qtree, TRACK_MOD &track_mod, double ds, double dsef, int ncross) {

    //variaveis locais

    STEP_MOD step_mod;


    dsef = 0.0;
    pengeom_mod.dstot = 0.0;
    ncross = 0;
    pengeom_mod.kslast = 0;

    for (int i = 1; i <= qsurf.nsurf; i++) {
        qtree.ksp[i] = 0; //side pointers laterais das superficies avaliadas
    }

    step_mod.mat0 = track_mod.mat; //material inicial

    if (track_mod.mat == 0)
        step_mod.dsres = 1e35; // no vacuo, as particulas voam livremente
    else step_mod.dsres = ds; // comprimento restante do caminho

    //A particula entra de fora do involucro

    if (track_mod.ibody == qtree.nbodys) {
        step_mod.kb1 = qtree.nbodys;
        stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
        if (step_mod.nsc == 0) {
            step_goto300(pengeom_mod, track_mod, qtree, step_mod, dsef);
            return;
        }
        step_mod.nsct = step_mod.nsc;
        step_mod.nst = *qtree.ksurf[step_mod.kb1, NXG];
        for (int ki = step_mod.nsct; ki >= 1; ki--) {
        // a particula atravessa uma superficie
            pengeom_mod.kslast = step_mod.is[ki];
            if (qtree.ksp[pengeom_mod.kslast] == 1)
                qtree.ksp[pengeom_mod.kslast] = 2;
            else
                qtree.ksp[pengeom_mod.kslast] = 1;
            step_mod.dsp = step_mod.s[ki];
            dsef = dsef + step_mod.dsp;
            pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsp;
            track_mod.x = track_mod.x + step_mod.dsp * track_mod.u;
            track_mod.y = track_mod.y + step_mod.dsp * track_mod.v;
            track_mod.z = track_mod.z + step_mod.dsp * track_mod.w;
            step_mod.nsc = step_mod.nsc - 1;
            if (step_mod.nsc > 0)
                for (int i = 1; i <= step_mod.nsc; i++) {
                    step_mod.s[i] = step_mod.s[i] - step_mod.dsp;
                }
            for (int kss = 1; kss <= step_mod.nst; kss++) {
                step_mod.ks1 = qtree.ksurf[step_mod.kb1, kss];
                step_mod.kf = qtree.kflag[step_mod.kb1, kss];
                if ((step_mod.kf < 3) && (qtree.ksp[step_mod.ks1] != step_mod.kf)) {
                    step_goto101(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
                    return;
                } 
            }
            step_goto100(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        }
    }

    //travessias de superficices
    step_mod.ibodyl = track_mod.ibody;
    if (pengeom_mod.lverb)
        step_mod.nerr = 0;
    step_goto102(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);


}

void step_goto100(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {
    // A particula entra no involucro
    steplb(step_mod.kb1, step_mod.ierr);
    //A particula entrou no submodulo
    if (step_mod.ierr == -1) {
        step_mod.kb1 = track_mod.ibody;
        stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
        step_goto100(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        return;
    }
    else {
        // A particula entra em um corpo material
        if (track_mod.mat != 0) {
            ncross = 1;
        }
        else {
            step_mod.kb1 = track_mod.ibody;
            stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
            step_goto200(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        }
    }
}

void step_goto101(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {
    step_goto300(pengeom_mod, track_mod, qtree, step_mod, dsef);
}

void step_goto102(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {

    step_mod.kb1 = track_mod.ibody;
    stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
    steplb(step_mod.kb1, step_mod.ierr);

    //evidencia os erros de arredondamentos
    if (step_mod.ierr != 0) {
        if (step_mod.nsc > 0) {
            //Quando uma superfície está muito próxima, movemos a partícula além dela
            if (step_mod.s[step_mod.nsc] < 1e-10) {
                pengeom_mod.kslast = step_mod.is[step_mod.nsc];
                if (qtree.ksp[pengeom_mod.kslast] == 1) {
                    qtree.ksp[pengeom_mod.kslast] = 2;
                }
                else {
                    qtree.ksp[pengeom_mod.kslast] = 1;
                }
                step_mod.dsp = step_mod.s[step_mod.nsc];
                track_mod.x = track_mod.x + step_mod.dsp * track_mod.u;
                track_mod.y = track_mod.y + step_mod.dsp * track_mod.v;
                track_mod.z = track_mod.z + step_mod.dsp * track_mod.w;
                if (track_mod.mat == step_mod.mat0) {
                    dsef = dsef + step_mod.dsp;
                    step_mod.dsres = step_mod.dsres - step_mod.dsp;
                }
                pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsp;
                step_mod.nsc = step_mod.nsc - 1;
                if (track_mod.ibody <= qtree.nbodys)
                    step_goto102(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);

            }
        }
        if (pengeom_mod.lverb) {
            step_mod.nerr = step_mod.nerr + 1;
            if ((qtree.nwarn < 100) && (track_mod.mat != 0)) {
                //Escrita de detalhes no arquivo de saida comandos write

                for (int kss = 1; kss <= qtree.ksurf[step_mod.kb1, NXG]; kss++) {
                    step_mod.ks = qtree.ksurf[step_mod.kb1, kss];
                    step_mod.kflo = qtree.kflag[step_mod.kb1, kss];
                    if (step_mod.kflo < 3) {
                        for (int ki = step_mod.nsc; ki >= 1; ki--) {
                            if (step_mod.ks == step_mod.is[ki]) {
                                step_mod.sw = step_mod.s[ki];
                                step_goto103(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
                            }
                        }
                        step_mod.sw = 0;
                        step_goto103(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
                    }

                }

                qtree.nwarn = qtree.nwarn + 1;
            }

        }
        if (track_mod.ibody <= qtree.nbodys)
            step_goto102(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
    }

    //A partícula permanece no mesmo material.

    if ((track_mod.mat != 0) && (step_mod.dsres < step_mod.s[step_mod.nsc])) {
        if (track_mod.mat == step_mod.mat0)
            dsef = dsef + step_mod.dsres;
        pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsres;
        track_mod.x = track_mod.x + step_mod.dsres * track_mod.u;
        track_mod.y = track_mod.y + step_mod.dsres * track_mod.v;
        track_mod.z = track_mod.z + step_mod.dsres * track_mod.w;
    }

    //Nova posição
    step_goto200(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);


}

void step_goto103(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {
    fsurf(step_mod.kb1);
    if (step_mod.kflo == qtree.ksp[step_mod.ks]) {
        //escreve na linha
    }
    else {
        //escreve na linha
        pengeom_mod.kslast = step_mod.ks;

    }

}



void step_goto200(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {

    if (step_mod.nsc == 0) {
        if (track_mod.mat == step_mod.mat0)
            dsef = dsef + step_mod.dsres;
        pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsres;
        track_mod.x = track_mod.x + step_mod.dsres * track_mod.u;
        track_mod.y = track_mod.y + step_mod.dsres * track_mod.v;
        track_mod.z = track_mod.z + step_mod.dsres * track_mod.w;
    }
    step_mod.nsct = step_mod.nsc;
    step_mod.matl = track_mod.mat;
    step_mod.ibodyl = track_mod.ibody;
    for (int ki = step_mod.nsct; ki >= 1; ki--) {
       // A etapa termina dentro do corpo.
        if (step_mod.dsres < step_mod.s[ki]) {
            if (track_mod.mat == step_mod.mat0)
                dsef = dsef + step_mod.dsres;
            pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsres;
            track_mod.x = track_mod.x + step_mod.dsres * track_mod.u;
            track_mod.y = track_mod.y + step_mod.dsres * track_mod.v;
            track_mod.z = track_mod.z + step_mod.dsres * track_mod.w;
            //A partícula atravessa uma superfície.
            pengeom_mod.kslast = step_mod.is[ki];
            if (qtree.ksp[pengeom_mod.kslast] == 1)
                qtree.ksp[pengeom_mod.kslast] == 2;
            else qtree.ksp[pengeom_mod.kslast] == 1;
            track_mod.x = track_mod.x + step_mod.dsp * track_mod.u;
            track_mod.y = track_mod.y + step_mod.dsp * track_mod.v;
            track_mod.z = track_mod.z + step_mod.dsp * track_mod.w;

            if (track_mod.mat == step_mod.mat0) {
                dsef = dsef + step_mod.dsp;
                step_mod.dsres = step_mod.dsres - step_mod.dsp;
            }

            pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsp;
            step_mod.nsc = step_mod.nsc - 1;
            if (step_mod.nsc > 0)
                for (int i = 1; i <= step_mod.nsc; i++) {
                    step_mod.s[i] = step_mod.s[i] - step_mod.dsp;
                }

            steplb(step_mod.kb1, step_mod.ierr);
            step_goto201(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        }
    }

    step_goto300(pengeom_mod, track_mod, qtree, step_mod, dsef);
}


void step_goto201(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {

    step_mod.kb1 = track_mod.ibody;
    if (step_mod.ierr == -1) {
        //A particula entra em um submodulo.
        stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
        steplb(step_mod.kb1, step_mod.ierr);
        step_goto201(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        return;
    }
    else 
        if (step_mod.ierr == 1) {
            //a partiula deixa o corpo ou o módulo
            if (track_mod.ibody <= 1) {
                stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
                steplb(step_mod.kb1, step_mod.ierr);
                step_goto201(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
                return;
            }
            else {
                // a particula sai do recinto
                if (track_mod.mat != step_mod.matl)
                    ncross = ncross + 1;
                step_goto300(pengeom_mod, track_mod, qtree, step_mod, dsef);
                return;
            }
        }
    // A partícula continua voando quando entra em uma região vazia
    if (track_mod.mat == 0) {
        if (step_mod.matl == step_mod.mat0)
            ncross = ncross + 1;
        step_mod.matl = 0;
        step_mod.dsres = 1e35;
        step_goto202(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
        return;
    }
    //A partícula continua voando quando entra em um novo corpo do mesmo material que não faz parte de um detector diferente
    else {
        if (track_mod.mat == step_mod.matl) {
            if (pengeom_mod.kdet[track_mod.ibody] == pengeom_mod.kdet[step_mod.ibodyl]) {
                step_goto202(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
                return;
            }
            else {
                ncross = ncross + 1;
                return;
            }
        }
        // e para quando penetra um novo corpo material ou um Detector.
        else {
            ncross = ncross + 1;
            return;
        }
    }
    step_goto202(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);   
}

void step_goto202(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, int& ncross, double& dsef) {
    stepsi(step_mod.kb1, step_mod.s, step_mod.is, step_mod.nsc);
    step_goto200(pengeom_mod, track_mod, qtree, step_mod, ncross, dsef);
    //final do laço DO


}


void step_goto300(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, STEP_MOD& step_mod, double& dsef) {
        //a particula sai do recinto
           step_mod.dsp = 1e36;
           track_mod.ibody = qtree.nbodys + 1;
           track_mod.mat = 0;
           if (track_mod.mat == step_mod.mat0)
               dsef = dsef + step_mod.dsp;
           pengeom_mod.dstot = pengeom_mod.dstot + step_mod.dsp;
           track_mod.x = track_mod.x + step_mod.dsp * track_mod.u;
           track_mod.y = track_mod.y + step_mod.dsp * track_mod.v;
           track_mod.z = track_mod.z + step_mod.dsp * track_mod.w;
           //return fim da funcao
}

void fsurf(int& kb1){

}

void steplb(int& kb1, int& ierr) {

    // Analise o corpo ou modulo atual

    qbody.





}

void stepsi(int& kb1, int* s, int* is, int& nsc) {


}*/




