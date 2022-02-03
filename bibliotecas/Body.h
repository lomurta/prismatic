#ifndef BODY_H
#define BODY_H

#include <string>

#include "Material.h"
#include "Surface_BM.h"
#include "DataTypes.h"
#include "Particle.h"



extern PENELOPE_MOD penelope_mod;
extern PENGEOM_MOD pengeom_mod;
extern TRACK_MOD track_mod;
extern QSURF qsurf;
extern QTREE qtree;

/**
  * class Body
  * 
  */
class Body
{
public:
  // Constructors/Destructors
  //  

  /**
   * Empty Constructor
   */
    Body();

  /**
   * Empty Destructor
   */
  virtual ~Body();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * Determina o corpo que contem as coordenadas de posição X, Y, Z.
   * 
   * Paramentros de Entrada
   * X, Y, Z: coordenadas da partícula
   * U, V, W: direção do movimento
   */




 //void locate(QSURF* qsurf, QTREE* qtree, PENGEOM_MOD* pengeom_mod, TRACK_MOD* track_mod);

 void locate(Particle& p);

 void locate();
 void fsurf(int ks, double& a, double& b, double& c);


 //void fsurf(int ks, double &a, double &b, double &c, PENGEOM_MOD* pengeom_mod, TRACK_MOD* track_mod, QSURF* qsurf, QTREE* qtree);
//  void locate();

  /**
   * Lida com a parte geometrica da simulação e o caminho percorrido pela particula.
   * 
   * Valores de Entrada
   * X, Y, Z: coordenadas do ponto inicial
   * U, V, W: cossenos diretores do deslocamento
   * IBODY: corpo onde o ponto inicial está localizado.
   * MAT: material do corpo IBODY
   * 
   * Parametros de Entrada:
   * DS: comprimento do caminho a percorrer
   * 
   * Parametros de Saida:
   * DSEF:  comprimento do caminho percorrido antes de deixar o material inicial ou
   * completar o salto.
   * NCROSS: 0 se toda o passo estiver contido no material inicial.
   * 
   * Valores de Saida
   * X, Y, Z: Coordenadas da posicão final
   * IBODY: corpo onde o ponto final está localizado
   * MAT: Material do corpo. MAT = 0, indica que a particula C escapou do sistema
   * estudado.
   * 
   * DSTOT: comprimento do caminho percorrido
   * KSLAST: Ultima superficie cruzada.
   */
 //void step(PENGEOM_MOD &pengeom_mod, PENELOPE_MOD &penelope_mod, QSURF &qsurf, QTREE &qtree, TRACK_MOD &track_md);
 //void step_goto300(PENGEOM_MOD& pengeom_mod, TRACK_MOD& track_mod, QTREE& qtree, double& dsp, double& dsef, int& mat0);
 //void steplb(int& kb1, int& ierr);
 //void stepsi(int& kb1, int* s,int* is,int& nsc);

protected:
  // Static Protected attributes
  //  

  // Protected attributes
  //  


  // Protected attribute accessor methods
  //  


  // Protected attribute accessor methods
  //

private:
  // Static Private attributes
  //  

  // Private attributes
  //  

  // Identificação do corpo
  int id;
  // Comprimento máximo do passo
  double dsMax;
  // Energia de absorçao local
  // (electron, photon, positron)
  double3 eAbsB;
  // Material que compõe o corpo
  Material material;
  // Superficies do corpo
  Surface_BM surfaces_BM;
  

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of id
   * Identificação do corpo
   * @param value the new value of id
   */
  void setId(int value)
  {
    id = value;
  }

  /**
   * Get the value of id
   * Identificação do corpo
   * @return the value of id
   */
  int getId()
  {
    return id;
  }

  /**
   * Set the value of dsMax
   * Comprimento máximo do passo
   * @param value the new value of dsMax
   */
  void setDsMax(double value)
  {
    dsMax = value;
  }

  /**
   * Get the value of dsMax
   * Comprimento máximo do passo
   * @return the value of dsMax
   */
  double getDsMax()
  {
    return dsMax;
  }

  /**
   * Set the value of eAbsB
   * Energia de absorçao local
   * (electron, photon, positron)
   * @param value the new value of eAbsB
   */
  void setEAbsB(double3 value)
  {
    eAbsB = value;
  }

  /**
   * Get the value of eAbsB
   * Energia de absorçao local
   * (electron, photon, positron)
   * @return the value of eAbsB
   */
  double3 getEAbsB()
  {
    return eAbsB;
  }

  /**
   * Set the value of material
   * Material que compõe o corpo
   * @param value the new value of material
   */
  void setMaterial(Material value)
  {
    material = value;
  }

  /**
   * Get the value of material
   * Material que compõe o corpo
   * @return the value of material
   */
  Material getMaterial()
  {
    return material;
  }

  /**
   * Set the value of surfaces_BM
   * Superficies do corpobor
   * @param value the new value of surfaces_BM
   */
  void setSurfaces_BM(Surface_BM value)
  {
    surfaces_BM = value;
  }

  /**
   * Get the value of surfaces_BM
   * Superficies do corpo
   * @return the value of surfaces_BM
   */


  /**
   * Get the value of surfaces_BM
   * Superficies do corpo
   * @return the value of surfaces_BM
   */
  Surface_BM getSurfaces_BM()
  {
      return surfaces_BM;
  }


  void initAttributes();

};

#endif // BODY_H
