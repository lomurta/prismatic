#ifndef PARTICLE_H
#define PARTICLE_H

#include "DataTypes.h"

/**
  * class Particle
  * 
  */

class Particle
{
public:

   /* PENELOPE_MOD pengeom_mod;
    TRACK_MOD track_mod;
    QSURF qsurf;
    QTREE qtree;
    QBODY qbody;
    STEP_MOD step_mod;*/


  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Particle();

  /**
   * Empty Destructor
   */
  virtual ~Particle();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //

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

  // Energia da particula
  double e;
  int wght;
  // Origem das partículas
  // ILB(1): 1 particulas primarias; 2 particula secundaria
  // ILB(2): tipo KParP particula pai, somente se ILB(1) > 1
  // ILB(3): mecanismo de interação ICOL que originou a particula, somente quando ILB(1) > 1
  // ILB(4): um valor nao zero identifica particulas emitidas a partir de eventos de relaxamento atomico.
  // ILB(5): rotulo utilizado pelo usuário 
  int5 ilb;
  // Coordenadas de posição em cm
  // r= (x,y,z)
  double3 position;
  // Direção cossenos de movimento
  // d=(u, v, w)
  double3 direction;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of e
   * Energia da particula
   * @param value the new value of e
   */
  void setE(double value)
  {
    e = value;
  }

  /**
   * Get the value of e
   * Energia da particula
   * @return the value of e
   */
  double getE()
  {
    return e;
  }

  /**
   * Set the value of wght
   * @param value the new value of wght
   */
  void setWght(int value)
  {
    wght = value;
  }

  /**
   * Get the value of wght
   * @return the value of wght
   */
  int getWght()
  {
    return wght;
  }

  /**
   * Set the value of ilb
   * Origem das partículas
   * ILB(1): 1 particulas primarias; 2 particula secundaria
   * ILB(2): tipo KParP particula pai, somente se ILB(1) > 1
   * ILB(3): mecanismo de interação ICOL que originou a particula, somente quando
   * ILB(1) > 1
   * ILB(4): um valor nao zero identifica particulas emitidas a partir de eventos de
   * relaxamento atomico.
   * ILB(5): rotulo utilizado pelo usuário
   * @param value the new value of ilb
   */
  void setIlb(int5 value)
  {
    ilb = value;
  }

  /**
   * Get the value of ilb
   * Origem das partículas
   * ILB(1): 1 particulas primarias; 2 particula secundaria
   * ILB(2): tipo KParP particula pai, somente se ILB(1) > 1
   * ILB(3): mecanismo de interação ICOL que originou a particula, somente quando
   * ILB(1) > 1
   * ILB(4): um valor nao zero identifica particulas emitidas a partir de eventos de
   * relaxamento atomico.
   * ILB(5): rotulo utilizado pelo usuário
   * @return the value of ilb
   */
  int5 getIlb()
  {
    return ilb;
  }

  /**
   * Set the value of position
   * Coordenadas de posição em cm
   * r= (x,y,z)
   * @param value the new value of position
   */
  void setPosition(double3 value)
  {
    position = value;
  }

  /**
   * Get the value of position
   * Coordenadas de posição em cm
   * r= (x,y,z)
   * @return the value of position
   */
  double3 getPosition()
  {
    return position;
  }

  /**
   * Set the value of direction
   * Direção cossenos de movimento
   * d=(u, v, w)
   * @param value the new value of direction
   */
  void setDirection(double3 value)
  {
    direction = value;
  }

  /**
   * Get the value of direction
   * Direção cossenos de movimento
   * d=(u, v, w)
   * @return the value of direction
   */
  double3 getDirection()
  {
    return direction;
  }

  void initAttributes();

};

#endif // PARTICLE_H
