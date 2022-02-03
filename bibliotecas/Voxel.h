
#ifndef VOXEL_H
#define VOXEL_H

#include <string>
#include "DataTypes.h"


/**
  * class Voxel
  * 
  */

class Voxel
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Voxel();

  /**
   * Empty Destructor
   */
  virtual ~Voxel();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * Escreve os resultados da distribuição de dose e a dose nos corpos
   */
  void sdose()
  {
  }

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

  // Energia depositada
  double eD;
  // Coordenadas de posição em cm do voxel
  // r= (x,y,z)
  double3 positon;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of eD
   * Energia depositada
   * @param value the new value of eD
   */
  void setED(double value)
  {
    eD = value;
  }

  /**
   * Get the value of eD
   * Energia depositada
   * @return the value of eD
   */
  double getED()
  {
    return eD;
  }

  /**
   * Set the value of positon
   * Coordenadas de posição em cm do voxel
   * r= (x,y,z)
   * @param value the new value of positon
   */
  void setPositon(double3 value)
  {
    positon = value;
  }

  /**
   * Get the value of positon
   * Coordenadas de posição em cm do voxel
   * r= (x,y,z)
   * @return the value of positon
   */
  double3 getPositon()
  {
    return positon;
  }

  void initAttributes();

};

#endif // VOXEL_H
