#ifndef MODULE_H
#define MODULE_H

#include <string>
#include "DataTypes.h"
#include "Body.h"
#include "Material.h"
#include "Surface_BM.h"



/**
  * class Module
  * 
  */

class Module
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Module();

  /**
   * Empty Destructor
   */
  virtual ~Module();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * Possiblita a clonagem de um módulo
   */
  void clone()
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

  // Identificação
  int id;
  // Corpos
  //Body* bodys[10];
  // Material
  Material material;
  // Módulos
  // Superfices
  Surface_BM surfaces_BM;
  // Fatores de Escala
  // X-Shift, Y-Shift, Z-Shift
  double3 shifts;
  // Ângulos Euler
  // Omega, Theta, Phi
  double3 angules;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of id
   * Identificação
   * @param value the new value of id
   */
  void setId(int value)
  {
    id = value;
  }

  /**
   * Get the value of id
   * Identificação
   * @return the value of id
   */
  int getId()
  {
    return id;
  }

  /**
   * Set the value of bodys
   * Corpos
   * @param value the new value of bodys
   */
 /* void setBodys(Body value)
  {
    //bodys = value;
  }

  /**
   * Get the value of bodys
   * Corpos
   * @return the value of bodys
   */
 /* Body getBodys()
  {
   // return *bodys;
  }
  */
  /**
   * Set the value of material
   * Material
   * @param value the new value of material
   */
  void setMaterial(Material value)
  {
    this->material = value;
  }

  /**
   * Get the value of material
   * Material
   * @return the value of material
   */
  Material getMaterial()
  {
    return this->material;
  }


  /**
   * Set the value of surfaces_BM
   * Superfices
   * @param value the new value of surfaces_BM
   */
  void setSurfaces_BM(Surface_BM value)
  {
    surfaces_BM = value;
  }

  /**
   * Get the value of surfaces_BM
   * Superfices
   * @return the value of surfaces_BM
   */
  Surface_BM getSurfaces_BM()
  {
    return surfaces_BM;
  }

  /**
   * Set the value of shifts
   * Fatores de Escala
   * X-Shift, Y-Shift, Z-Shift
   * @param value the new value of shifts
   */
  void setShifts(double3 value)
  {
    shifts = value;
  }

  /**
   * Get the value of shifts
   * Fatores de Escala
   * X-Shift, Y-Shift, Z-Shift
   * @return the value of shifts
   */
  double3 getShifts()
  {
    return shifts;
  }

  /**
   * Set the value of angules
   * Ângulos Euler
   * Omega, Theta, Phi
   * @param value the new value of angules
   */
  void setAngules(double3 value)
  {
    angules = value;
  }

  /**
   * Get the value of angules
   * Ângulos Euler
   * Omega, Theta, Phi
   * @return the value of angules
   */
  double3 getAngules()
  {
    return angules;
  }

  void initAttributes();

};

#endif // MODULE_H
