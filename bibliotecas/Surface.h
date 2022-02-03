
#ifndef SURFACE_H
#define SURFACE_H

#include <string>
#include "DataTypes.h"


/**
  * class Surface
  * 
  */

class Surface
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Surface();

  /**
   * Empty Destructor
   */
  virtual ~Surface();

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

  // Forma Explícita
  // Termo Independente
  double a0;
  // Forma Explícita
  // AX, AY, AZ
  double3 a1;
  // Forma Explícita
  // AXX, AYY, AZZ
  double3 a2;
  // Forma Implícita
  // i1, i2, i3, i4, i5
  double5 indices;
  // Vetor Deslocamento
  // X-Shift, Y-Shift, Z-Shift
  double3 shifts;
  // Fatores de Escala
  // X-Scale, Y-Scale, Z-Scale
  double3 scales;
  // Angulos Euler
  // OMEGA, THETA, PHI
  double3 angules;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of a0
   * Forma Explícita
   * Termo Independente
   * @param value the new value of a0
   */
  void setA0(double value)
  {
    a0 = value;
  }

  /**
   * Get the value of a0
   * Forma Explícita
   * Termo Independente
   * @return the value of a0
   */
  double getA0()
  {
    return a0;
  }

  /**
   * Set the value of a1
   * Forma Explícita
   * AX, AY, AZ
   * @param value the new value of a1
   */
  void setA1(double3 value)
  {
    a1 = value;
  }

  /**
   * Get the value of a1
   * Forma Explícita
   * AX, AY, AZ
   * @return the value of a1
   */
  double3 getA1()
  {
    return a1;
  }

  /**
   * Set the value of a2
   * Forma Explícita
   * AXX, AYY, AZZ
   * @param value the new value of a2
   */
  void setA2(double3 value)
  {
    a2 = value;
  }

  /**
   * Get the value of a2
   * Forma Explícita
   * AXX, AYY, AZZ
   * @return the value of a2
   */
  double3 getA2()
  {
    return a2;
  }

  /**
   * Set the value of indices
   * Forma Implícita
   * i1, i2, i3, i4, i5
   * @param value the new value of indices
   */
  void setIndices(double5 value)
  {
    indices = value;
  }

  /**
   * Get the value of indices
   * Forma Implícita
   * i1, i2, i3, i4, i5
   * @return the value of indices
   */
  double5 getIndices()
  {
    return indices;
  }

  /**
   * Set the value of shifts
   * Vetor Deslocamento
   * X-Shift, Y-Shift, Z-Shift
   * @param value the new value of shifts
   */
  void setShifts(double3 value)
  {
    shifts = value;
  }

  /**
   * Get the value of shifts
   * Vetor Deslocamento
   * X-Shift, Y-Shift, Z-Shift
   * @return the value of shifts
   */
  double3 getShifts()
  {
    return shifts;
  }

  /**
   * Set the value of scales
   * Fatores de Escala
   * X-Scale, Y-Scale, Z-Scale
   * @param value the new value of scales
   */
  void setScales(double3 value)
  {
    scales = value;
  }

  /**
   * Get the value of scales
   * Fatores de Escala
   * X-Scale, Y-Scale, Z-Scale
   * @return the value of scales
   */
  double3 getScales()
  {
    return scales;
  }

  /**
   * Set the value of angules
   * Angulos Euler
   * OMEGA, THETA, PHI
   * @param value the new value of angules
   */
  void setAngules(double3 value)
  {
    angules = value;
  }

  /**
   * Get the value of angules
   * Angulos Euler
   * OMEGA, THETA, PHI
   * @return the value of angules
   */
  double3 getAngules()
  {
    return angules;
  }

  void initAttributes();

};

#endif // SURFACE_H
