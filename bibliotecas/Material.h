#ifndef MATERIAL_H
#define MATERIAL_H

#include "DataTypes.h"

/**
  * class Material
  * 
  */

class Material
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Material();

  /**
   * Empty Destructor
   */
  virtual ~Material();

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

  // Número do Material conforme tabelas disponivel pelo PENELOPE
  int number;
  // Deflexão angular média
  float c1;
  // Perda fracionaria máxima de energia
  float c2;
  // Perda de enregia para colisoes inelasticas duras
  float wcc;
  // Perda de energia para emissao de bremstrahlung
  float wcr;
  // Energia de absorção
  // (electron, photon, positron)
  double3 eAbs;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of number
   * Número do Material conforme tabelas disponivel pelo PENELOPE
   * @param value the new value of number
   */
  void setNumber(int value)
  {
    number = value;
  }

  /**
   * Get the value of number
   * Número do Material conforme tabelas disponivel pelo PENELOPE
   * @return the value of number
   */
  int getNumber()
  {
    return number;
  }

  /**
   * Set the value of c1
   * Deflexão angular média
   * @param value the new value of c1
   */
  void setC1(float value)
  {
    c1 = value;
  }

  /**
   * Get the value of c1
   * Deflexão angular média
   * @return the value of c1
   */
  float getC1()
  {
    return c1;
  }

  /**
   * Set the value of c2
   * Perda fracionaria máxima de energia
   * @param value the new value of c2
   */
  void setC2(float value)
  {
    c2 = value;
  }

  /**
   * Get the value of c2
   * Perda fracionaria máxima de energia
   * @return the value of c2
   */
  float getC2()
  {
    return c2;
  }

  /**
   * Set the value of wcc
   * Perda de enregia para colisoes inelasticas duras
   * @param value the new value of wcc
   */
  void setWcc(float value)
  {
    wcc = value;
  }

  /**
   * Get the value of wcc
   * Perda de enregia para colisoes inelasticas duras
   * @return the value of wcc
   */
  float getWcc()
  {
    return wcc;
  }

  /**
   * Set the value of wcr
   * Perda de energia para emissao de bremstrahlung
   * @param value the new value of wcr
   */
  void setWcr(float value)
  {
    wcr = value;
  }

  /**
   * Get the value of wcr
   * Perda de energia para emissao de bremstrahlung
   * @return the value of wcr
   */
  float getWcr()
  {
    return wcr;
  }

  /**
   * Set the value of eAbs
   * Energia de absorção
   * (electron, photon, positron)
   * @param value the new value of eAbs
   */
  void setEAbs(double3 value)
  {
    eAbs = value;
  }

  /**
   * Get the value of eAbs
   * Energia de absorção
   * (electron, photon, positron)
   * @return the value of eAbs
   */
  double3 getEAbs()
  {
    return eAbs;
  }

  void initAttributes();

};

#endif // MATERIAL_H
