
#ifndef ENERGY_SPECTRUM_H
#define ENERGY_SPECTRUM_H

#include <string>
#include "DataTypes.h"


/**
  * class Energy_Spectrum
  * 
  */

class Energy_Spectrum
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Energy_Spectrum();

  /**
   * Empty Destructor
   */
  virtual ~Energy_Spectrum();

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

  // Energia da extremidade inferior
  double ei;
  // Probabilidade relativa
  double pi;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of ei
   * Energia da extremidade inferior
   * @param value the new value of ei
   */
  void setEi(double value)
  {
    ei = value;
  }

  /**
   * Get the value of ei
   * Energia da extremidade inferior
   * @return the value of ei
   */
  double getEi()
  {
    return ei;
  }

  /**
   * Set the value of pi
   * Probabilidade relativa
   * @param value the new value of pi
   */
  void setPi(double value)
  {
    pi = value;
  }

  /**
   * Get the value of pi
   * Probabilidade relativa
   * @return the value of pi
   */
  double getPi()
  {
    return pi;
  }

  void initAttributes();

};

#endif // ENERGY_SPECTRUM_H
