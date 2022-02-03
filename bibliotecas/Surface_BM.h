#ifndef SURFACE_BM_H
#define SURFACE_BM_H

#include <string>
#include "Surface.h"

/**
  * class Surface_BM
  * 
  */

class Surface_BM
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Surface_BM();

  /**
   * Empty Destructor
   */
  virtual ~Surface_BM();

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

  // Superficies do corpo
  Surface surface;
  // -1: coordenadas menores
  //  1: coordenadas maiores
  int sidePointer;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of surface
   * Superficies do corpo
   * @param value the new value of surface
   */
  void setSurface(Surface value)
  {
    surface = value;
  }

  /**
   * Get the value of surface
   * Superficies do corpo
   * @return the value of surface
   */
  Surface getSurface()
  {
    return surface;
  }

  /**
   * Set the value of sidePointer
   * -1: coordenadas menores
   *  1: coordenadas maiores
   * @param value the new value of sidePointer
   */
  void setSidePointer(int value)
  {
    sidePointer = value;
  }

  /**
   * Get the value of sidePointer
   * -1: coordenadas menores
   *  1: coordenadas maiores
   * @return the value of sidePointer
   */
  int getSidePointer()
  {
    return sidePointer;
  }

  void initAttributes();

};

#endif // SURFACE_BM_H
