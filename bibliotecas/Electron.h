#ifndef ELECTRON_H
#define ELECTRON_H

#include "Particle.h"


/**
  * class Electron
  * 
  */

class Electron : public Particle
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Electron();

  /**
   * Empty Destructor
   */
  virtual ~Electron();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * Colis�o El�stica
   * 
   * Parametros Entrada
   * A, B: parametros de distribui��o angular
   * RNDC: probabilidae de corte
   * 
   * Parametros de Sa�da
   * RMU: deflex�o angular
   */
  void EELa()
  {
  }


  /**
   * Colis�o Inel�stica
   * 
   * Paramentros de Entrada
   * E: energia do eletron (eV)
   * M: Material onde os eletrons se propagam
   * DELTA: corre��o do efeito de densidade de Fermi
   * 
   * Parametros de Sa�da
   * DE: energia perdida (eV)
   * EP: energia do eletron espalhado (eV)
   * CDT: cosseno do angulo polar de espalhamento
   * ES: energia do eletron secund�rio emitido (eV)
   * CDTS: cosseno polar de dire��o do eletron secund�rio
   * IOSC: indice do oscilador que foi ionizado.
   * 
   */
  void EINa()
  {
  }


  /**
   * Efeito Bremsstrahlung
   * 
   * Parametros de Entrada
   * E: energia do eletron
   * M: material no qual o eletron est� percorrendo
   * 
   */
  void EBRa()
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


  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //

};

#endif // ELECTRON_H
