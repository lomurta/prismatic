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
   * Colisão Elástica
   * 
   * Parametros Entrada
   * A, B: parametros de distribuição angular
   * RNDC: probabilidae de corte
   * 
   * Parametros de Saída
   * RMU: deflexão angular
   */
  void EELa()
  {
  }


  /**
   * Colisão Inelástica
   * 
   * Paramentros de Entrada
   * E: energia do eletron (eV)
   * M: Material onde os eletrons se propagam
   * DELTA: correção do efeito de densidade de Fermi
   * 
   * Parametros de Saída
   * DE: energia perdida (eV)
   * EP: energia do eletron espalhado (eV)
   * CDT: cosseno do angulo polar de espalhamento
   * ES: energia do eletron secundário emitido (eV)
   * CDTS: cosseno polar de direção do eletron secundário
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
   * M: material no qual o eletron está percorrendo
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
