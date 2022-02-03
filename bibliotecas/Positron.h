

#ifndef POSITRON_H
#define POSITRON_H

#include "Particle.h"
#include "Particle.h"


/**
  * class Positron
  * 
  */

class Positron : public Particle
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Positron();

  /**
   * Empty Destructor
   */
  virtual ~Positron();

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
  void PELa()
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
  void PINa()
  {
  }


  /**
   * Aniquila��o de P�sitrons
   * 
   * Parametros de Entrada
   * E: energia do positron
   * M: material onde o positron esta percorrendo
   * 
   * Parametros de Saida
   * E1: energia do foton incidente
   * E2: energia do foton incidente
   * CDT1: cosseno de dire��o poloar do foton 1
   * CDT2: cosseno de dire��o polar do foton 2
   */
  void PANa()
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

#endif // POSITRON_H
