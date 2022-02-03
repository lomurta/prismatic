


#ifndef PHOTON_H
#define PHOTON_H

#include "Particle.h"
#include "Particle.h"


/**
  * class Photon
  * 
  */

class Photon : public Particle
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Photon();

  /**
   * Empty Destructor
   */
  virtual ~Photon();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * Espalhamento Compton (Incoerente)
   * 
   * Parametros de Entrada
   * E: energia do foton incidente (eV)
   * M: material onde os fotons se propagam
   * 
   * Parametros de Sa�da
   * DE: perda de energia (eV)
   * EP: energia do foton espalhado (eV)
   * CDT: cosseno do angulo de espalhamento polar
   * ES: energia do foton emitido (eV)
   * CDTS: cosseno polar de dire�o do eletron
   * IZZ: numero atomico do atomo onde ocorreu a dispers�o
   * ISH: camada de eletron atomica que foi ionizada
   * 
   */
  void GCOa()
  {
  }


  /**
   * Espalhamento Rayleigh (Coerente)
   * 
   * Parametros de Entrada
   * E: energia do foton incidente (eV)
   * M: material onde os fotons se propagam
   * 
   * Parametros de Sa�da
   * CDT: cosseno do angulo de espalhamento polar
   * 
   * 
   */
  void GRAa()
  {
  }


  /**
   * Efeito FotoEletrico
   * 
   * Parametros de Saida
   * ES: energia cinetica do foton
   * IZZ: numero atomico do atomo onde ocorreu a absor��o
   * ISH: camada de eletrons at�mica que foi ionizada
   */
  void GPHa()
  {
  }


  /**
   * Produ��o de Pares
   * 
   * Parametros de Saida
   * EE: energia cinetica do eletron
   * CDTE: cosseno de dire��o polar do eletron
   * EP: energia cinetica do positron
   * CDTP: cosseno de dire��o polar do positron
   * IZZ: numero atomico do atomo onde ocorreu a absor��o
   * ISH: camada de eletron atomica que foi ionizada
   */
  void GPPa()
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

#endif // PHOTON_H
