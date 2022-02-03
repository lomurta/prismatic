#ifndef MAIN_H
#define MAIN_H

#include <string>
#include "DataTypes.h"
#include "Module.h"
#include "Energy_Spectrum.h"
#include "Body.h"

#include "Particle.h"
#include "Voxel.h"


/**
  * class Main
  * 
  */

//const int NS = 10000;
//const int NB = 5000;
//const int NXG = 250;
//const int MAXMAT = 10;






class Main
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  Main();

  /**
   * Empty Destructor
   */
  virtual ~Main();

  // Static Public attributes
  //  

  // Public attributes
  //  


  // Public attribute accessor methods
  //  


  // Public attribute accessor methods
  //  



  /**
   * L� o arquivo de entrada e inicializa o PENELOPE e o Pengeom
   */
  void pmrdr()
  {
  }


  /*
   * Entrada de dados de materiais e inicializa��o de rotinas de simula��o
   * 
   * Cada material � definido atraves de um arquivo de entrada que � criado pelo
   * programa Material.exe usando informa��oes contidas no banco de dados.
   * 
   * Argumentos de Entrada
   * EMAX: Energia maxima das particulas ( energia cinetica para eletrons e
   * positrons)
   * NMATER: n�mero de materias da geometria
   * PMFILE: matriz de strings de com os nomes dos materiais utilizados na simula��o.
   * O arquivo PMFILE(M) cont�m os dados de intera��o de radia��o para o material M.
   * ( ou seja, a ordem � importante)
   * IWR: unidade de saida
   * INFO: determina a quantidade de informa��es a serem gravadas. 1, 2 ou 3.
   * 
   * Para os calculos prelimires � preciso j� ter a informa��o dos dados de absor�ao
   * EABS e os parametros de controle da simula��o C1, C2, WCC e WCR.
   */
  void peinit()
  {
  }


  /*
   * L� o arquivo de geometria e configura os arrays para rastrear as particulas.
   * 
   * Argumentos de Entrada
   * PARINP: Array contendo paramentros opcionais
   * NPINP: indice do mario componente de PARINP
   * IRD: unidade de entrada do arquivo de defini��o de geometria
   * IWR: unidade de sa�da com arquivo que descreve a geometria real usada na
   * simula��o
   * 
   * Argumentos de sa�da
   * NMATG: numero de materiais diferentes em corpos inteiros.
   * NBOD: Numero de corpos definidos
   * 
   */
  void geomin()
  {
  }


  /**
   * Simula um nova historia de particula e registra a pontua��o de quantidades
   * relevantes
   */
  void shower()
  {
  }


  /**
   * Inicializa a pilha de part�culas secundarias onde seus estados iniciais s�o
   * armazenadas.
   * Deve ser chamada antes de iniciar a simula��o de cada particula primaria
   */
  void cleans()
  {
  }


  /**
   * � chamada antes de iniciar uma nova particula ( primaria ou secundaria)j e
   * tambem quando uma particula cruza sua interface
   * For�a o pr�ximo evento a ser um evento artificial macia.
   */
  void start()
  {
  }


  /**
   * Determina o comprimento DS da part�cula para o proximo evento de intera��o
   * O parametro de entrada dsMax define o comprimento maxmo permitido para eletrons
   * e positrons, para fotons nao tem efeito.
   * 
   * Paramentro de Entrada
   * dsMax: comprimenoto maximo permitido
   * 
   * Paramentro de S�ida
   * DS: comrpimeito do seguimento
   */
  void jump()
  {
  }


  /**
   * Simula o evento de inter��o, calcula novos valores de energia, dire��o de
   * movimento e armazena os estados iniciais das particulas secundarias geradas, se
   * houver.
   * 
   * Parametros de saida
   * de: Energia depositada pela part�cula no material
   * iCol: Tipo de intera��o sofrida pela part�cula
   */
  void knock()
  {
  }


  /**
   * Calcula os novos cossenos diretores da part�cula apos uma determina colis�o com
   * espalhamento em angulos polares e azimutal.
   * 
   * Parametros de Entrada:
   * U, V, W: cossenos diretores iniciais
   * CDT: coseno do angulo de espalhamento polar
   * DF: angulo de espalhamento azimutal
   * 
   * Parametros de sa�da
   * U, V, W: novos cossenos diretores
   */
  void direct()
  {
  }


  /**
   * Armazena o estado inicial de uma nova particula secundaria
   * � necessario informar os paramentros de
   * Energia: E
   * Posicao: X, Y, Z
   * Dire��o: U, V, W
   * Peso da Particula: WGHT
   * Tipo de Particula: KPAR
   * Rotulo da Particula: ILB
   */
  void stores()
  {
  }


  /**
   * Fornece o estado inicial de uma paticula secundaria e a remove da pilha.
   * 
   * Parametro de Saida
   * left: � o numero de particulas secundarias que parmaneceram na pilha na hora da
   * chamada
   */
  void secPar()
  {
  }


  /**
   * Fornece o tempo de execu��o em segundos
   */
  void timer()
  {
  }


  /**
   */
  void rand()
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

  // Titulo
  std::string title;
  // 1: Eletron
  // 2: Foton
  // 3: Positron
  int kPar;
  // 1: MonoEnergetico
  // 2: PoliEnergetico
  int beamType;
  // Energia do Feixo Mono Energetico
  double sEnerg;
  // Espectro de Energia do Feixe
  Energy_Spectrum spectr;
  // Numero de particulas simuladas
  long int nSimSh;
  // Tempo da Simula��o
  long int time;
  // Coordenadas do centro de origem da fonte
  // SX0, SY0, SZ0
  double3 sPosit;
  // Feixe conico/circular
  // Theta, Phi, Alpha
  double3 sCone;
  // Feixe Retangular/quadratico
  // ThetaL, ThetaU, PhiL, PhiU
  double4 sRecta;
  // ISEED1 e ISSE2
  int2 rSeed;
  // XL, XU, XBN
  float3 gridX;
  // YL, YU, YBN
  float3 gridY;
  // ZL, ZU, ZBN
  float3 gridZ;
  // Nome do arquivo Dump
  std::string dumpTo;
  // Periodo de dumping em segundos
  int dumpP;
  // Arquivo de restaura��o DUMP
  std::string resume;
  // Numero m�ximo de materias. Limitado pela capacidade de memoria da GPU
  int maxMat;
  // Tamanho maximo das plihas de eltrons, fotons e positrons = 1000
  // Ao execeder esse valor, o eletro, foton ou positron menos energetico � trocado pela nova particula a ser armazenada.
  int nms;
  // Corpos da simula��o
  // modulos da simula��o

  // Pilha de Photons
  Particle photons;
  // Pilha de Eletrons
  Particle electrons;
  // Pilha de Positrons
  Particle positrons;
  // Vetor de Voxels criados com base em GRIDX, GRIDY, GRIDZ e GRIDBN
  Voxel voxels;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /*
   * Set the value of title
   * Titulo
   * @param value the new value of title
   */
  void setTitle(std::string value)
  {
    title = value;
  }

  /**
   * Get the value of title
   * Titulo
   * @return the value of title
   */
  std::string getTitle()
  {
    return title;
  }

  /**
   * Set the value of kPar
   * 1: Eletron
   * 2: Foton
   * 3: Positron
   * @param value the new value of kPar
   */
  void setKPar(int value)
  {
    kPar = value;
  }

  /**
   * Get the value of kPar
   * 1: Eletron
   * 2: Foton
   * 3: Positron
   * @return the value of kPar
   */
  int getKPar()
  {
    return kPar;
  }

  /**
   * Set the value of beamType
   * 1: MonoEnergetico
   * 2: PoliEnergetico
   * @param value the new value of beamType
   */
  void setBeamType(int value)
  {
    beamType = value;
  }

  /**
   * Get the value of beamType
   * 1: MonoEnergetico
   * 2: PoliEnergetico
   * @return the value of beamType
   */
  int getBeamType()
  {
    return beamType;
  }

  /**
   * Set the value of sEnerg
   * Energia do Feixo Mono Energetico
   * @param value the new value of sEnerg
   */
  void setSEnerg(double value)
  {
    sEnerg = value;
  }

  /**
   * Get the value of sEnerg
   * Energia do Feixo Mono Energetico
   * @return the value of sEnerg
   */
  double getSEnerg()
  {
    return sEnerg;
  }

  /**
   * Set the value of spectr
   * Espectro de Energia do Feixe
   * @param value the new value of spectr
   */
  void setSpectr(Energy_Spectrum value)
  {
    spectr = value;
  }

  /**
   * Get the value of spectr
   * Espectro de Energia do Feixe
   * @return the value of spectr
   */
  Energy_Spectrum getSpectr()
  {
    return spectr;
  }

  /**
   * Set the value of nSimSh
   * Numero de particulas simuladas
   * @param value the new value of nSimSh
   */
  void setNSimSh(long int value)
  {
    nSimSh = value;
  }

  /**
   * Get the value of nSimSh
   * Numero de particulas simuladas
   * @return the value of nSimSh
   */
  long int getNSimSh()
  {
    return nSimSh;
  }

  /**
   * Set the value of time
   * Tempo da Simula��o
   * @param value the new value of time
   */
  void setTime(long int value)
  {
    time = value;
  }

  /**
   * Get the value of time
   * Tempo da Simula��o
   * @return the value of time
   */
  long int getTime()
  {
    return time;
  }

  /**
   * Set the value of sPosit
   * Coordenadas do centro de origem da fonte
   * SX0, SY0, SZ0
   * @param value the new value of sPosit
   */
  void setSPosit(double3 value)
  {
    sPosit = value;
  }

  /**
   * Get the value of sPosit
   * Coordenadas do centro de origem da fonte
   * SX0, SY0, SZ0
   * @return the value of sPosit
   */
  double3 getSPosit()
  {
    return sPosit;
  }

  /**
   * Set the value of sCone
   * Feixe conico/circular
   * Theta, Phi, Alpha
   * @param value the new value of sCone
   */
  void setSCone(double3 value)
  {
    sCone = value;
  }

  /**
   * Get the value of sCone
   * Feixe conico/circular
   * Theta, Phi, Alpha
   * @return the value of sCone
   */
  double3 getSCone()
  {
    return sCone;
  }

  /**
   * Set the value of sRecta
   * Feixe Retangular/quadratico
   * ThetaL, ThetaU, PhiL, PhiU
   * @param value the new value of sRecta
   */
  void setSRecta(double4 value)
  {
    sRecta = value;
  }

  /**
   * Get the value of sRecta
   * Feixe Retangular/quadratico
   * ThetaL, ThetaU, PhiL, PhiU
   * @return the value of sRecta
   */
  double4 getSRecta()
  {
    return sRecta;
  }

  /**
   * Set the value of rSeed
   * ISEED1 e ISSE2
   * @param value the new value of rSeed
   */
  void setRSeed(int2 value)
  {
    rSeed = value;
  }

  /**
   * Get the value of rSeed
   * ISEED1 e ISSE2
   * @return the value of rSeed
   */
  int2 getRSeed()
  {
    return rSeed;
  }

  /**
   * Set the value of gridX
   * XL, XU, XBN
   * @param value the new value of gridX
   */
  void setGridX(float3 value)
  {
    gridX = value;
  }

  /**
   * Get the value of gridX
   * XL, XU, XBN
   * @return the value of gridX
   */
  float3 getGridX()
  {
    return gridX;
  }

  /**
   * Set the value of gridY
   * YL, YU, YBN
   * @param value the new value of gridY
   */
  void setGridY(float3 value)
  {
    gridY = value;
  }

  /**
   * Get the value of gridY
   * YL, YU, YBN
   * @return the value of gridY
   */
  float3 getGridY()
  {
    return gridY;
  }

  /**
   * Set the value of gridZ
   * ZL, ZU, ZBN
   * @param value the new value of gridZ
   */
  void setGridZ(float3 value)
  {
    gridZ = value;
  }

  /**
   * Get the value of gridZ
   * ZL, ZU, ZBN
   * @return the value of gridZ
   */
  float3 getGridZ()
  {
    return gridZ;
  }

  /**
   * Set the value of dumpTo
   * Nome do arquivo Dump
   * @param value the new value of dumpTo
   */
  void setDumpTo(std::string value)
  {
    dumpTo = value;
  }

  /**
   * Get the value of dumpTo
   * Nome do arquivo Dump
   * @return the value of dumpTo
   */
  std::string getDumpTo()
  {
    return dumpTo;
  }

  /**
   * Set the value of dumpP
   * Periodo de dumping em segundos
   * @param value the new value of dumpP
   */
  void setDumpP(int value)
  {
    dumpP = value;
  }

  /**
   * Get the value of dumpP
   * Periodo de dumping em segundos
   * @return the value of dumpP
   */
  int getDumpP()
  {
    return dumpP;
  }

  /**
   * Set the value of resume
   * Arquivo de restaura��o DUMP
   * @param value the new value of resume
   */
  void setResume(std::string value)
  {
    resume = value;
  }

  /**
   * Get the value of resume
   * Arquivo de restaura��o DUMP
   * @return the value of resume
   */
  std::string getResume()
  {
    return resume;
  }

  /**
   * Set the value of maxMat
   * Numero m�ximo de materias. Limitado pela capacidade de memoria da GPU
   * @param value the new value of maxMat
   */
  void setMaxMat(int value)
  {
    maxMat = value;
  }

  /**
   * Get the value of maxMat
   * Numero m�ximo de materias. Limitado pela capacidade de memoria da GPU
   * @return the value of maxMat
   */
  int getMaxMat()
  {
    return maxMat;
  }

  /**
   * Set the value of nms
   * Tamanho maximo das plihas de eltrons, fotons e positrons = 1000
   * Ao execeder esse valor, o eletro, foton ou positron menos energetico � trocado
   * pela nova particula a ser armazenada.
   * @param value the new value of nms
   */
  void setNms(int value)
  {
    nms = value;
  }

  /**
   * Get the value of nms
   * Tamanho maximo das plihas de eltrons, fotons e positrons = 1000
   * Ao execeder esse valor, o eletro, foton ou positron menos energetico � trocado
   * pela nova particula a ser armazenada.
   * @return the value of nms
   */
  int getNms()
  {
    return nms;
  }

  
  /**
   * Set the value of modules
   * modulos da simula��o
   * @param value the new value of modules
   */



  /**
   * Set the value of photons
   * Pilha de Photons
   * @param value the new value of photons
   */
  void setPhotons(Particle value)
  {
    photons = value;
  }

  /**
   * Get the value of photons
   * Pilha de Photons
   * @return the value of photons
   */
  Particle getPhotons()
  {
    return photons;
  }

  /**
   * Set the value of electrons
   * Pilha de Eletrons
   * @param value the new value of electrons
   */
  void setElectrons(Particle value)
  {
    electrons = value;
  }

  /**
   * Get the value of electrons
   * Pilha de Eletrons
   * @return the value of electrons
   */
  Particle getElectrons()
  {
    return electrons;
  }

  /**
   * Set the value of positrons
   * Pilha de Positrons
   * @param value the new value of positrons
   */
  void setPositrons(Particle value)
  {
    positrons = value;
  }

  /**
   * Get the value of positrons
   * Pilha de Positrons
   * @return the value of positrons
   */
  Particle getPositrons()
  {
    return positrons;
  }

  /**
   * Set the value of voxels
   * Vetor de Voxels criados com base em GRIDX, GRIDY, GRIDZ e GRIDBN
   * @param value the new value of voxels
   */
  void setVoxels(Voxel value)
  {
    voxels = value;
  }

  /**
   * Get the value of voxels
   * Vetor de Voxels criados com base em GRIDX, GRIDY, GRIDZ e GRIDBN
   * @return the value of voxels
   */
  Voxel getVoxels()
  {
    return voxels;
  }

  void initAttributes();

};

#endif // MAIN_H
