
#ifndef DOUBLE4_H
#define DOUBLE4_H


/**
  * class double4
  * 
  */

class double4
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  double4();

  /**
   * Empty Destructor
   */
  virtual ~double4();

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

  double d1;
  double d2;
  double d3;
  double d4;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of d1
   * @param value the new value of d1
   */
  void setD1(double value)
  {
    d1 = value;
  }

  /**
   * Get the value of d1
   * @return the value of d1
   */
  double getD1()
  {
    return d1;
  }

  /**
   * Set the value of d2
   * @param value the new value of d2
   */
  void setD2(double value)
  {
    d2 = value;
  }

  /**
   * Get the value of d2
   * @return the value of d2
   */
  double getD2()
  {
    return d2;
  }

  /**
   * Set the value of d3
   * @param value the new value of d3
   */
  void setD3(double value)
  {
    d3 = value;
  }

  /**
   * Get the value of d3
   * @return the value of d3
   */
  double getD3()
  {
    return d3;
  }

  /**
   * Set the value of d4
   * @param value the new value of d4
   */
  void setD4(double value)
  {
    d4 = value;
  }

  /**
   * Get the value of d4
   * @return the value of d4
   */
  double getD4()
  {
    return d4;
  }

  void initAttributes();

};

#endif // DOUBLE4_H
