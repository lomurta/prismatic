
#ifndef INT2_H
#define INT2_H


/**
  * class int2
  * 
  */

class int2
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  int2();

  /**
   * Empty Destructor
   */
  virtual ~int2();

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

  int i1;
  int i2;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of i1
   * @param value the new value of i1
   */
  void setI1(int value)
  {
    i1 = value;
  }

  /**
   * Get the value of i1
   * @return the value of i1
   */
  int getI1()
  {
    return i1;
  }

  /**
   * Set the value of i2
   * @param value the new value of i2
   */
  void setI2(int value)
  {
    i2 = value;
  }

  /**
   * Get the value of i2
   * @return the value of i2
   */
  int getI2()
  {
    return i2;
  }

  void initAttributes();

};

#endif // INT2_H
