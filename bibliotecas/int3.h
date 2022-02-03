

#ifndef INT3_H
#define INT3_H


/**
  * class int3
  * 
  */

class int3
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  int3();

  /**
   * Empty Destructor
   */
  virtual ~int3();

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
  int i3;

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

  /**
   * Set the value of i3
   * @param value the new value of i3
   */
  void setI3(int value)
  {
    i3 = value;
  }

  /**
   * Get the value of i3
   * @return the value of i3
   */
  int getI3()
  {
    return i3;
  }

  void initAttributes();

};

#endif // INT3_H
