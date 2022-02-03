
#ifndef FLOAT5_H
#define FLOAT5_H


/**
  * class float5
  * 
  */

class float5
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  float5();

  /**
   * Empty Destructor
   */
  virtual ~float5();

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

  float f1;
  float f2;
  float f3;
  float f4;
  float f5;

  // Private attribute accessor methods
  //  


  // Private attribute accessor methods
  //  


  /**
   * Set the value of f1
   * @param value the new value of f1
   */
  void setF1(float value)
  {
    f1 = value;
  }

  /**
   * Get the value of f1
   * @return the value of f1
   */
 float getF1()
  {
    return f1;
  }

  /**
   * Set the value of f2
   * @param value the new value of f2
   */
  void setF2(float value)
  {
    f2 = value;
  }

  /**
   * Get the value of f2
   * @return the value of f2
   */
  float getF2()
  {
    return f2;
  }

  /**
   * Set the value of f3
   * @param value the new value of f3
   */
  void setF3(float value)
  {
    f3 = value;
  }

  /**
   * Get the value of f3
   * @return the value of f3
   */
  float getF3()
  {
    return f3;
  }

  /**
   * Set the value of f4
   * @param value the new value of f4
   */
  void setF4(float value)
  {
    f4 = value;
  }

  /**
   * Get the value of f4
   * @return the value of f4
   */
  float getF4()
  {
    return f4;
  }

  /**
   * Set the value of f5
   * @param value the new value of f5
   */
  void setF5(float value)
  {
    f5 = value;
  }

  /**
   * Get the value of f5
   * @return the value of f5
   */
  float getF5()
  {
    return f5;
  }

  void initAttributes();

};

#endif // FLOAT5_H
