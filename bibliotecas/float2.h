
#ifndef FLOAT2_H
#define FLOAT2_H


/**
  * class float2
  * 
  */

class float2
{
public:
  // Constructors/Destructors
  //  


  /**
   * Empty Constructor
   */
  float2();

  /**
   * Empty Destructor
   */
  virtual ~float2();

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

  void initAttributes();

};

#endif // FLOAT2_H
