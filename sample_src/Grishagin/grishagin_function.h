#ifndef __GRISHAGIN_FUNCTION_H__
#define __GRISHAGIN_FUNCTION_H__

#define PROPERTY(T, N)    \
  T Get ## N() const;     \
  void Set ## N(T value);

namespace vagrisfunction {

  class GrishaginFunction {

  private:
    int mFunctionNumber;
    unsigned char icnf[45];
    double af[7][7], bf[7][7], cf[7][7], df[7][7];

    double rndm20(unsigned char k[]);
    void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);

  public:
    GrishaginFunction();
    ~GrishaginFunction();

    PROPERTY(int, FunctionNumber);

    double Calculate(const double* y) const;
    double CalculateXDerivative(const double* y) const;
    double CalculateYDerivative(const double* y) const;

    void GetOptimumCoordinates(double* y) const;
    double GetOptimalValue() const;
    void GetDomainBounds(double* lb, double* ub) const;
    int GetDimension() const;
  };
}

#endif
