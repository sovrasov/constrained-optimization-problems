#ifndef CONSTRAINED_PROBLEM_HPP
#define CONSTRAINED_PROBLEM_HPP

#include <vector>

template <class FType>
class ConstrainedProblem
{
protected:

  FType* mPObjective;
  std::vector<FType*> mPConstraints;
  std::vector<double> mConstraintsParams;

public:
  explicit ConstrainedProblem(FType* objective,
                              const std::vector<FType*>& constraints,
                              const std::vector<double>& parameters)
  {
    mPConstraints = constraints;
    mPObjective = objective;
    mConstraintsParams = parameters;
  }

  double GetFunctionRHS(int fNumber) const
  {
    return mConstraintsParams[fNumber];
  }

  int GetConstraintsNumber() const
  {
    return static_cast<int>(mPConstraints.size());
  }

  double CalculateFunction(const double* y, int fNumber)
  {
    if(fNumber != static_cast<int>(mPConstraints.size()))
      return mPConstraints[fNumber]->Calculate(y) - mConstraintsParams[fNumber];
    else
      return mPObjective->Calculate(y);
  }
};

#endif
