#ifndef CONSTRAINED_PROBLEM_HPP
#define CONSTRAINED_PROBLEM_HPP

#include <vector>

template <class FType>
  class ConstrainedProblem
{
protected:

public:
  explicit ConstrainedProblem(FType* objective,
                              std::vector<FType*>& constraints,
                              std::vector<double>& parameters)
  {

  }

  double GetFunctionRHS(int fNumber) const
  {
    return 0.;
  }

  int GetConstraintsNumber() const
  {
    return 0;
  }

  double EvaluateFunction(const double* y, int fNumber)
  {
    return 0.;
  }
};

#endif
