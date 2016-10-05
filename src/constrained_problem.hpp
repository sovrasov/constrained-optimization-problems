#ifndef CONSTRAINED_PROBLEM_HPP
#define CONSTRAINED_PROBLEM_HPP

#include <vector>

template <class FType>
class ConstrainedProblem
{
protected:

  FType* mPObjective;
  std::vector<FType*> mPConstraints;
  std::vector<double> mConstraintsRHS;

public:
  explicit ConstrainedProblem(FType* objective,
                              const std::vector<FType*>& constraints,
                              const std::vector<double>& constraintsRHS) :
    mPObjective(objective), mPConstraints(constraints), mConstraintsRHS(constraintsRHS)
  {
  }

  double GetFunctionRHS(int fNumber) const
  {
    return mConstraintsRHS[fNumber];
  }

  int GetDimension() const
  {
    return static_cast<int>(mPObjective->GetDimension());
  }

  void GetBounds(double *lb, double* ub)
  {
    mPObjective->GetBounds(lb, ub);
  }

  int GetConstraintsNumber() const
  {
    return static_cast<int>(mPConstraints.size());
  }

  double CalculateFunction(const double* y, int fNumber)
  {
    if(fNumber != static_cast<int>(mPConstraints.size()))
      return mPConstraints[fNumber]->Calculate(y) - mConstraintsRHS[fNumber];
    else
      return mPObjective->Calculate(y);
  }
};

#endif
