#ifndef PROBLEM_GENERATOR_HPP
#define PROBLEM_GENERATOR_HPP

#include <vector>
#include "constrained_problem.hpp"

enum GenerateMode { RHS = 0, DELTA = 1 }

template <class FType>
class ConstrainedProblemGenerator
{
protected:

  FType* mPObjective;
  std::vector<FType*> mPConstraints;
  std::vector<double> mConstraintsParams;
  std::vector<bool>   mNeedTuneParam;

  double EvaluateRHS(FType* function, double delta)
  {
    return 0.;
  }

public:

  void AddConstraint(FType* function, double parameter, int mode = DELTA)
  {
    mPConstraints.push_back(function);
    mConstraintsParams.push_back(parameter);
    if(mode == DELTA)
      mNeedTuneParam.push_back(true);
    else
      mNeedTuneParam.push_back(false);
    };
  }

  void SetObjective(FType* function)
  {
    mPObjective = function;
  }

  ConstrainedProblem<FType> GenerateProblem()
  {
    for(size_t i = 0; i < mPConstraints.size(); i++)
    {
      if(mNeedTuneParam[i])
        mConstraintsParams[i] = EvaluateRHS(mPConstraints[i], mConstraintsParams[i]);
    }

    return ConstrainedProblem(mPObjective, mPConstraints, mConstraintsParams);
  }
};

#endif
