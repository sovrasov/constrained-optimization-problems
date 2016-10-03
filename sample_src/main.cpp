#include "constrained_problem.hpp"
#include "problem_generator.hpp"
#include "GKLS/gkls_function.hpp"
#include "Grishagin/grishagin_function.hpp"

#include <iostream>

class GKLSDProxy
{
protected:
  gklsfunction::GKLSFunction* mPFunction;

public:
  GKLSDProxy(gklsfunction::GKLSFunction* function) :
    mPFunction(function)
  {}

  double Calculate(const double* y)
  {
    return mPFunction->EvaluateDFunction(y);
  }

  void GetDomainBounds(double* lb, double* ub)
  {
    mPFunction->GetDomainBounds(lb, ub);
  }

  double GetOptimalValue() const
  {
    return mPFunction->GetGlobalMinimumValue();
  }

  int GetDimension() const
  {
    return mPFunction->GetDimension();
  }
};

int main(int argc, char** argv)
{
  gklsfunction::GKLSFunction objective, constraint;
  objective.SetFunctionClass(gklsfunction::Simple, 2);
  constraint.SetFunctionClass(gklsfunction::Simple, 2);
  objective.SetFunctionNumber(1);
  constraint.SetFunctionNumber(2);

  ConstrainedProblemGenerator<GKLSDProxy> generator;
  GKLSDProxy pObjective(&objective), pConstraint(&constraint);
  generator.SetObjective(&pObjective);
  generator.AddConstraint(&pConstraint, 0.01);

  ConstrainedProblem<GKLSDProxy> problem = generator.GenerateProblem();
  std::cout << "GKLS RHS = " << problem.GetFunctionRHS(0) << "\n";

  vagrisfunction::GrishaginFunction gObjective, gConstraint;
  gObjective.SetFunctionNumber(1);
  gConstraint.SetFunctionNumber(2);

  ConstrainedProblemGenerator<vagrisfunction::GrishaginFunction> gGenerator;
  gGenerator.SetObjective(&gObjective);
  gGenerator.AddConstraint(&gConstraint, 0.9);

  ConstrainedProblem<vagrisfunction::GrishaginFunction> gProblem = gGenerator.GenerateProblem();
  std::cout << "Grishagin RHS = " << gProblem.GetFunctionRHS(0) << "\n";

  return 0;
}
