#include "constrained_problem.hpp"
#include "problem_generator.hpp"
#include "GKLS/gkls_function.hpp"
#include "Grishagin/grishagin_function.hpp"

#include <iostream>

int main(int argc, char** argv)
{
  gklsfunction::GKLSFunction objective, constraint;
  objective.SetFunctionClass(gklsfunction::Hard, 2);
  objective.SetType(gklsfunction::TD);
  constraint.SetFunctionClass(gklsfunction::Hard, 2);
  constraint.SetType(gklsfunction::TD);
  objective.SetFunctionNumber(1);
  constraint.SetFunctionNumber(2);

  ConstrainedProblemGenerator<gklsfunction::GKLSFunction> generator;
  generator.SetObjective(&objective);
  generator.AddConstraint(&constraint, 0.01);

  ConstrainedProblem<gklsfunction::GKLSFunction> problem = generator.GenerateProblem();
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
