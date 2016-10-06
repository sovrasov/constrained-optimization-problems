#include "constrained_problem.hpp"
#include "problem_generator.hpp"
#include "GKLS/gkls_function.hpp"
#include "Grishagin/grishagin_function.hpp"

#include <iostream>

int main(int argc, char** argv)
{
  gkls::GKLSFunction objective, constraint;
  objective.SetFunctionClass(gkls::Hard, 2);
  objective.SetType(gkls::TD);
  constraint.SetFunctionClass(gkls::Hard, 2);
  constraint.SetType(gkls::TD);
  objective.SetFunctionNumber(1);
  constraint.SetFunctionNumber(2);

  ConstrainedProblemGenerator<gkls::GKLSFunction> generator;
  generator.SetObjective(&objective);
  generator.AddConstraint(&constraint, 0.01);

  ConstrainedProblem<gkls::GKLSFunction> problem = generator.GenerateProblem();
  std::cout << "GKLS constraint RHS = " << problem.GetFunctionRHS(0) << "\n";

  vargish::GrishaginFunction gObjective, gConstraint;
  gObjective.SetFunctionNumber(2);
  gConstraint.SetFunctionNumber(1);

  ConstrainedProblemGenerator<vargish::GrishaginFunction> gGenerator;
  gGenerator.SetObjective(&gObjective);
  gGenerator.AddConstraint(&gConstraint, 0.1);

  ConstrainedProblem<vargish::GrishaginFunction> gProblem = gGenerator.GenerateProblem();
  std::cout << "Grishagin constraint RHS = " << gProblem.GetFunctionRHS(0) << "\n";

  return 0;
}
