#include "constrained_problem.hpp"
#include "problem_generator.hpp"
#include "GKLS/gkls_function.hpp"
#include "Grishagin/grishagin_function.hpp"

#include <iostream>

int main(int argc, const char** argv)
{
  gkls::GKLSFunction objective, constraint;
  objective.SetFunctionClass(gkls::Hard, 3);
  objective.SetType(gkls::TD);
  constraint.SetFunctionClass(gkls::Hard, 3);
  constraint.SetType(gkls::TD);
  objective.SetFunctionNumber(1);
  constraint.SetFunctionNumber(2);

  ConstrainedProblemGenerator<gkls::GKLSFunction> generator;
  generator.SetObjective(&objective);
  generator.AddConstraint(&constraint, 0.01);

  double timer = omp_get_wtime();
  ConstrainedProblem<gkls::GKLSFunction> problem = generator.GenerateProblem();
  std::cout << "GKLS constraint RHS = " << problem.GetFunctionRHS(0) << "\n";

  vagrish::GrishaginFunction gObjective, gConstraint;
  gObjective.SetFunctionNumber(2);
  gConstraint.SetFunctionNumber(1);

  ConstrainedProblemGenerator<vagrish::GrishaginFunction> gGenerator;
  gGenerator.SetObjective(&gObjective);
  gGenerator.AddConstraint(&gConstraint, 0.1);

  ConstrainedProblem<vagrish::GrishaginFunction> gProblem = gGenerator.GenerateProblem();
  std::cout << "Grishagin constraint RHS = " << gProblem.GetFunctionRHS(0) << "\n";

  return 0;
}
