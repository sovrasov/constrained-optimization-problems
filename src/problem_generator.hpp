#ifndef PROBLEM_GENERATOR_HPP
#define PROBLEM_GENERATOR_HPP

#include <vector>
#include <cmath>
#include <algorithm>

#include "constrained_problem.hpp"

enum GenerateMode { RHS = 0, DELTA = 1 };

template <class FType>
class ConstrainedProblemGenerator
{
protected:

  FType* mPObjective;
  std::vector<FType*> mPConstraints;
  std::vector<double> mConstraintsParams;
  std::vector<bool>   mNeedTuneParam;

  double CalculateRHS(FType* function, double delta, const double* objectiveMin)
  {
    double rhs;
    double Epsilon = 0.01;
    int m = 100;
    unsigned dimension = function->GetDimension();
    double hmin = function->GetOptimalValue();

    double hmax = hmin;
    double d = 0;
    //многомерная решетка, в узлах - испытания
    int* size = new int[dimension];//кол-во колво точек на размерность
    double* step = new double[dimension];//шаг по каждой размерности
    int sumn = 1;//число испытаний

    double* a = new double[dimension];
    double* b = new double[dimension];
    function->GetDomainBounds(a, b);
    for (unsigned i = 0; i < dimension; i++)
    {
      d = (b[i] - a[i]);
      size[i] = (int)ceil(d / Epsilon) + 1;
      step[i] = d / (size[i] - 1);
      sumn *= (size[i]);
    }
    double* f = new double[sumn];//значение функции
    double* y = new double[dimension];

    for (int i = 0; i < sumn; i++)
    {
      double w;
      int z = i;
      //Вычисляе координаты точек испытания
      for (unsigned j = 0; j < dimension; j++)
      {
        w = z % size[j];//определяем номер узла на i-ой размерности
        y[j] = a[j] + w * step[j];//левая граница + номер узла на i-ой размерности * шаг на i-ой размерности
        z = z / size[j];//для вычисления номера узла на следующей размерности
      }
      //проводим испытание
      f[i] = function->Calculate(y);
      hmax = std::max(f[i], hmax);
      hmin = std::min(f[i], hmin);
    }

    double* h1 = new double[m];
    double* h2 = new double[m];
    int* p = new int[m];
    int* s = new int[m];

    double deltah = (hmax - hmin) / m;

    for (int i = 0; i < m; i++)
    {
      h1[i] = hmin + i * deltah;
      h2[i] = hmin + (i + 1) * deltah;
      p[i] = 0;
      s[i] = 0;
    }

    for (int i = 0; i < sumn; i++)
      for (int j = 0; j < m; j++)
        if ((f[i] >= h1[j]) && (f[i] <= h2[j]))
        {
          p[j] ++;
          break;
        }

    s[0] = p[0];
    for (int i = 1; i < m; i++)
    {
      s[i] = s[i - 1] + p[i];
    }

    double smax = s[m - 1];
    double g = delta * smax;
    for (int i = 0; i < m; i++)
    {
      if (s[i] >= g)
      {
        rhs = h2[i];
        break;
      }
    }

    double dm = delta;
    if (dm == 0)
      dm += 0.1;
    dm = dm * (hmax - hmin);

    double criticalValue = function->Calculate(objectiveMin);

    if (rhs < criticalValue)
      rhs = criticalValue + dm;

    delete[] size;
    delete[] step;
    delete[] a;
    delete[] b;
    delete[] f;
    delete[] y;

    delete[] h1;
    delete[] h2;
    delete[] p;
    delete[] s;

    return rhs;
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
  }

  void SetObjective(FType* function)
  {
    mPObjective = function;
  }

  ConstrainedProblem<FType> GenerateProblem()
  {
    std::vector<double> objectiveArgMin(mPObjective->GetDimension());
    mPObjective->GetOptimumCoordinates(&objectiveArgMin.front());

    for(unsigned i = 0; i < mPConstraints.size(); i++)
    {
      if(mNeedTuneParam[i])
        mConstraintsParams[i] = CalculateRHS(mPConstraints[i], mConstraintsParams[i], &objectiveArgMin.front());
    }

    return ConstrainedProblem<FType>(mPObjective, mPConstraints, mConstraintsParams);
  }
};

#endif
