#ifndef PROBLEM_GENERATOR_HPP
#define PROBLEM_GENERATOR_HPP

#include <vector>
#include <iostream>
#include <cmath>
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

  double EvaluateRHS(FType* function, double delta, double objectiveMin)
  {
    double q;
    //double hmin, double minInGlobal, double Epsilon, int m
    double Epsilon = 0.01; int m = 100;
    int dimension = function->GetDimension();
    double hmin = function->GetOptimalValue();

    double hmax = hmin;
    double d = 0;
    //многомерная решетка, в узлах - испытания,
    int* size = new int[dimension];//кол-во колво точек на размерность
    double* step = new double[dimension];//шаг по каждой размерности
    int sumn = 1;//число испытаний

    double* a = new double[dimension];
    double* b = new double[dimension];
    function->GetDomainBounds(a, b);
    for (unsigned i = 0; i < dimension; i++)
    {
      d = (b[i] - a[i]);
      size[i] = (int)ceil(d / Epsilon) + 1;//(int)ceil(d / Epsilon) + 1;
      step[i] = d / (size[i] - 1);
      sumn *= (size[i]);
    }
    double* f = new double[sumn];//значение функции
    double* y = new double[dimension];
    //#pragma omp parallel for num_threads(parameters.NumThread)
    for (int i = 0; i < sumn; i++)
    {
      double w;
      int z = i;
      //Вычисляе координаты точек испытания
      for (unsigned j = 0; j < dimension; j++)
      {
        w = z % size[j];//определяем номер узла на i-ой размерносте
        y[j] = a[j] + w * step[j];//левая граница + номер узла на i-ой размерносте * шаг на i-ой размерносте
        z = z / size[j];//для вычисления номера узла на следующей оазмерносте
      }
      //проводим испытание
      f[i] = function->Calculate(y);
      if (f[i] > hmax)
        hmax = f[i];
      if (f[i] < hmin)
        hmin = f[i];
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
    {
      for (int j = 0; j < m; j++)
      {
        if ((f[i] >= h1[j]) && (f[i] <= h2[j]))
        {
          p[j] ++;
          break;
        }
      }
    }

    s[0] = p[0];
    for (int i = 1; i < m; i++)
    {
      s[i] = s[i - 1] + p[i];
    }

    double smax = s[m - 1];
    double g = delta * smax;
    double globalMin = objectiveMin;
    for (int i = 0; i < m; i++)
    {
      if (s[i] >= g)
      {
        q = h2[i];
        break;
      }
    }

    //printf("In constraint %d q = %lf\n", constraintNum, q);
    double dm = delta;
    if (dm == 0)
      dm += 0.1;
    dm = dm * (hmax - hmin);

    if (q < globalMin)
    {
      q = globalMin + dm;
      //printf("In constraint %d replace q  = %lf\n", constraintNum, q);
    }

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

    return q;
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
    double objectiveMinVal = mPObjective->GetOptimalValue();
    for(unsigned i = 0; i < mPConstraints.size(); i++)
    {
      if(mNeedTuneParam[i])
        mConstraintsParams[i] = EvaluateRHS(mPConstraints[i], mConstraintsParams[i], objectiveMinVal);
    }

    return ConstrainedProblem<FType>(mPObjective, mPConstraints, mConstraintsParams);
  }
};

#endif
