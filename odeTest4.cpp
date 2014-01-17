#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CExpRKMethod.h"

// Parameters:
// 1. double mT
// 2. double mY
// 3. size_t mRootNum
// 4. double mRootValue
//typedef void (*pEvalR)(const double * , const double * ,
//		       const size_t *, double * );

//Need static?
// Parameters:
// 1. double mT
// 2. double mY
// 3. double Yp
//typedef void (*pEvalF)(const double *, 
//		       const double *, double * );

void deriv(const double *, const double *,
	   double *);

void event(const double *, const double *,
	   const size_t *, double *);


int main()
{
  CExpRKMethod ode45;
  double y[4] = {0, 200, 0, 0};
  double r[10] = {0.235, .1895, .5223, .0023, .6829, .31108, .13829, .22208, .1713, .8002};
  
  ode45.mDim = 4;
  ode45.mY   = y;
  ode45.mRootNum   = 3;
  ode45.mEventFunc = &event;
  ode45.mDerivFunc = &deriv;
  ode45.mHybrid    = true;
  ode45.mODEState  = 0;

  int rId=0;
  y[3] = log(r[rId]);

  ode45.mT    = 0;
  ode45.mTEnd = 100;
  std::cout.precision(18);
  while(1)
    {
      ode45.integrate();
      if(ode45.mODEState == 3)
	{
	  if(ode45.mRootId == -1)
	    {
	      --y[1]; ++y[2];
	      rId  = (rId+1) % 10;
	      y[3] = log(r[rId]);
	      ode45.mODEState = 1;

	      std::cout << "t: " << ode45.mT << "  Slow " << y[2] << std::endl;
	    }
	  else if((ode45.mRootId == 0) || (ode45.mRootId == 1))
	    {
	      ode45.mODEState = 2;
	      //std::cout << "t: " << ode45.mT << "  ID " << ode45.mRootId << std::endl;
	    }
	  else
	    {
	      ode45.mODEState = 1;
	      y[0] = 0; y[1] = 200; y[2] = 0;
	      std::cout << "t: " << ode45.mT << "  ID " << ode45.mRootId << std::endl;
	    }
	}
      else if(ode45.mODEState == 4)
	{
	  std::cout << "Finish Integration!" << std::endl;
	  std::cout << "t: " << ode45.mT << " ";
	  for (int j=0; j<ode45.mDim; ++j)
	    std::cout << y[j] << " ";	    
	  std::cout << std::endl;
	  break;
	}
      else if(ode45.mODEState == -2)
	{
	  std::cout << "Error" << std::endl;
	  break;
	}
    }
  return 0;
}


void deriv(const double *t, const double *y,
	   double *yp)
{
  double y1 = y[0], y2 = y[1], y3 = y[2];
  double k1 = 0.5,  k2 = 1,    k3 = 0.01;
  yp[0] = -k1*y1 + k2*y2;
  yp[1] =  k1*y1 - k2*y2;
  yp[2] =  0;
  yp[3] =  k3*y2;

  return;
}

void event(const double *t, const double *y,
	   const size_t *rn, double *root)
{
  root[0] = y[0] - 10;
  root[1] = y[0] - 10 - 1e-9;
  root[2] = y[2] - 30;
  return;
}
