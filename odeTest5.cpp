#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CExpRKMethod.h"

#define MMOL 6.02214129e20
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
  double y[4] = {1.5*MMOL, 1.0*MMOL, 0.0, -1.0};

  ode45.mDim = 4;
  ode45.mY   = y;
  ode45.mDerivFunc = &deriv;
  ode45.mEventFunc = &event;
  ode45.mRootNum   = 1;
  ode45.mHybrid    = false;
  ode45.mODEState  = 0;

  ode45.mAbsTol = 1e-9;
  ode45.mRelTol = 1e-8;

  ode45.mT    = 0;
  ode45.mTEnd = 300;
  std::cout.precision(15);

  while(1)
    {
      ode45.integrate();
      if(ode45.mODEState == 3)
	{
	  std::cout << "t: " << ode45.mT << " ";
	  for (int j=0; j<ode45.mDim-1; ++j)
	    std::cout << y[j]/MMOL << " ";	    
	  std::cout << std::endl;

	  y[0] = 1.5*MMOL;
	  y[1] = 1.0*MMOL;
	  y[2] = 0;
	  ode45.mODEState = 1;
	}
      else if(ode45.mODEState == 4)
	{
	  std::cout << "Finish Integration!" << std::endl;
	  std::cout << "t: " << ode45.mT << " ";
	  for (int j=0; j<ode45.mDim-1; ++j)
	    std::cout << y[j]/MMOL << " ";	    
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
  double k1 = 2.3,  k2 = 0.8,  k3 = 0.02;
  yp[0] = -k1*y1 + k2*y2;
  yp[1] =  k1*y1 - k2*y2 - k3*y2;
  yp[2] =  k3*y2;
  yp[3] =  0;

  return;
}

//typedef void (*pEvalR)(const double * , const double * ,
//		       const size_t *, double * );
void event(const double *t, const double *y, const size_t *n, double *v)
{
  v[0] = y[2]/MMOL - 1.0;
  return;
}
