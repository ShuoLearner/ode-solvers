#include <iostream>
#include          "CExpRKMethod.h"

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


int main()
{
  CExpRKMethod ode45;
  double y[3] = {0, 200, 0};
  
  ode45.mDim = 3;
  ode45.mY   = y;
  ode45.mDerivFunc = &deriv;
  ode45.mHybrid    = false;
  ode45.mODEState  = 0;
  ode45.mAbsTol    = 1e-3;
  ode45.mRelTol    = 1e-6;

  for(int i=0; i<10; i++)
    {
      ode45.mT    = i*10;
      ode45.mTEnd = (i+1)*10;
      
      ode45.integrate();
      if(ode45.mODEState == 4)
	{
	  ode45.mODEState = 1;
	  std::cout << "t: " << ode45.mT << " ";
	  for (int j=0; j<ode45.mDim; ++j)
	    std::cout << y[i] << " ";
	    
	  std::cout << std::endl;
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
  yp[1] =  k1*y1 - (k2+k3)*y2;
  yp[2] =  k3*y2;

  return;
}
