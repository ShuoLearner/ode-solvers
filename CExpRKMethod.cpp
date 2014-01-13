#include "CExpRKMethod.h"
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//*********************************//
//* Constructor and Destructor    *//
//*********************************//

/*
 * Default Constructor
 */
CExpRKMethod::CExpRKMethod()
{
  // Default error tolerances
  mAbsTol = 1e-12;
  mRelTol = 1e-10;

  // Default function pointers
  mDerivFunc = NULL;
  mEventFunc = NULL;
  
  // Default system state
  mDim    = 0;
  mT      = 0;
  mTEnd   = 0;
  mTNew   = 0;
  mY      = NULL;
  mYNew   = NULL;
  mK      = NULL;
  mFinish = false;

  // Default coefficients
  mP       = 0;
  mStage   = 0;
  mOrderY  = 0;
  mOrderYp = 0;

  // Default step size related 
  mh         = 0;
  mhMin      = 0;
  mhMax      = 0;
  mFac       = 0;
  mFacMin    = 0;
  mFacMax    = 0;
  mFacMaxRej = 0;

  // Default solver status
  mStatis         = false;
  mhNoFailed      = false;
  mHasEvent       = false;
  mODEState       = 0;
  mODEStateRecord = 0;

  // Default root finder related
  mRootId       = -1;
  mRootNum      = 0;
  mTimeRange    = NULL;
  mRootValue    = NULL;
  mRootValueOld = NULL;

  clearQueue(mRootQueue);

  // Default statistic variables
  mStepNum   = 0;
  mAcceptNum = 0;
  mRejectNum = 0;
  mfEvalNum  = 0;


  // Default tempt variable
  mZ1 = NULL;
  mZ2 = NULL;
  mZ3 = NULL;

  std::cout << "Finish Constructur" << std::endl;
}

/*
 * Default Destructor
 */
CExpRKMethod::~CExpRKMethod()
{
  if (mDerivFunc)
    mDerivFunc = NULL;

  if (mEventFunc)
    mEventFunc = NULL;

  if (mYNew)
    {
      delete [] mYNew;
      mYNew = NULL;
    }

  if(mK)
    {
      for(int i=mStage; i>=0; i--)
	delete [] mK[i];
      
      delete [] mK;
      mK = NULL;
    }

  if(mRootValueOld)
    {
      delete [] mRootValueOld;
      mRootValueOld = NULL;
    }

  if(mRootValue)
    {
      delete [] mRootValue;
      mRootValue = NULL;
    }

  if(mTimeRange)
    {
      delete [] mTimeRange;
      mTimeRange = NULL;
    }


  if(mZ1)
    {
      delete [] mZ1;
      mZ1 = NULL;
    }

  if(mZ2)
    {
      delete [] mZ2;
      mZ2 = NULL;
    }

  if(mZ3)
    {
      delete [] mZ3;
      mZ3 = NULL;
    }

  return;
}


//*************************//
//* Main Process Function *//
//*************************//
void CExpRKMethod::integrate()
{
  //====================//
  // 1 check mODEstate  //
  //====================//
  checkODEState();
  std::cout << "After ODEState Check" << std::endl;
  if(mODEState == 1)//Restart
    {
      mFinish = false;
      allocateSpace();
      setInitialStepSize();
      mDerivFunc(&mT, mY, mK[0]);//record derivative to mK

      //std::cout << "Finish ODEState==1" << std::endl;
    }
  else if (mODEState == 3) // has event
    {
      std::cout << "In State mODEState == 3" << std::endl;
      //If events queue isn't empty, deal with the next event, and return
      if (!mRootQueue.empty())
	{
	  std::cout << "Still has root" << std::endl;
	  calculateRootState();
	  return;
	}
      else
	advanceStep();
    }
  else if(mODEState == -2)//has error
    return;

  //=============//
  // 2 Main Loop //
  //=============//
  while(!mFinish)
    {      
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // (1) Check Whether mT is close to mTEnd //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      if (1.1*mh >= (mTEnd-mT))
	{
	  mh = mTEnd - mT;
	  mFinish = true;
	  std::cout << "mT is close to mTEnd" << std::endl;
	}

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // (2) Set Some Parameters before One Step //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      mhNoFailed = true;      

      //~~~~~~~~~~~~~~~~~~~~~~~//
      // (3) Continue One Step //
      //~~~~~~~~~~~~~~~~~~~~~~~//
      while(true)
	{
	  // (i) Do One Single Step
	  doOneStep();
 
	  // (ii) Update Statistic Record
	  mfEvalNum += mStage;
	  
	  // (iii) Estimate Error
	  double error = estimateError();
	  //std::cout << "Get Error: " << error << std::endl;

	  //(iv) Update step size mh
	  if (error > 1.0) // Step Rejected
	    {
	      mhNoFailed = false;
	      mRejectNum++;
	      mh *= 0.5; // Use half step size h
	      if (mh < mhMin)
		{
		  mODEState = -2;
		  mODEStateRecord = mODEState;
		  std::cout << "Failure at t=" << mT << std::endl;
		  std::cout << "Unable to meet integration tolerances without reducing the step size below the smallest value!" << std::endl;
		  return;
		}
	    }
	  else // Step Accept
	    {

	      mhMin = deps(mTNew) * 16.0;

	      //std::cout << "Step Accepted" << std::endl;
	      mAcceptNum++;
	      
	      double fac = pow(1/error, 1.0/(mP+1.0));
	      if (!mhNoFailed) //previous step is rejected
		mh *= dmin(mFacMaxRej, dmax(mFacMin, mFac*fac));
	      else
		mh *= dmin(mFacMax, dmax(mFacMin, mFac*fac));

	      break;
	    }
	} // end while
      mStepNum++;

      //~~~~~~~~~~~~~~~~~~//
      // (4) Check Events //
      //~~~~~~~~~~~~~~~~~~//
      if (mHasEvent)
	findRoots();

      if (mHybrid)
	findSlowReaction();

      if (!mRootQueue.empty()) //has events
	{
	  //do sorting
	  calculateRootState();
	  return;
	}
	 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // (5) Integration Finished //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~//
      if(mFinish)
	{
	  mODEState = 4;
	  mT = mTEnd;
	  for(int i=0; i<mDim; i++)
	    mY[i] = mYNew[i];

	  return;
	}

      //~~~~~~~~~~~~~~~~~~~~~~//
      // (6) Advance New Step //
      //~~~~~~~~~~~~~~~~~~~~~~//
      advanceStep();
      
      //getchar();
    }

  return;
}


/*
 * checkODEState()
 * Check the state attribute mODEState.
 * If mODEState=0, first call of this method, initialization are required 
 * If mODEState=1, restart this method, method key parameters should be checked
 * If mODEState=2, continue from last step which has an event (mODEState=3).
 *                 First check if events left. If not, continue, else return with next event
 * Else, method has error.
 */
void CExpRKMethod::checkODEState()
{
  if (mODEState == 0) // need initialization
    {
      // initialize coefficients, mY and mK, 
      // in case that dimension is updated 
      // for a new problem
      initialize();
      if (mODEState == -2)
	{
	  mODEStateRecord = -2;
	  return;
	}
    }
  else if (mODEState == 1)
    {
      // In such case, parameters must be changed by
      // users, such as error tolerance and mInitY, 
      // exept the dimension of system.
      checkParameter();
      if (mODEState == -2)
	{
	  mODEStateRecord = -2;
	  return;
	}
    }
  else if(mODEState == 2)
    {
      //Deal with errors 
      if (mODEStateRecord != 3) 
	{
	  if (mODEStateRecord == -2)
	    std::cout << "Errors happen when continuing integration. Parameters should be check!" << std::endl;

	  if (mODEStateRecord == 4)
	    std::cout << "Integration has been finished. Please restart this solver by setting new time span!" << std::endl;

	  mODEState = -2;
	  mODEStateRecord = mODEState;
	  return;
	}
      else
	{
	  mODEState = 3;
	  mODEStateRecord = mODEState;
	  return;
	}
    }
  else
    {
      std::cout << "mODEState should be set as 0, 1 or 2!" << std::endl;
      return;
    }

 return;
}


/*
 * doOneStep()
 * This is a virtual function for RK Method which calculates approximated derivatives,
 * new Y and new T.
 */
void CExpRKMethod::doOneStep()
{
  double t;
  // (1) Calculate mK[1] to mK[mStage-1]
  for(int s=1; s<mStage; s++)
    {
      t = mT + mh*mC[s];//tmp time

      for(int i=0; i<mDim; i++)// tmp Y
	mZ1[i] = mY[i];

      for(int i=0; i<s; i++) //tmp Y + Yp*h
	{
	  double a = mA[s][i] * mh;
	  for (int j=0; j<mDim; j++)
	    mZ1[i] += mK[i][j] * a;
	}

      mDerivFunc(&t, mZ1, mK[s]);
    }


  // (2) New Time, mTNew
  size_t s = mStage-1;
  mTNew = mT + mh;

  if(mFinish)
    {
      mTNew = mTEnd;
      mh = mTEnd - mT;
    }
  
  // (3) New Y, mYNew
  for(int i=0; i<mDim; i++)
    mYNew[i] = mY[i];
  
  //std::cout << "mT = " << mT << ", mh = " << mh << std::endl;
  for(int s=0; s<mStage; s++)
    {
      double b = mB[s] * mh;
      for (int i=0; i<mDim; i++)
	mYNew[i] += b * mK[s][i];
    }

  // (4) New Yp, recording it into last row of mK
  mDerivFunc(&mTNew, mYNew, mK[mStage]);

  return;
}


/*
 * estimateError()
 * Function that calculate error in terms of algorithms in book ""
 * Chapter II, Automatic Step Size Control, pp 167-169.
 */
double CExpRKMethod::estimateError()
{
  // (1) Calculate |ynew - ynew*| in terms of mE
  for (int i=0; i<mDim; i++)
    mZ2[i] = 0;

  for(int s=0; s<mStage+1; s++)
    {
      double e = mE[s] * mh;
      for(int i=0; i<mDim; i++)
	mZ2[i] += e * mK[s][i];
    }


  // (2) Calculate Standard sc=Atol + max(|y|,|ynew|)*Rtol
  for(int i=0; i<mDim; i++)
    mZ3[i] = mAbsTol + dmax(dabs(mY[i]), dabs(mYNew[i]))*mRelTol;

  
  // (3) Calculate Error
  double error = 0, tmp;
  for (int i=0; i<mDim; i++)
    {
      tmp = mZ2[i]/mZ3[i];
      error += tmp*tmp;
    }

  error = sqrt(error/mDim);

  return error;
}



/*
 * advanceStep()
 * Set new mT, mY and mK[0]
 */
void CExpRKMethod::advanceStep()
{
  mT = mTNew;
  for(int i=0; i<mDim; i++)
    mY[i] = mYNew[i];

  for(int i=0; i<mDim; i++)
    mK[0][i] = mK[mStage][i];

  if(mEventFunc)
    {
      for(int i=0; i<mRootNum; i++)
	mRootValue[i] = mRootValueOld[i];
    }

  return;
}


//***************************************//
//* Functions for System Initialization *//
//***************************************//
/*
 * initialize()
 * Initialize statistic record variables, coefficients for RK method, mODEState 
 * and check event function mEventFunc.
 */
void CExpRKMethod::initialize()
{
  checkParameter();
  if (mODEState == -2)
    return;

  setStatRecord();
  setCoeff();

  mODEState = 1;

  if (!mEventFunc) //no event 
    {
      mHasEvent = false;
      mRootNum = 0;
    }
  else
    {
      mHasEvent     = true;
      mRootValueOld = new double[mRootNum];
      mRootValue    = new double[mRootNum];

      //calculate mRootValueOld
      (*mEventFunc)(&mT, mY, &mRootNum, mRootValueOld);

      clearQueue(mRootQueue);

      mTimeRange = new SRange[mRootNum];
    }

  return;
}


void CExpRKMethod::allocateSpace()
{  
  // ----(1)----
  if(mYNew)
    delete [] mYNew;

  mYNew = new double[mDim];

  // ----(2)----
  size_t size = (mDim>mRootNum) ? mDim : mRootNum;
  size = (size>(MAX_STAGE+2)) ? size : (MAX_STAGE+2);
  if (mZ1)
    delete [] mZ1;
  
  mZ1 = new double[size];

  if (mZ2)
    delete [] mZ2;
  
  mZ2 = new double[size];

  if (mZ3)
    delete [] mZ3;
  
  mZ3 = new double[size];

  return;
}


void CExpRKMethod::setCoeff()
{
  mP       = 4;
  mStage   = 6;
  mOrderY  = 4;
  mOrderYp = 3;

  //----Set mA----
  double A[6][6] = {
    {          0,            0,           0,         0,            0,  0},
    {       1./5,            0,           0,         0,            0,  0},
    {      3./40,        9./40,           0,         0,            0,  0},
    {     44./45,      -56./15,       32./9,         0,            0,  0},
    {19372./6561, -25360./2187, 64448./6561, -212./729,            0,  0},
    { 9017./3168,     -355./33, 46732./5247,   49./176, -5103./18656,  0}
  };
  
  for(int r=0; r<mStage; r++)
    {
      for (int c=0; c<mStage; c++)
	mA[r][c] = A[r][c];
    }

  //----Set mC----
  double C[6] = {0, 1./5, 3./10, 4./5, 8./9, 1.};
  for(int c=0; c<mStage; c++)
    mC[c] = C[c];

  //----Set mB----
  double B[6] = {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84};
  for(int c=0; c<mStage; c++)
    mB[c] = B[c];

  //----Set mE----
  double E[7] = {71./57600, 0, -71./16695, 71./1920, -17253./339200, 22./525, -1./40};
    for(int c=0; c<mStage+1; c++)
    mE[c] = E[c];

  //----Set mI----
  double I[7][4] = {
    {1.,  -183./64,     37./12,   -145./128},
    {0,          0,          0,           0},
    {0,  1500./371, -1000./159,   1000./371},
    {0,   -125./32   , 125./12,    -375./64},
    {0, 9477./3392,  -729./106, 25515./6784},
    {0,     -11./7,      11./3,     -55./28},
    {0,       3./2,         -4,        5./2}
  };

  for(int r=0; r<mStage+1; r++)
    {
      for(int c=0; c<mOrderY; c++)
	mI[r][c] = I[r][c];
    }

  //----Set mK----
  mK = new double*[mStage+1];
  for (int r=0; r<mStage+1; r++)
    mK[r] = new double[mDim];


  return;
}


/*
 * setStatRecord()
 * Function is used to set statistc record related variables to be 0
 */
void CExpRKMethod::setStatRecord()
{
  mStepNum   = 0;
  mAcceptNum = 0;
  mRejectNum = 0;
  mfEvalNum  = 0;
}



/*
 * setInitialStepSize()
 * Function is used to set the initial step size mh. Algorithm which is applied 
 * is the one given in Book "Solving Ordinary Differential Equitions I", page
 * 169. Vector norm is the infinity norm picking the element having maximum 
 * absolute value.
 */
void CExpRKMethod::setInitialStepSize()
{
  // (1) First, set parameters, related step size control
  double absT   = dabs(mT);

  mhMin      = 16.0 * deps(mT);
  mhMax      = dabs(mTEnd-mT) / 10;
  mFac       = 0.8;
  mFacMin    = 0.1;
  mFacMax    = 4;
  mFacMaxRej = 1.0;
  
  double d0, d1, d2, h0, h1;
  
  // (2) Calculate h0
  d0 = infNorm(mDim, mY); 
  
  mDerivFunc(&mT, mY, mZ1);//mZ1 is y'(t)
  d1 = infNorm(mDim, mZ1);

  if ((d0<1.0e-5) || (d1<1.0e-5))
    h0 = 1.0e-6;
  else
    h0 = 0.01*(d0/d1);

  // (3) Calculate h1
  for(size_t i=0; i<mDim; i++)
    mZ3[i] = mY[i] + h0*mZ1[i]; // mZ3 is y(t+h0)

  double tCp = mT;
  mT += h0;
  mDerivFunc(&mT, mZ3, mZ2);// mZ2 is y'(t+h0)
  for(size_t i=0; i<mDim; i++)
    mZ3[i] = (mZ1[i]-mZ2[i])/h0;  // mZ3 is y''(t+h0)

  d2 = infNorm(mDim, mZ3);

  double tmp = dmax(d1, d2);
  if (tmp <= 1.0e-15)
    h1 = dmax(1.0e-6, h0*1.0e-3);
  else
    h1 = pow(0.01/tmp, 1.0/(mP+1.0));

  // (4) Calculate h
  mh = dmin(100*h0, h1);

  mT = tCp;
  return;

}


//************************//
// Root Finder Functions *//
//************************//

void CExpRKMethod::interpolation(const double tInterp, double *yInterp)
{
  double tmp = (tInterp-mT)/mh;
  double S[MAX_STAGE];

  S[0] = tmp;
  for(int i=1; i<mOrderY; i++)
    S[i] = S[i-1]*tmp;

  for(int d=0; d<mDim; d++)
    {
      yInterp[d] = mY[d];
      
      for(int s=0; s<mOrderY; s++)
	{
	  tmp = 0;
	  
	  for(int j=0; j<mStage+1; j++)
	    tmp += mI[j][s] * mK[j][d];
	    
	  yInterp[d] += mh * tmp * S[s];
	}
    }
  
  return;
}

void CExpRKMethod::findRoots()
{
  // 1. Calculate New Root Value for mRootValue
  (*mEventFunc)(&mT, mY, &mRootNum, mRootValue);

  // 2. Check Sign Change
  bool hasEvent = false;
  for (int i=0; i<mRootNum; i++)
    {
      if (mRootValue[i]*mRootValueOld[i] <= 0)
	{
	  mTimeRange[i].inRange = true;
	  mTimeRange[i].tLeft   = mT;
	  mTimeRange[i].tRight  = mTNew;
	  mTimeRange[i].vLeft   = mRootValueOld[i];
	  mTimeRange[i].vRight  = mRootValue[i];

	  hasEvent = true;
	}
      else
	mTimeRange[i].inRange = false;
    }

  if (!hasEvent)
    return;

  std::cout << "Find root Here" << std::endl;
  
  // 3. Find Roots
  int maxIter = 100;
  double tol, delta, step;
  double tL, tR, vL, vR, tTry;
  SRoot root;

  double *rArray = mZ2;
  double *yTry = mZ1;

  tol = dmax(deps(mT), deps(mTNew)) * 128;
  tol = dmin(tol, dabs(mTNew-mT));

  for (int r=0; r<mRootNum; r++)
    {
      if (!mTimeRange[r].inRange)
	continue;

      // (1) Set Initial Value
      tL = mTimeRange[r].tLeft;
      tR = mTimeRange[r].tRight;
      vL = mTimeRange[r].vLeft;
      vR = mTimeRange[r].vRight;
      
      // (2) Regula Falsi Iteration 
      for(int i=0; i<maxIter; i++)
	{
	  // Find Root
	  bool findRoot = true;
	  delta = tR - tL;
	  if ((dabs(delta)<tol) || (vL == 0.0))
	    root.t = tL;
	  else if(vR == 0.0)
	    root.t  = tR;
	  else
	    findRoot = false;

	  if(findRoot)
	    {
	      root.id = r;
	      mRootQueue.push(root);
	      break;
	    }
	 
	  // Go Next Step
	  step = -vL/(vR-vL) * delta;
	  step = dmax(0.5*tol, dmin(step, delta-0.5*tol));
	  tTry = tL + step;
	
	  interpolation(tTry, yTry);
	  (*mEventFunc)(&tTry, yTry, &mRootNum, rArray);
	  
	  // Update root r
	  if(rArray[r]*vL >=0) // same sign as vL
	    {
	      mTimeRange[r].tLeft = tTry;
	      mTimeRange[r].vLeft = rArray[r];
	    }
	  else
	    {
	      mTimeRange[r].tRight = tTry;
	      mTimeRange[r].vRight = rArray[r];
	    }

	  // Update Other root information
	  for(int j=r+1; j<mRootNum; j++)
	    {
	      if(mTimeRange[j].inRange && (tTry>mTimeRange[j].tLeft) && (tTry<mTimeRange[j].tRight))
		{
		  if(mTimeRange[j].vLeft*rArray[j]>=0)
		    {
		      mTimeRange[j].tLeft = tTry;
		      mTimeRange[j].vLeft = rArray[j];
		    }
		  else
		    {
		      mTimeRange[j].tRight = tTry;
		      mTimeRange[j].vRight = rArray[j];
		    }
		}
	    }// end update root information

	}// for one root
    }

  return;
}


void CExpRKMethod::findSlowReaction()
{
  // Check slow reactin existance
  if (mY[mDim-1]*mYNew[mDim-1]>0)//now slow reaction
    return;

  double *tArray = mZ1;
  double *yArray = mZ2;
  double *yInter = mZ3;
  size_t cnt = 1;

  // Record t and y
  tArray[0] = mT;
  yArray[0] = mY[mDim-1];
  for (int i=1; i<=mStage; i++)
    {
      if (mC[i]>0 && mC[i]<1)
	{
	  tArray[cnt] = mT + mC[i]*mh;
	  interpolation(tArray[cnt], yInter);
	  yArray[cnt] = yInter[mDim-1];
	  cnt++;
	}
    }
  tArray[cnt] = mTNew;
  yArray[cnt] = mYNew[mDim-1];
  cnt++;

  // check whether yArray[i] is close to 0
  SRoot root;

  double tol = EPS;
  for(int i=0; i<cnt; ++i)
    {
      if(dabs(yArray[i]) <= tol)
	{
	  root.id = mRootNum;
	  root.t  = tArray[i];
	  mRootQueue.push(root);
	  return;
	}
    }

  // Do Inverse Interpolation
  double t = 0, localT;
  for(int i=0; i<cnt; i++)
    {
      localT = 1;
      for(int j=0; j<cnt; j++)
	{
	  if(i != j)
	    localT *= yArray[j] / (yArray[j]-yArray[i]);
	}

      t += localT * tArray[i]; 
    }

  root.id = mRootNum;
  root.t  = t;
  mRootQueue.push(root);
  return;
}



void CExpRKMethod::calculateRootState()
{
  if (mRootQueue.empty())
    {
      std::cout << "No root in mRootQueue! Check Code......" << std::endl;
      mODEState = -2;
      return;
    }

  SRoot root = mRootQueue.front();
  mRootQueue.pop();

  mRootId = root.id;
  mT      = root.t;
  interpolation(root.t, mY);

  mODEState = 3;
  return;
}



//*****************************//
//* Parameters Check Function *//
//*****************************//

void CExpRKMethod::checkParameter()
{
  size_t state = mODEState;

  mODEState = -2;

  if(mDim <= 0) //check dimension
    {
      std::cout << "Dimension of system should be POSITIVE" << std::endl;
      return;
    }

  if((mAbsTol<0) || (mRelTol<0) || (mAbsTol+mRelTol==0))
    {
      std::cout << "Error Tolerances should be nonnegativ and at least one should be positive" << std::endl;
      return;
    }

  if(mTEnd <= mT)
    {
      std::cout << "In this solver, we just support positive integration where stop time should be larger than start time!" << std::endl;
      return;
    }

  if(!mY)
    {
      std::cout << "mY should be set before integration!" << std::endl;
      return;
    }

  if(!mDerivFunc)
    {
      std::cout << "Function that calculates direvatives should be set!" << std::endl;
      return;
    }

  if(mEventFunc && (mRootNum<=0))
    {
      std::cout << "Root Number should be an positive integer!" << std::endl;
      return;
    }

  mODEState = state;
  return;

}


//***************************//
//* Other Helpful Functions *//
//***************************//
double CExpRKMethod::infNorm(const size_t &len, const double *y)
{
  double result, tmp;
  result =(y[0]>=0)?y[0]:-y[0];

  for(size_t i=1; i<len; i++)
    {
      tmp = (y[i]>=0)?y[i]:-y[i];
      if(tmp>result)
	result = tmp;
    }

  return result;
}

double CExpRKMethod::dmax(const double &x1, const double &x2)
{
  return (x1>x2)?x1:x2;
}

double CExpRKMethod::dmin(const double &x1, const double &x2)
{
  return (x1>x2)?x2:x1;
}

double CExpRKMethod::dabs(const double &x)
{
  return (x>0)?x:(-1*x);
}

double CExpRKMethod::deps(const double &x)
{
  return (x==0)? EPS0 : dabs(x)*EPS;
}


void CExpRKMethod::clearQueue(std::queue<SRoot> &qRoot)
{
  std::queue<SRoot> emp;
  std::swap(qRoot, emp);
  return;
}
