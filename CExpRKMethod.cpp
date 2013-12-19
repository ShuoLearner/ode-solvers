#include "CExpRKMethod.h"


//*********************************//
//* Constructor and Destructor    *//
//*********************************//

/*
 * Default Constructor
 */
void CExpRKMethod::CExpRKMethod()
{
  // Default error tolerances
  mAbsTol = 1e-12;
  mRelTol = 1e-10;

  // Default function pointers
  mDerivFunc = NULL;
  mEventFunc = NULL;
  
  // Default system state
  mDim  = 0;
  mT    = 0;
  mTEnd = 0;
  mTNew = 0;
  mY    = NULL;
  mYNew = NULL;
  mK    = NULL;

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
  mhFailed        = false;
  mHasEvent       = false;
  mODEState       = 0;
  mODEStateRecord = 0;

  // Default root finder related
  mRootNum = 0;
  mRootQueue.resize(0);
  mState = NULL;

  // Default statistic variables
  mStepNum   = 0;
  mAcceptNum = 0;
  mfEvalNum  = 0;
}

/*
 * Default Destructor
 */
voidCExpRKMethod::~CExpRKMethod()
{
  if (mDerivFunc)
    mDerivFunc = NULL;

  if (mEventFunc)
    mEventFunc = NULL;

  if (mY)
    {
      delete [] mY;
      mY = NULL;
    }

  if (mYNew)
    {
      delete [] mYNew;
      mYNew = NULL;
    }


  if(mK)
    {
      for(int i=mDim-1; i>=0; i--)
	delete [] mK[i];
      
      delete [] mK;
      mK = NULL;
    }
  return;
}


//*************************//
//* Main Process Function *//
//*************************//
void CExpRKMethod::integrate()
{
  //====================//
  // 1 check mODEstate =//
  //====================//
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
  if (mODEState == 1)
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

      setInitialY();
      setInitialStepSize();
    }
  else if(mODEState == 2)
    {
      //Deal with errors 
      if (mODESateRecord != 3) 
	{
	  if (mODEStateRecord == -2)
	    std::cout << "Errors happen when continuing integration. Parameters should be check!" << std::endl;

	  if (mODEStateRecord == 4)
	    std::cout << "Integration has been finished. Please restart this solver by setting new time span!" << std::endl;

	  mODEState = -2;
	  mODEStateRecord = mODEState;
	  return;
	}


      //If events queue isn't empty, deal with the next event, and return
      if (!mRootQueue.empty())
	{

	}
      
      //else do a new step.

    }
  else
    {
      std::cout << "mODEState should be set as 0 or 1!" << std::endl;
      return;
    }
}


//***************************************//
//* Functions for System Initialization *//
//***************************************//
void CExpRKMethod::initialize()
{
  checkParameter();
  if (mODEState == -2)
    return;

  setInitY();
  setCoeff();
  setInitialSetpSize();

  mODEState = 1;
  mHasInitialized = true;

  if (!mEventFunc) //no event 
    mHasEvent = false;
  else
    mHasEvent = true;

  return;
}


void CExpRKMethod::setInitialY()
{
  if (mY)
    delete mY;

  mY = new double[mDim];
 
  std::vector<double>::iterator it = mInitY.begin();
  const std::vector<double>::iterator itEnd = mInitY.end();

  for (int i=0; it<itEnd; it++, i++)
      mY[i] = *it;
  
  return;
}


void CExpRKMethod::setCoeff()
{
  mP       = 4;
  mStage   = 7;
  mOrderY  = 4;
  mOrderYp = 3;

  //----Set mA----
  double A[7][7] = {
    {          0,            0,           0,         0,            0,      0, 0},
    {       1./5,            0,           0,         0,            0,      0, 0},
    {      3./40,        9./40,           0,         0,            0,      0, 0},
    {     44./45,      -56./15,       32./9,         0,            0,      0, 0},
    {19372./6561, -25360./2187, 64448./6561, -212./729,            0,      0, 0},
    { 9017./3168,     -355./33, 46732./5247,   49./176, -5103./18656,      0, 0},
    {    35./384,            0,   500./1113,  125./192,  -2187./6784, 11./84, 0}
  };
  
  for(int r=0; r<mStage; r++)
    {
      for (int c=0; c<mStage; c++)
	mA[r][c] = A[r][c];
    }

  //----Set mC----
  double C[7] = {0, 1./5, 3./10, 4./5, 8./9, 1, 1};
  for(int c=0; c<mStage; c++)
    mC[c] = C[c];

  //----Set mB----
  double B[7] = {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84, 0};
  for(int c=0; c<mStage; c++)
    mB[c] = B[c];


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

  for(int r=0; r<mStage; r++)
    {
      for(int c=0; c<mOrderY; c++)
	mI[r][c] = I[r][c];
    }

  //----Set mK----
  mK = new int*[mDim];
  for (int r=0; r<mDim; r++)
    mK[r] = new int[mStage];


  return;
}



void CExpRKMethod::setInitialStepSize()
{



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

  if(mInitY.size()!=mDim)
    {
      std::cout << "Dimension of initial state should be the same as mDim which has been set!" << std::endl;
      return;
    }

  if(!mDerivFunc)
    {
      std::cout << "Function that calculates direvatives should be set!" << std::endl;
      return;
    }

  mODEState = state;
  return;

}
