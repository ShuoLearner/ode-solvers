#include "CExpRKMethod.h"


//*********************************//
//* Constructor and Destructor    *//
//*********************************//
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
  mA = NULL;
  mB = NULL;
  mC = NULL;
  mE = NULL;
  mI = NULL;

  mP     = 0;
  mStage = 0;

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
  mODEState       = -2;
  mHasEvent       = false;
  mHasInitialized = false;

  // Default root finder related
  mRootNum = 0;
  mRootQueue.resize(0);
  mState = NULL;

  // Default statistic variables
  mStepNum   = 0;
  mAcceptNum = 0;
  mfEvalNum  = 0;
}

voidCExpRKMethod::~CExpRKMethod()
{
  if (mDerivFunc)
    mDerivFunc = NULL;

  if (mEventFunc)
    mEventFunc = NULL;

  if (mY)
    {
      delete mY;
      mY = NULL;
    }

  if (mYNew)
    {
      delete mYNew;
      mYNew = NULL;
    }

  if (mK)
    {
      delete mK;
      mK = NULL;
    }

  if (mA)
    {
      delete mA;
      mA = NULL;
    }

  if (mB)
    {
      delete mB;
      mB = NULL;
    }

  if (mC)
    {
      delete mC;
      mC = NULL;
    }

  if (mE)
    {
      delete mE;
      mE = NULL;
    }

  if (mI)
    {
      delete mI;
      mI = NULL;
    }
}


//*************************//
//* Main Process Function *//
//*************************//
void CExpRKMethod::integrate()
{



}


//***************************************//
//* Functions for System Initialization *//
//***************************************//
void CExpRKMethod::initialize()
{
  


}

void CExpRKMethod::setInitialStepSize()
{



}
