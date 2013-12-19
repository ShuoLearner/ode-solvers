#include <iostream>
#include <queue>


#ifndef CEXPRKMETHOD
#define CEXPRKMETHOD


#define MAX_STAGE 8
#define EPS       2.220446049250313e-16
#define EPS0      4.940656458412465e-100

//Define function pointer of functions for computing derivative and 
//event function
typedef void (*pEvalR)(const size_t, const double * , const double * ,
		       const size_t *, double * );

//Need static?
// Parameters:
// 1. size_t mDim
// 2. double mT
// 3. double mY
// 4. double Yp
typedef void (*pEvalF)(const size_t *, const double *, 
		       const double *, double * );


struct SRoot
{
  size_t id;
  double t;
};



class CExpRKMethod:
{
  /*============Functions============*/
 public:

  //*********************************//
  //* Constructor and Destructor    *//
  //*********************************//

  /*
   * Default constructor
   */
  CExpRKMethod();

  /*
   * Default destructor
   */
  ~CExpRKMethod();

  //*************************//
  //* Main Process Function *//
  //*************************//
  void integrate();



  //***************************************//
  //* Functions for System Initialization *//
  //***************************************//
  void initialize();



  //***********************************//
  //* Functions for step size control *//
  //***********************************//




  //*****************************//
  //* Function for Root Finder  *//
  //*****************************//


  /*============Attributes============*/

  //*****************************************//
  //* Attributs that should be set by users *//
  //*****************************************//
 public:
  /*
   * mDim, dimension of this problem
   */
  size_t mDim;


  /*
   * mAbsTol, absolute error tolerance 
   * mRelTol. relative error tolerance
   */ 
  double mAbsTol;
  double mRelTol;

  /*
   * mT, current time
   */
  double mT;


  /*
   * mTEnd, terminal time this solver will reach
   */
  double mTEnd;
  
  /*
   * mInitY, a vector to store the initial value of system
   *         Before use, user should clean it
   */ 
  std::vector<double> mInitY;

  /*
   * mDerivFunc, function pointer of function calculating
   *    derivatives
   */
  pEvalF mDerivFunc;


  /*
   * mEventFunc, function pointer of function calculating
   *    event values
   */
  pEvalR mEventFunc;

  
  /*
   * mStatic, a boolin variable 
   * mStatic == 1, write statistic results into a txt file
   * mStatic == 0, do not output statistic results
   */
  bool mStatis;

  /*
   * mODEState, an int varialbe, recording the state of the solver
   * Input:
   *   mODEState == 0, ODE solver is called firstly, need to be initialized
   *   mODEState == 1, ODE solver starts a new integration, without initialization
   *   mODEState == 2, ODE solver continues integration, coming back from event
   *                   has been found
   *
   * Output:
   *   mODEState == -2, some errors happened
   *   mODEState == 3, ODE solver stops at time t < tEnd, indicating having events
   *   mODEState == 4, ODE solver finishes integration at t == tEnd;
   */
  int mODEState;

  
 private:
  //**********************************************//
  //* Variables recording system states          *//
  //**********************************************//
  
  /*
   * mODEStateRecord, recording state of mODEState 
   * before return, in usage of checking correctness
   * of mODEState setting at the initial step
   */
  int mODEStateRecord;
  

  /*
   * mTNew, new time in the next step
   */
  double mTNew;

   /*
   * mY, a double pointer pointing to an array recording 
   *     system values at privous step
   */
  double *mY;

  /*
   * mYNew, a double pointer pointing to an arrya recording
   *        system values at new step
   */
  double *mYNew;

  //***********************************************************//
  //* Some coefficients should be set for process of one step *//
  //* calculation, according to formula                       *//
  //*    ki = f(t_n + c_i*h, y)n + h*(sum a_ij*k_j))          *//
  //*    y_n+1 = y_n + h*(sum b_i*k_i)                        *//
  //*                                                         *//
  //* Absolute local error of one step is calculated as       *//
  //*    absErr = h * (sum e_i*k_i)                           *//
  //***********************************************************//

  /*
   * mP, the order this solver provides
   */
  size_t mP;

  /*
   * mStage, the stage of this method
   */
  size_t mStage;

  /*
   * mA, a double two dimension array, recording 
   *     coefficients a_ij
   */
  double mA[MAX_STAGE][MAX_STAGE];

  /*
   * mB, a double array, recording coefficients b_i
   */
  double mB[MAX_STAGE];

  /*
   * mC, a double array, recording coefficients c_i
   */
  double mC[MAX_STAGE];

  /*
   * mE, a double array, recording coefficients e_i
   *     for absolute error calculation
   */
  double mE[MAX_STAGE];

  /*
   * mK, a double pointer of a two dimension array, recording 
   *     approximated derivatives 
   */
  double **mK;

  //*********************************************************//
  //* Step size h, and step size control related parameters *//
  //*********************************************************//
  
  /*
   * mh, step size
   */
  double mh;

  /*
   * mhMin, minimum step size
   */
  double mhMin;

  /*
   * mhMax, maximum step size
   */
  double mhMax;

  /*
   * mFac, a number as a factor for step size control
   */
  double mFac;

  /*
   * mFacMin, minimum factor value
   */
  double mFacMin;

  /*
   * mFacMax, maximum factor value
   */
  double mFacMax;

  /*
   * mFacMaxRej, maxmum factor value after a rejected step
   */
  double mFacMaxRej;

  //*********************************************************//
  //* Some state records, for usage of ODE solver control   *//
  //*********************************************************//
  
  /*
   * mhFailed, a boolean variable 
   * mhFailed == true, success after a reject step
   * mhFailed == false, previous step is accepted
   */
  bool mhFailed;

  /*
   * mHasEvent, a boolean variable
   * mHasEvent == true,  mEventFunc != NULL
   * mHasEvent == false, mEventFunc == NULL
   */
  bool mHasEvent;

  
  //********************************************//
  //* Integration process statistic variables  *//
  //********************************************//
  
  /*
   * mStepNum, a size_t variable, recording how many steps are executed
   */
  size_t mStepNum;

  /*
   * mAcceptNum, a size_t varialbe, recording how many successful steps
   */
  size_t mAcceptNum;

  /*
   * mfEvalNum, a size_t variable, recording how many times pDerivFunc 
   *            are called
   */
  size_t mfEvalNum;

  //******************************************//
  //* Variables for Root Finding functions   *//
  //******************************************//

 public:
  /*
   * mRootNum, a size_t variable, number of roots
   */
  size_t mRootNum;

 private:
  /*
   * mRootQueue, a queue of struct SRoot, which recording
   * root index and corresponding time, in a time ascending 
   * order
   */
  std::queue<SRoot> mRootQueue;

  /*
   * mState, a double array pointer, recording the state 
   *        of the system
   *        mState[0] = t;
   *        mState[1] ~ mState[dim] = y
   */
  double *mState;

  /*
   * mI, a two dimension double array, for interpolation 
   */
  double mI[MAX_STAGE][MAX_STAGE];

  /*
   * mOrderY, the order of Y interpolation can achieve
   */
  size_t mOrderY;

  /*
   * mOrderYp, the order of Y prime interpolation can achieve
   */
  size_t mOrderYp;

  
  //************************//
  //* Some Other Attributs *//
  //************************//
  double *mZ1, *mZ2, *mZ3;

};


#endif 
