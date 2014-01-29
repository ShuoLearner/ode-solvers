#include <iostream>

#ifndef CEXPRKMETHOD
#define CEXPRKMETHOD


#define MAX_STAGE 8
#define EPS       2.220446049250313e-16
#define EPS0      4.940656458412465e-300

//Define function pointer of functions for computing derivative and 
//event function
// Parameters:
// 1. double mT
// 2. double mY
// 3. size_t mRootNum
// 4. double mRootValue
typedef void (*pEvalR)(const double * , const double * ,
		       const size_t *, double * );

//Need static?
// Parameters:
// 1. double mT
// 2. double mY
// 3. double Yp
typedef void (*pEvalF)(const double *, 
		       const double *, double * );


struct SRoot
{
  size_t id;
  double t;
};


int compare(const void *r1, const void *r2);

class CExpRKMethod
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


  void checkODEState();


  void doOneStep();

  
  double estimateError();


  void advanceStep();




  //***************************************//
  //* Functions for System Initialization *//
  //***************************************//
  void initialize();


  void allocateSpace();


  void setCoeff();


  void setStatRecord();



  //***********************************//
  //* Functions for step size control *//
  //***********************************//

  void setInitialStepSize();


  //*****************************//
  //* Function for Root Finder  *//
  //*****************************//

  void interpolation(const double, double *);


  void copyData();


  void findRoots();


  double rootFindBySecant(const size_t id);


  double rootFindByBisection(const size_t id);


  double rootFindByFalsi(const size_t id);


  void findSlowReaction();


  void calculateRootState();

  //*****************************//
  //* Parameters Check Function *//
  //*****************************//

  void checkParameter();


  //***************************//
  //* Other Helpful Functions *//
  //***************************//
  double infNorm(const size_t &, const double*);


  double dmax(const double&, const double&);


  double dmin(const double&, const double&);


  double dabs(const double&);


  double deps(const double&);


  void clearQueue();

  
  bool queueIsEmpty();

  
  void queuePop();


  void queuePush(const SRoot&);


  /*============Attributes============*/

  //*****************************************//
  //* Attributs that should be set by users *//
  //*****************************************//
 public:
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~Input Parameters~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
   * mDim, dimension of this problem
   */
  size_t mDim;

  /*
   * mRootNum, a size_t variable, number of roots
   */
  size_t mRootNum;

  /*
   * mAbsTol, absolute error tolerance 
   * mRelTol. relative error tolerance
   */ 
  double mAbsTol;
  double mRelTol;

  /*
   * mTEnd, terminal time this solver will reach
   */
  double mTEnd;

  /*
   * mDerivFunc, function pointer of function calculating
   *             derivatives
   */
  pEvalF mDerivFunc;

  /*
   * mEventFunc, function pointer of function calculating
   *             event values
   */
  pEvalR mEventFunc;


  /*
   * mHybrid, a boolean variable
   * mHybrid == false, a regular ODE solver
   * mHybrid == true,  providing inverse interpolation for hybrid method
   */
  bool mHybrid;

  
  /*
   * mStatic, a boolean variable 
   * mStatic == 1, write statistic results into a txt file
   * mStatic == 0, do not output statistic results
   */
  bool mStatis;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~Output Parameters~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
   * mRootId, a size_t variable indicating which root is found
   *          0 <= mRootId <= mRootNum-1
   */
  int mRootId;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~Input and Output Parameters~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*
   * mT, current time
   */
  double mT;

   /*
   * mY, a double pointer pointing to an array recording 
   *     system values at privous step
   */
  double *mY;

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
   * mTCp, copy of mT
   */

  double mTCp;

  /*
   * mTNew, new time in the next step
   */
  double mTNew;

  /*
   * mYNew, a double pointer pointing to an array recording
   *        system values at new step
   */
  double *mYNew;

  /*
   * mYCp, a double pointer pointing to an array recording 
   *        mY
   */
  double *mYCp;

  /*
   * mFinish, a boolean variable, indicating whether integration has 
   *          been finished or not
   */
  bool mFinish;


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
   *     approximated derivatives (mStage*mDim)
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
   * mhNoFailed, a boolean variable 
   * mhNoFailed == true, success after a reject step
   * mhNoFailed == false, previous step is accepted
   */
  bool mhNoFailed;

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
   * mRejectNum, a size_t varialbe, recording how many rejected steps
   */
  size_t mRejectNum;

  /*
   * mfEvalNum, a size_t variable, recording how many times mDerivFunc 
   *            are called
   */
  size_t mfEvalNum;

  /*
   * mrEvalNum, a size_t variable, recording how many times mEventFunc 
   *            are called
   */
  size_t mrEvalNum;


  //******************************************//
  //* Variables for Root Finding functions   *//
  //******************************************//

 private:
  /*
   * mRootQueue, a queue of struct SRoot, which recording
   * root index and corresponding time, in a time ascending 
   * order
   */
  SRoot *mRootQueue;

  size_t mQueueLen;
  size_t mQueueSite;

  /*
   * mI, a two dimension double array, for interpolation 
   */
  double mI[MAX_STAGE][MAX_STAGE];

  /*
   * mRootValueOld, a double array pointer, recording root values
   *                at previous step
   */
  double *mRootValueOld;

  /*
   * mRootValue, a double array pointer, recording root value at
   *             current time mT
   */
  double *mRootValue;
  

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
