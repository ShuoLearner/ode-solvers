#include <iostream>
#include <queue>


#ifndef CEXPRKMETHOD
#define CEXPRKMETHOD

//Define function pointer of functions for computing derivative and 
//event function
typedef void (*pEvalR)(const size_t, const double * , const double * ,
		       const size_t *, double * );

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
   * mY, a double pointer pointing to an array recording 
   *     system values at privous step
   */
  double *mY;

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

  
 private:
  //**********************************************//
  //* Variables recording system states          *//
  //**********************************************//
 
  /*
   * mTNew, new time in the next step
   */
  double mTNew;

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
   * mA, a double pointer of a two dimension array, recording 
   *     coefficients a_ij
   */
  double **mA;

  /*
   * mB, a double array pointer, recording coefficients b_i
   */
  double *mB;

  /*
   * mC, a double array pointer, recording coefficients c_i
   */
  double *mC;

  /*
   * mE, a double array pointer, recording coefficients e_i
   *     for absolute error calculation
   */
  double *mE;

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
   * mODEState, an int varialbe, recording the state of the solver
   * Input:
   *   mODEState == 0, ODE solver is called firstly, need to be initialized
   *   mODEState == 1, ODE solver is called continuously, without initialization
   *
   * Output:
   *   mODEState == -1, some errors happened
   *   mODEState == 1, ODE solver stops at time t < tEnd, indicating having events
   *   mODEState == 2, ODE solver finishes integration at t == tEnd;
   */
  int mODEState;

  /*
   * mHasEvent, a boolean variable
   * mHasEvent == true,  pEventFunc != NULL
   * mHasEvent == false, pEventFunc == NULL
   */
  bool mHasEvent;

  /*
   * mHasInitialized, a boolean variable
   * After user calls integrate(), algorithm should first
   * check whether mHasInitialized == true. If not, 
   * outputs errors
   */
  bool mHasInitialized;

  
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
   * mI, a two dimension double array pointer, for interpolation 
   */
  double **mI;

};


#endif 
