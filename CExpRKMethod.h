#include <iostream>


//Define function pointer of functions for computing derivative and 
//event function
typedef void (*pEvalR)(const size_t, const double * , const double * ,
		       const size_t *, double * );

typedef void (*pEvalF)(const size_t *, const double *, 
		       const double *, double * );


class CExpRKMethod:
{
  /*============Functions============*/
 public:
  CExpRKMethod();

  ~CExpRKMethod();


  /*============Attributes============*/

  //--------ODE Solver Elements--------
 public:
  /*
   * dim, dimension of this problem
   */
  size_t dim;


  /*
   * absTol, absolute error tolerance 
   * relTol. relative error tolerance
   */ 
  double absTol;
  double relTol;

  /*
   * t, current time
   */
  double t;


  /*
   * tEnd, terminal time this solver will reach
   */
  double tEnd;
  
  /*
   * y, a double pointer pointing to an array recording 
   *    system values at privous step
   */
  double *y;

  /*
   * pDerivFunc, function pointer of function calculating
   *    derivatives
   */
  pEvalF pDerivFunc;


  /*
   * pEventFunc, function pointer of function calculating
   *    event values
   */
  pEvalR pEventFunc;

  
  /*
   * static, a boolin variable 
   * static == 1, write statistic results into a txt file
   * static == 0, do not output statistic results
   */
  bool statis;

  
 private:
  //**********************************************//
  //* Variables recording system states          *//
  //**********************************************//
 
  /*
   * tNew, new time in the next step
   */
  double tNew;

  /*
   * yNew, a double pointer pointing to an arrya recording
   *       system values at new step
   */
  double *yNew;

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
   * A, a double pointer of a two dimension array, recording 
   *    coefficients a_ij
   */
  double **A;

  /*
   * B, a double array pointer, recording coefficients b_i
   */
  double *B;

  /*
   * C, a double array pointer, recording coefficients c_i
   */
  double *C;

  /*
   * E, a double array pointer, recording coefficients e_i
   *    for absolute error calculation
   */
  double *E;

  /*
   * K, a double pointer of a two dimension array, recording 
   *    approximated derivatives 
   */
  double **K;

  //*********************************************************//
  //* Step size h, and step size control related parameters *//
  //*********************************************************//
  
  /*
   * p, the order this solver provides
   */
  size_t p;

  /*
   * h, step size
   */
  double h;

  /*
   * hMin, minimum step size
   */
  double hMin;

  /*
   * hMax, maximum step size
   */
  double hMax;

  /*
   * fac, a number as a factor for step size control
   */
  double fac;

  /*
   * facMin, minimum factor value
   */
  double facMin;

  /*
   * facMax, maximum factor value
   */
  double facMax;

  /*
   * facMaxRej, maxmum factor value after a rejected step
   */
  double facMaxRej;

  //*********************************************************//
  //* Some state records, for usage of ODE solver control   *//
  //*********************************************************//
  
  /*
   * hFailed, a boolean variable 
   * hFailed == true, success after a reject step
   * hFailed == false, previous step is accepted
   */
  bool hFailed;
  
  /*
   * odeState, an int varialbe, recording the state of the solver
   * Input:
   *   odeState == 0, ODE solver is called firstly, need to be initialized
   *   odeState == 1, ODE solver is called continuously, without initialization
   *
   * Output:
   *   odeState == -1, some errors happened
   *   odeState == 1, ODE solver stops at time t < tEnd, indicating having events
   *   odeState == 2, ODE solver finishes integration at t == tEnd;
   */
  int odeState;

  /*
   * hasEvent, a boolean variable
   * hasEvent == true,  pEventFunc != NULL
   * hasEvent == false, pEventFunc == NULL
   */
  bool hasEvent;

  
  //********************************************//
  //* Integration process statistic variables  *//
  //********************************************//
  
  /*
   * nStep, a size_t variable, recording how many steps are executed
   */
  size_t nStep;

  /*
   * nAccept, a size_t varialbe, recording how many successful steps
   */
  size_t nAccept;

  /*
   * nfEval, a size_t variable, recording how many times pDerivFunc 
   *         are called
   */
  size_t nfEval;

  //******************************************//
  //* Variables for Root Finding functions   *//
  //******************************************//


};
