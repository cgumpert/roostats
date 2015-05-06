#include "RooAbsData.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/SamplingDistribution.h"
using namespace RooStats;

namespace CG_Statistics
{
  // verbosity level (higher value means more verbose)
  enum VERBOSITY {eSILENT = 0, eERROR = 1, eWARNING = 2, eINFO = 3, eDEBUG = 4};
  extern VERBOSITY iVERBOSITY;
  
  // calculate Feldman-Cousins interval
#ifndef CG_EXPERIMENTAL
  HypoTestInverterResult* GetFCInterval(RooAbsData& data,               // dataset
					ModelConfig& mc,                // model definition
					double dConf = 0.683,           // confidence level
					double dLow = -1e6,             // minimum of interval to scan (will be set to max(dLow,poi->getMin())
					double dHigh = 1e6,             // maximum of interval to scan (will be set to min(dHigh,poi->getMax())
					double dStep = 0.05,            // step size to scan interval (determines number of points to test = (dHigh - dLow)/dStep + 1)
					unsigned int iToys = 10000      // number of toys to calculate CL(s+b) at each point
					);
#else
  HypoTestInverterResult* GetFCInterval(RooAbsData& data,                 // dataset
					ModelConfig& mc,                  // model definition
					double dConf = 0.683,             // confidence level
					unsigned int iMaxIterations = 3,  // maximum number of iterations
					double dLow = -1e6,               // minimum of interval to scan (will be set to max(dLow,poi->getMin())
					double dHigh = 1e6,               // maximum of interval to scan (will be set to min(dHigh,poi->getMax())
					unsigned int iPoints = 11,        // number of points used to sample intervals
					unsigned int iToys = 10000        // number of toys to calculate CL(s+b) at each point
					);
#endif // CG_EXPERIMENTAL  

  // calculate confidence interval based on profiled likelihood function
  LikelihoodInterval* GetLikelihoodInterval(RooAbsData& data,           // dataset
					    ModelConfig& mc,            // model definition
					    double dConf = 0.683        // confidence level
					    );

  // calculate significance
  HypoTestResult* GetSignificance(RooAbsData& data,                     // dataset
				  ModelConfig& mc,                      // model definition
				  int iToys = -1                        // number of toys for calculating significance (-1 = asymptotic formulae)
				  );

  // calculate upper limit
  HypoTestInverterResult* GetUpperLimit(RooAbsData& data,                 // dataset
					ModelConfig& mc,                  // model definition
					bool bUseCLs = true,              // do CL(s) instead of CL(s+b) upper limits
					double dConf = 0.95,              // confidence level 
					double dLow = -1e6,               // minimum of interval to scan (will be set to max(dLow,poi->getMin())
					double dHigh = 1e6,               // maximum of interval to scan (will be set to min(dHigh,poi->getMax())
					double dPrecision = 0.01,         // desired relative precision on calculated limit
					unsigned int iMaxIterations = 20, // maximum number of iterations to find the limit
					int iToys = -1                    // number of toys for calculating significance (-1 = asymptotic formulae)
					);

  // get sampling distributions
  SamplingDistribution* GetSamplingDist(RooAbsData& data,
					ModelConfig& mc,
					int iToys = 1000);
}
