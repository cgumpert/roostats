#include "RooAbsData.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
using namespace RooStats;

namespace CG_Statistics
{
  LikelihoodInterval* GetLikelihoodInterval(RooAbsData& data,
					    ModelConfig& mc,
					    double dConf)
  {
    // initialise profile likelihood calculator
    ProfileLikelihoodCalculator plc(data,mc);
    plc.SetConfidenceLevel(dConf);

    return plc.GetInterval();
  }
}
