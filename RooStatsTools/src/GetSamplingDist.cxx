#include "Math/ProbFunc.h"

#include "RooAbsData.h"
//#include "RooAbsPdf.h"
//#include "RooFitResult.h"
//#include "RooRealVar.h"
//#include "RooLinkedListIter.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/FrequentistCalculator.h"
using namespace RooStats;

namespace CG_Statistics
{
  SamplingDistribution* GetSamplingDist(RooAbsData& data,
					ModelConfig& mc,
					int iToys)
  {
    // need second model
    ModelConfig* bModel = (ModelConfig*)mc.Clone("bModel");

    // initialise frequentist calculator
    FrequentistCalculator fcalc(data,*bModel,mc);
    fcalc.SetToys(iToys,0);

    // use profile likelihood as test statistic
    ProfileLikelihoodTestStat profll(*mc.GetPdf());
    profll.SetOneSidedDiscovery(true);

    ToyMCSampler *toymcs = (ToyMCSampler*)fcalc.GetTestStatSampler();
    toymcs->SetTestStatistic(&profll);

    if (!mc.GetPdf()->canBeExtended())
      toymcs->SetNEventsPerToy(1);
      
    HypoTestResult* pResult = fcalc.GetHypoTest();
    delete bModel;

    return pResult->GetNullDistribution();
  }
}
