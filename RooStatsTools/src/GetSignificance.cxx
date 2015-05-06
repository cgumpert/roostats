#include "Math/ProbFunc.h"

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooLinkedListIter.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/FrequentistCalculator.h"
using namespace RooStats;

namespace CG_Statistics
{
  HypoTestResult* GetSignificance(RooAbsData& data,
				  ModelConfig& mc,
				  int iToys)
  {
    if(iToys > 0)
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

      return pResult;
    }
    else
    {
      RooAbsPdf* pPDF = mc.GetPdf();
      assert(pPDF);

      // make sure POI can float
      RooLinkedListIter it = mc.GetParametersOfInterest()->iterator();
      RooRealVar *myarg;
      while((myarg = (RooRealVar *)it.Next()))
	myarg->setConstant(0);
      
      // unconditional fit
      RooFitResult* r = pPDF->fitTo(data,
				    Extended(pPDF->canBeExtended()),
				    Save(),
				    PrintLevel(0));
      r->Print("v");
      double dUncondNll = r->minNll();
      delete r; r = 0;

      // conditional fit

      // load null hypothesis values and set POIs to constant
      mc.LoadSnapshot();
      it = mc.GetParametersOfInterest()->iterator();
      while((myarg = (RooRealVar *)it.Next()))
	myarg->setConstant(1);

      r = pPDF->fitTo(data,
		      Extended(pPDF->canBeExtended()),
		      Save(),
		      PrintLevel(0));
      r->Print("v");
      double dCondNll = r->minNll();
      delete r; r = 0;

      return new HypoTestResult("AsymptoticSignificance",ROOT::Math::normal_cdf_c(sqrt(2 * (dCondNll - dUncondNll))),0);
    }
  }
}
