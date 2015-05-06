#include <iostream>

#include "RooMsgService.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverterResult.h"
using namespace RooStats;

// custom include(s)
#include "RooStatsTools.h"

namespace CG_Statistics
{
  double EvaluateQMu(RooAbsData& data,RooAbsPdf& pdf,RooRealVar& poi,double dMuTest)
  {
    // calculate q_mu: arxiv:1007.1727v3 equation 14                                                                                                     
    RooFitResult* r = 0;

    // make unconditional fit                            
    poi.setConstant(false);
    r = pdf.fitTo(data,PrintLevel(-1),Save(),Extended(pdf.canBeExtended()));
    const double dUncondNLL = r->minNll();
    delete r; r = 0;
    
    // check for mu > mu^hat
    if(poi.getVal() > dMuTest)
    {
      return 0.0;
    }
    // do conditional fit
    else
    {                        
      poi.setVal(dMuTest);
      poi.setConstant(true);
      r = pdf.fitTo(data,PrintLevel(-1),Save(),Extended(pdf.canBeExtended()));
      const double dCondNLL = r->minNll();
      delete r; r = 0;
      poi.setConstant(false);

      return 2 * (dCondNLL - dUncondNLL);
    }
  }

  double GetSigma(RooAbsData& data,ModelConfig& mc,double dMuTest)
  {
    RooAbsPdf* pPDF = mc.GetPdf();
    assert(pPDF);

    RooRealVar* pPOI = (RooRealVar*)mc.GetParametersOfInterest()->first();
    assert(pPOI);

    // store global observables                                                                                                                  
    const RooArgSet* pGlobalObservables = mc.GetGlobalObservables();  
    RooArgSet* allVars = pPDF->getVariables();
    RooArgSet globObs;
    if(pGlobalObservables)
      pGlobalObservables->snapshot(globObs);

    // build Asimov dataset and set global observables                                                                                           
    RooArgSet globs;
    pPOI->setVal(dMuTest);
    RooAbsData* pAsimovData = AsymptoticCalculator::MakeAsimovData(data,mc,*pPOI,globs);
    *allVars = globs;
  
    // evaluate q_mu_A for a test value of mu                                                                                                    
    double dMuEval = std::min(pPOI->getMax(),dMuTest + pPOI->getError());
    pPOI->setVal(dMuEval);

    // arxiv:1007.1727v3 equation 54                                                                                                             
    double dSigma = (dMuEval - dMuTest) / sqrt(EvaluateQMu(*pAsimovData,*pPDF,*pPOI,dMuEval));

    delete pAsimovData;
    assert(dSigma > 0);

    // reset global observables                                                                                                                  
    *allVars = globObs;
    delete allVars;

    return dSigma;
  }

  HypoTestInverterResult* GetUpperLimit(RooAbsData& data,
					ModelConfig& mc,
					bool bUseCLs,
					double dConf,
					double dLow,
					double dHigh,
					double dPrecision,
					unsigned int iMaxIterations,
					int iToys)
  {
    RooMsgService::instance().setGlobalKillBelow(ERROR);
    AsymptoticCalculator::SetPrintLevel(-1);
    
    // using toys
    if(iToys > 0)
    {
      return 0;
    }
    // using asymptotic formulae
    else
    {
      // get PDF
      RooAbsPdf* pPDF = mc.GetPdf();
      assert(pPDF);
    
      // get parameter of interes
      RooRealVar* pPOI = (RooRealVar*)mc.GetParametersOfInterest()->first();
      assert(pPOI);

      // make sure POI can go negative
      double dOldMin = pPOI->getMin();
      if(dOldMin >= 0)
	pPOI->setMin(-1 * pPOI->getMax());

      // set scan range
      dLow  = std::max(dLow,pPOI->getMin());
      dHigh = std::min(dHigh,pPOI->getMax());

      // run bi-section to find CLs upper limit
      HypoTestInverterResult* result = new HypoTestInverterResult("AsymptoticCLs",*pPOI,dConf);
      result->UseCLs(bUseCLs);
      unsigned int iIteration = 0;
      double dRelDiff = 1;
      double dMuTest, CLsb, CLb, dSigma, sqrt_q_mu;
      do
      {
	++iIteration;
	// set new tested value
	dMuTest = 0.5 * (dLow + dHigh);

	// getr sqrt(q_mu)
	sqrt_q_mu = sqrt(EvaluateQMu(data,*pPDF,*pPOI,dMuTest));

	// get CL(s+b) according to arxiv:1007.1727v3 equation 59
	CLsb = ROOT::Math::normal_cdf_c(sqrt_q_mu,1);

	// apply correction for CL(s)
	if(bUseCLs)
	{
	  // get sigma(mu^prime) at mu^prime = 0
	  dSigma = GetSigma(data,mc,0);
  
	  // get CL(b) according to arxiv:1007.1727v3 equation 57 with mu^prime = 0
	  CLb = ROOT::Math::normal_cdf_c(sqrt_q_mu - dMuTest/dSigma,1);
	  assert(CLb > 0);
	}
	else
	  CLb = 1.0;

	if(CLsb/CLb > 1 - dConf)
	  dLow = dMuTest;
	else
	  dHigh = dMuTest;

	dRelDiff = fabs(1 - dConf - CLsb/CLb) / (1 - dConf);

	if(iVERBOSITY >= eINFO)
	  std::cout << pPOI->GetName() << " = " << dMuTest << ": CL(s+b) = " << CLsb << " and CL(b) = " << CLb << " --> CL(s) = " << CLsb/CLb << std::endl;

	result->Add(dMuTest,HypoTestResult("AsymptoticCLs",CLb,CLsb));
      }
      while((iIteration < iMaxIterations) && (dRelDiff > dPrecision));

      if(iVERBOSITY >= eINFO)
	std::cout << "found limit after " << iIteration << " iterations" << std::endl;
    
      // reset minimum of POI
      pPOI->setMin(dOldMin);

      return result;
    }
  }
}
