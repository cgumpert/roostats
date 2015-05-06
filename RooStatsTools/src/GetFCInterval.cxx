#include <algorithm>
#include <utility>
#include <set>
#include <iterator>

#include "RooAbsData.h"
#include "RooRealVar.h"
using namespace RooFit;

#include "RooStats/ModelConfig.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/ToyMCSampler.h"
using namespace RooStats;

// custom include(s)
#include "RooStatsTools.h"

namespace CG_Statistics
{
#ifndef CG_EXPERIMENTAL
  HypoTestInverterResult* GetFCInterval(RooAbsData& data,
					ModelConfig& mc,
					double dConf,
					double dLow,
					double dHigh,
					double dStep,
					unsigned int iToys)
#else
  HypoTestInverterResult* GetFCInterval(RooAbsData& data,
					ModelConfig& mc,
					double dConf,
					unsigned int iMaxIterations,
					double dLow,
					double dHigh,
					unsigned int iPoints,
					unsigned int iToys)
#endif // CG_EXPERIMENTAL    
  {
    // get parameter of interest
    RooRealVar* poi = (RooRealVar*)mc.GetParametersOfInterest()->first();

    // check input
    assert(dLow < dHigh);
    assert(dLow < poi->getMax());
    assert(dHigh > poi->getMin());
    assert(iToys > 1/(1 - dConf));
    
    // build profile likelihood test statistics
    ProfileLikelihoodTestStat plts(*mc.GetPdf());

    ModelConfig* bModel = (ModelConfig*)mc.Clone("bModel");
    poi->setVal(0);
    bModel->SetSnapshot(*poi);

    // build the frequentist calculator
    FrequentistCalculator fc(data,*bModel,mc);
    fc.SetToys(iToys,0);
  
    // initialise hypo test inverter
    HypoTestInverter calc(fc,poi,1-dConf);
    if(iVERBOSITY >= eINFO)
      calc.SetVerbose(2);

    // configrue toy MC sample
    ToyMCSampler* toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
    toymcs->SetTestStatistic(&plts);
    if (!mc.GetPdf()->canBeExtended())
      toymcs->SetNEventsPerToy(1);

    // set scan range
    dLow  = std::max(dLow,poi->getMin());
    dHigh = std::min(dHigh,poi->getMax());

#ifndef CG_EXPERIMENTAL
    // run fixed scan
    unsigned int iPoints = (unsigned int)((dHigh - dLow)/dStep + 0.5) + 1;
    calc.SetFixedScan(iPoints,dLow,dHigh);

    // get result
    HypoTestInverterResult* r = calc.GetInterval();
#else
    std::vector<std::pair<double,double> > vScanRanges;
    vScanRanges.push_back(std::make_pair(dLow,dHigh));

    // number of iterations performed
    unsigned int iIterations = 1;
    // list of scanned points with associated p-values
    std::set<std::pair<double,double> > sPoints;
    // obtained results from the scans
    HypoTestInverterResult* r = 0;
    // number of points already scanned
    int iScannedPoints = 0;
    // number of points to use in a scan
    int iCurPoints = iPoints;
    do
    {
      if(iVERBOSITY >= eINFO)
	std::cout << "iteration: " << iIterations << " with " << vScanRanges.size() << " range(s)" << std::endl;
      
      // do fixed scans in region of interests
      for(auto& range : vScanRanges)
      {
	if(iVERBOSITY >= eINFO)
	  std::cout << "scan [" << range.first << " ... " << range.second << "] with " << iCurPoints << " points" << std::endl;
	
	calc.RunFixedScan(iCurPoints,range.first,range.second);
      }

      // get current result
      r = calc.GetInterval();
      assert(r);
      
      // get list of recently scanned points with associated p-values
      for(int iIndex = iScannedPoints; iIndex < r->ArraySize(); ++iIndex)
	sPoints.insert(std::make_pair(r->GetXValue(iIndex),r->CLsplusb(iIndex)));
      assert(!sPoints.empty());

      if(iVERBOSITY >= eDEBUG)
      {
	for(auto it = sPoints.begin(); it != sPoints.end(); ++it)
	  std::cout << (*it).first << ": " << (*it).second << std::endl;
      }
      
      // update number of scanned points
      iScannedPoints = r->ArraySize();

      // clear regions of interest
      vScanRanges.clear();

      // update regions of interest
      // search for adjacent points A and B which have p-values on different sides of the target p-value
      // then add the range [A ... B] to the regions of interest
      // as we already have the p-values at A and B, the actual new range to scan is
      // step = abs(A - B)/ (iPoints - 1) --> [A + step ... B - step] in iPoints - 2 steps (assuming A < B)
      for(auto it = sPoints.begin(); it != std::prev(sPoints.end()); ++it)
      {
	auto next = std::next(it);
	if(((*it).second - (1 - dConf)) * ((*next).second - (1 - dConf)) < 0)
	{
	  dLow = (*it).first + ((*next).first - (*it).first) / (iPoints - 1);
	  dHigh = (*next).first - ((*next).first - (*it).first) / (iPoints - 1);
	  vScanRanges.push_back(std::make_pair(dLow,dHigh));
	}
      }

      
      // delete clone
      delete r;

      // reduce number of scanned points by 2 after first iteration as new range boundaries will already be scanned from previous iteration
      if(iIterations == 1)
	iCurPoints -= 2;
      
      ++iIterations;
    } while(iIterations <= iMaxIterations);

    // get final result
    r = calc.GetInterval();
#endif // CG_EXPERIMENTAL    

    
    // clean up
    delete bModel;

    return r;
  }
}
