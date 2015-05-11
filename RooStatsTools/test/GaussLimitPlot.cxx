// system include(s)
#include <iostream>
#include "unistd.h"

// ROOT include(s)
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TF1.h"

// RooFit include(s)
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

// RooStats include(s)
#include "RooStats/LikelihoodInterval.h"

// custom include(s)
#include "RooStatsTools.h"

using namespace RooFit;
using namespace RooStats;
using namespace CG_Statistics;

void RunGaussLimits(const double xMin,const double xMax,const unsigned int iPoints,const double conf)
{
  // check input
  assert(conf > 0);
  assert(conf < 1);
  assert(xMin < xMax);
  assert(xMin > -5);
  assert(xMax < 5);
  
  // set up simple gaussian model with model config
  RooWorkspace w("gauss");
  w.factory("Gaussian:gaus(x[0,-10,10],mean[0,-8,8],width[1])");
  
  ModelConfig mcGaus("gauss",&w);
  mcGaus.SetPdf("gaus");
  mcGaus.SetObservables("x");
  mcGaus.SetParametersOfInterest("mean");
  w.import(mcGaus);

  // Feldman-Cousins interval
  TGraph* grFC_up = new TGraph(iPoints + 1);
  TGraph* grFC_down = new TGraph(iPoints + 1);
  // likelihood interval without bound
  TGraph* grLogL_up = new TGraph(iPoints + 1);
  TGraph* grLogL_down = new TGraph(iPoints + 1);
  // likelihood interval with bound
  TGraph* grLogL_bound_up = new TGraph(iPoints + 1);
  TGraph* grLogL_bound_down = new TGraph(iPoints + 1);
  // CLs upper limit
  TGraph* grCLs = new TGraph(iPoints + 1);
  // CL(s+b) upper limit
  TGraph* grCLsb = new TGraph(iPoints + 1);
  
  // helpers
  double xObs;
  RooArgSet rObservables(*w.var("x"));
  RooDataSet* pData = 0;
  RooRealVar* pPOI = (RooRealVar*)(mcGaus.GetParametersOfInterest()->first());
  LikelihoodInterval* pInterval = 0;
  HypoTestInverterResult* pCLs = 0;
  HypoTestInverterResult* pCLsb = 0;
  HypoTestInverterResult* pFCResult = 0;
  for(unsigned int i = 0; i <= iPoints; ++i)
  {
    std::cout << "\rstart " << i+1 << " of " << iPoints+1 << " points";
    std::cout.flush();
    
    // create dummy dataset
    xObs = xMin + i * (xMax - xMin)/iPoints;
    pData = new RooDataSet("data","data",rObservables);
    w.var("x")->setVal(xObs);
    pData->add(rObservables);

    // get likelihood interval without bound
    pInterval = GetLikelihoodInterval(*pData,mcGaus,conf);
    grLogL_up->SetPoint(i,xObs,pInterval->UpperLimit(*pPOI));
    grLogL_down->SetPoint(i,xObs,pInterval->LowerLimit(*pPOI));
    delete pInterval;
    
    // get likelihood interval with bound
    pPOI->setMin(0);
    pInterval = GetLikelihoodInterval(*pData,mcGaus,conf);
    grLogL_bound_up->SetPoint(i,xObs,pInterval->UpperLimit(*pPOI));
    grLogL_bound_down->SetPoint(i,xObs,pInterval->LowerLimit(*pPOI));
    pPOI->setMin(-5);
    delete pInterval;
    
    //get Feldman-Cousin interval
    pPOI->setMin(0);
#ifndef CG_EXPERIMENTAL    
    pFCResult = GetFCInterval(*pData,mcGaus,conf,xObs - 3,xObs + 4,0.1,10000);
#else
    pFCResult = GetFCInterval(*pData,mcGaus,conf,3,-100,100,5,50000);
#endif // CG_EXPERIMENTAL    
    grFC_up->SetPoint(i,xObs,pFCResult->UpperLimit());
    grFC_down->SetPoint(i,xObs,TMath::AreEqualRel(pFCResult->LowerLimit(),pFCResult->UpperLimit(),1e-4) ? 0 : pFCResult->LowerLimit());
    pPOI->setMin(-5);
    delete pFCResult;

    // get CL(s+b) and CLs upper limit
    pCLsb = GetUpperLimit(*pData,mcGaus,false,conf,xObs, xObs + 4);
    pCLs = GetUpperLimit(*pData,mcGaus,true,conf,xObs, xObs + 4);
    grCLsb->SetPoint(i,xObs,pCLsb->UpperLimit());
    grCLs->SetPoint(i,xObs,pCLs->UpperLimit());
    delete pCLsb;
    delete pCLs;
    
    // clean up
    delete pData;
  }
  std::cout << std::endl;

  // beautfiy graphs
  grLogL_up->SetLineWidth(2);
  grLogL_up->SetLineColor(kCyan);
  grLogL_up->GetXaxis()->SetTitle("x_{obs}");
  grLogL_up->GetYaxis()->SetTitle("mean");
  grLogL_up->SetTitle("Standard Gaussian Problem");
  grLogL_up->SetLineWidth(2);
  grLogL_up->SetLineColor(kBlue);
  grLogL_down->SetLineWidth(2);
  grLogL_down->SetLineColor(kBlue);
  grLogL_bound_up->SetLineWidth(2);
  grLogL_bound_up->SetLineColor(kOrange);
  grLogL_bound_down->SetLineWidth(2);
  grLogL_bound_down->SetLineColor(kOrange);
  grFC_up->SetLineWidth(2);
  grFC_up->SetLineColor(kCyan);
  grFC_down->SetLineWidth(2);
  grFC_down->SetLineColor(kCyan);
  grCLs->SetLineWidth(2);
  grCLs->SetLineColor(kRed);
  grCLsb->SetLineWidth(2);
  grCLsb->SetLineColor(kGreen+3);

  // analytical solutions
  TF1* fClassical = new TF1("classical",Form("x + %.3f",ROOT::Math::normal_quantile(conf,1)),xMin,xMax);
  TF1* fCLs = new TF1("CLs",Form("x + ROOT::Math::normal_quantile(1 - %.5f*ROOT::Math::normal_cdf(x,1),1)",1-conf),xMin,xMax);
  TGraph* grFCHigh = 0;
  TGraph* grFCLow = 0;
  if(conf == 0.95)
  {
    grFCHigh = new TGraph(31);
    grFCHigh->SetPoint(0,-3,0.42);
    grFCHigh->SetPoint(1,-2.8,0.45);
    grFCHigh->SetPoint(2,-2.6,0.48);
    grFCHigh->SetPoint(3,-2.4,0.52);
    grFCHigh->SetPoint(4,-2.2,0.56);
    grFCHigh->SetPoint(5,-2.0,0.62);
    grFCHigh->SetPoint(6,-1.8,0.68);
    grFCHigh->SetPoint(7,-1.6,0.76);
    grFCHigh->SetPoint(8,-1.4,0.86);
    grFCHigh->SetPoint(9,-1.2,0.97);
    grFCHigh->SetPoint(10,-1.0,1.10);
    grFCHigh->SetPoint(11,-0.8,1.25);
    grFCHigh->SetPoint(12,-0.6,1.41);
    grFCHigh->SetPoint(13,-0.4,1.58);
    grFCHigh->SetPoint(14,-0.2,1.77);
    grFCHigh->SetPoint(15,0.0,1.96);
    grFCHigh->SetPoint(16,0.2,2.16);
    grFCHigh->SetPoint(17,0.4,2.36);
    grFCHigh->SetPoint(18,0.6,2.56);
    grFCHigh->SetPoint(19,0.8,2.76);
    grFCHigh->SetPoint(20,1.0,2.96);
    grFCHigh->SetPoint(21,1.2,3.16);
    grFCHigh->SetPoint(22,1.4,3.36);
    grFCHigh->SetPoint(23,1.6,3.56);
    grFCHigh->SetPoint(24,1.8,3.76);
    grFCHigh->SetPoint(25,2.0,3.96);
    grFCHigh->SetPoint(26,2.2,4.16);
    grFCHigh->SetPoint(27,2.4,4.36);
    grFCHigh->SetPoint(28,2.6,4.56);
    grFCHigh->SetPoint(29,2.8,4.76);
    grFCHigh->SetPoint(30,3.0,4.96);

    grFCLow = new TGraph(8);
    grFCLow->SetPoint(0,1.6,0.);
    grFCLow->SetPoint(1,1.8,0.16);
    grFCLow->SetPoint(2,2.0,0.35);
    grFCLow->SetPoint(3,2.2,0.53);
    grFCLow->SetPoint(4,2.4,0.69);
    grFCLow->SetPoint(5,2.6,0.84);
    grFCLow->SetPoint(6,2.8,0.99);
    grFCLow->SetPoint(7,3.0,1.14);
  }

  fClassical->SetLineStyle(kDashed);
  fClassical->SetLineColor(kBlack);
  fCLs->SetLineStyle(kDashed);
  fCLs->SetLineColor(kBlack);
  grFCLow->SetLineStyle(kDashed);
  grFCLow->SetLineColor(kBlack);
  grFCHigh->SetLineStyle(kDashed);
  grFCHigh->SetLineColor(kBlack);
  
  TLegend* leg = new TLegend(0.18,0.58,0.48,0.88);
  leg->SetFillColor(kWhite);
  // leg->AddEntry(grQMu,"classical (q_{#mu})","l");
  // leg->AddEntry(grQMuTilde,"#tilde{q}_{#mu}","l");
  leg->AddEntry(grCLs,"CL_{s}","l");
  leg->AddEntry(grCLsb,"CL_{s+b} (q_{#mu})","l");
  leg->AddEntry(grLogL_down,"#Delta log(L) (t_{#mu})","l");
  leg->AddEntry(grLogL_bound_down,"#Delta log(L), #mu #geq 0","l");
  leg->AddEntry(grFC_down,"Feldman-Cousins (#tilde{t}_{#mu})","l");
  
  leg->Draw("same");
  TCanvas* c = new TCanvas("c","Standard Gaussian Problem",800,600);
  c->SetGridx(1);
  c->SetGridy(1);
  c->SetTickx(1);
  c->SetTicky(1);
  grLogL_up->Draw("al same");
  grLogL_down->Draw("l same");
  grLogL_bound_up->Draw("l same");
  grLogL_bound_down->Draw("l same");
  grFC_up->Draw("l same");
  grFC_down->Draw("l same");
  grCLs->Draw("l same");
  grCLsb->Draw("l same");

  fClassical->Draw("c same");
  fCLs->Draw("c same");
  grFCHigh->Draw("c same");
  grFCLow->Draw("c same");
  leg->Draw("same");
  
  c->SaveAs("GaussModelLimits.png");
  
  delete c;
  delete grFC_up;
  delete grFC_down;
  delete grCLs;
  delete grLogL_up;
  delete grLogL_down;
  delete grLogL_bound_up;
  delete grLogL_bound_down;  
}

int main(int argc, char** argv)
{
  RooMsgService::instance().setGlobalKillBelow(ERROR);

  // options to run test
  double xMin          = -3;
  double xMax          = 3;
  unsigned int iPoints = 61;
  double conf          = 0.95;
  VERBOSITY verb       = eSILENT;

  // parse options
  int i;
  while((i = getopt(argc,argv,"l:u:p:c:v:h")) != -1)
  {
    switch(i)
    {
    case 'l':
      xMin = atof(optarg);
      break;
    case 'u':
      xMax = atof(optarg);
      break;
    case 'p':
      iPoints = atoi(optarg);
      break;
    case 'c':
      conf = atof(optarg);
      break;
    case 'v':
      verb = (VERBOSITY)atoi(optarg);
      break;
    case 'h':
    case '?':
      std::cout << "usage: ./GaussLimitPlot -l <LOWER> -u <UPPER> -i <POINTS> -c <CONFIDENCE> -v <VERBOSITY>" << std::endl;
      std::cout << std::endl;
      std::cout << "options:" << std::endl;
      std::cout << "-l LOWER  : lower bound of observed values (default: -3)" << std::endl;
      std::cout << "-u UPPER  : upper bound of observed values (default: 3)" << std::endl;
      std::cout << "-p POINTS : number of points to scan (default: 61)" << std::endl;
      std::cout << "-c CONF   : confidence level 0 < CONF < 1 (default: 0.95)" << std::endl;
      std::cout << "-v VERB   : verbosity level (0 ... silent to 4 .. debug mode) (default: 0)" << std::endl;
      std::cout << "-h        : print this help message" << std::endl;
      return 0;
    default:
      std::cerr << "got unknown option '" << i << "'" << std::endl;
      return 1;
    }
  }

  std::cout << std::endl;
  std::cout << "=============================" << std::endl;
  std::cout << "| running GaussLimiPlot with:" << std::endl;
  std::cout << "|" << std::endl;
  std::cout << "| range = [" << xMin << " ... " << xMax << "]" << std::endl;
  std::cout << "| points = " << iPoints << std::endl;
  std::cout << "| confidenve level = " << conf * 100 << "%" << std::endl;
  std::cout << "=============================" << std::endl;
  std::cout << std::endl;
  
  RunGaussLimits(xMin,xMax,iPoints,conf);

  return 0;
}
