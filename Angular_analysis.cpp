#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooStats/SPlot.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include <string>
#include <iostream>

using namespace RooFit;
using namespace RooStats;

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define BINNING 60
#define BINNING2 20
#define WORKSPACE_NAME "B0_PV_M fit"
#define MONTECARLO_PART_REC_BKG "Bu2K1EE"
#define SIG_COLOR 2
#define SIG_MARKSTYL 3
#define SIG_MONT_COLOR 7
#define SIG_MONT_MARKSTYL 23
#define H1_COLOR kBlue
#define H1_FILL 3001
#define H1_PLOT_OPTIONS "BAR"
/*
#define BKG_COLOR 7
#define BKG_MARKSTYL  23
#define BKG2_COLOR 3
#define BKG2_MARKSTYL 34
#define BKG3_COLOR 6
#define BKG3_MARKSTYL 28
*/
void Angular_analysis()
{
//-- LOAD DATA
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("Angular_analysis.cpp","");
  dir.ReplaceAll("/./","/");
  TString workspaceFile = Form("%sDATA/Workspace_B0_PV_M.root",dir.Data());
  TFile *FILE = new TFile(workspaceFile,"read");
  RooWorkspace * ws = (RooWorkspace*)FILE->Get(WORKSPACE_NAME);
// LOAD DATASET
  RooDataSet* sPlotdata = (RooDataSet*) ws->data("dataSetRead");
  RooDataSet* dataSetMonte_sig = (RooDataSet*) ws->data("dataSetMonte_sig");
  RooDataSet* dataSetMonte_partrec_bkg = (RooDataSet*) ws->data("dataSetMonte");
// LOAD VARIABLES
  RooRealVar* cosThetaK = ws->var("cosThetaK");
  RooRealVar* cosThetaL = ws->var("cosThetaL");
  RooRealVar* phi = ws->var("phi");
  RooRealVar* nsig = ws->var("nsig");
  RooRealVar* nbkg_combi = ws->var("nbkg_combi");
  RooRealVar* nbkg_partrec = ws->var("nbkg_partrec");


// -- END OF LOAD DATA
  RooDataSet * SIGNAL = new RooDataSet(sPlotdata->GetName(),sPlotdata->GetTitle(),sPlotdata,*sPlotdata->get(),0,"nsig_sw");
  RooDataSet * BACKGROUND = new RooDataSet(sPlotdata->GetName(),sPlotdata->GetTitle(),sPlotdata,*sPlotdata->get(),0,"nbkg_combi_sw");
  RooDataSet * BACKGROUND2 = new RooDataSet(sPlotdata->GetName(),sPlotdata->GetTitle(),sPlotdata,*sPlotdata->get(),0,"nbkg_partrec_sw");

//Rescale parameters
  Float_t rescale_bkg=1./nbkg_partrec->getValV()*dataSetMonte_partrec_bkg->numEntries();
  Float_t rescale_sig=1./nsig->getValV()*dataSetMonte_sig->numEntries();
  TH1 *temph1;
//SIGNAL cosThetaK
  RooPlot* cosThetaK_sig = cosThetaK->frame(Title("Signal cosThetaK")) ;
  temph1= dataSetMonte_sig->createHistogram("Monte_sig_cosThetaK",*cosThetaK);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  cosThetaK_sig->addTH1(temph1,H1_PLOT_OPTIONS);
  //dataSetMonte_sig->plotOn(cosThetaK_sig,Name("Monte_sig"));
  SIGNAL->plotOn(cosThetaK_sig,Binning(BINNING2),Name("sPlot_sig"),Rescale(rescale_sig),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));
//BACKGROUND cosThetaK
  RooPlot* cosThetaK_partrecBKG = cosThetaK->frame(Title("Partially reconstructed background cosThetaK")) ;
  temph1= dataSetMonte_partrec_bkg->createHistogram("Monte_bkg_cosThetaK",*cosThetaK);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  cosThetaK_partrecBKG->addTH1(temph1,H1_PLOT_OPTIONS);
  BACKGROUND2->plotOn(cosThetaK_partrecBKG,Binning(BINNING2),Name("sPlot_partrec_bkg"),Rescale(rescale_bkg),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));

//SIGNAL cosThetaL
  RooPlot* cosThetaL_sig = cosThetaL->frame(Title("Signal cosThetaL")) ;
  temph1= dataSetMonte_sig->createHistogram("Monte_bkg_cosThetaK",*cosThetaL);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  cosThetaL_sig->addTH1(temph1,H1_PLOT_OPTIONS);
  SIGNAL->plotOn(cosThetaL_sig,Binning(BINNING2),Rescale(rescale_sig),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));

//BACKGROUND cosThetaK
  RooPlot* cosThetaL_partrecBKG = cosThetaL->frame(Title("Partially reconstructed background cosThetaL")) ;
  temph1= dataSetMonte_partrec_bkg->createHistogram("Monte_bkg_cosThetaL",*cosThetaL);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  cosThetaL_partrecBKG->addTH1(temph1,H1_PLOT_OPTIONS);
  BACKGROUND2->plotOn(cosThetaL_partrecBKG,Binning(BINNING2),Rescale(rescale_bkg),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));

//SIGNAL phi
  RooPlot* phi_sig = phi->frame(Title("Signal phi"));
  temph1= dataSetMonte_sig->createHistogram("Monte_sig_phi",*phi);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  phi_sig->addTH1(temph1,H1_PLOT_OPTIONS);
  SIGNAL->plotOn(phi_sig,Binning(BINNING2),Rescale(rescale_sig),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));

//BACKGROUND phi
  RooPlot* phi_partrecBKG = phi->frame(Title("Partially reconstructed background phi"));
  temph1= dataSetMonte_partrec_bkg->createHistogram("Monte_bkg_phi",*phi);
  temph1->SetFillColor(H1_COLOR);
  temph1->SetFillStyle(H1_FILL);
  phi_partrecBKG->addTH1(temph1,H1_PLOT_OPTIONS);
  BACKGROUND2->plotOn(phi_partrecBKG,Binning(BINNING2),Rescale(rescale_bkg),DataError(RooAbsData::SumW2),MarkerStyle(SIG_MARKSTYL),MarkerColor(SIG_COLOR));


  TCanvas* c = new TCanvas("All Splots","All Splots",800,400);
  c->Divide(2,3);


  c->cd(1);
  cosThetaK_sig->SetAxisRange(0,120,"Y");
  cosThetaK_sig->GetYaxis()->SetTitleOffset(1.6);
  cosThetaK_sig->Draw();
//LEGEND SIGNAL
  TLegend *leg = new TLegend(0.74,0.70,.99,.99);
  //DRAW BEFORE LEGEND OR IT WILL NOT SHOW
  leg->SetFillColor(kWhite);
  leg->AddEntry(temph1,"Montecarlo signal/background","f");
  leg->AddEntry("sPlot_sig","sPlot signal/background", "P");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry((TObject*)0, "Number of bins = " STR(BINNING2), "");
  leg->Draw();
  c->cd(2);
  cosThetaK_partrecBKG->SetAxisRange(-10,25,"Y");
  cosThetaK_partrecBKG->Draw();
  cosThetaK_partrecBKG->GetYaxis()->SetTitleOffset(1.6) ;
  cosThetaK_partrecBKG->Draw();

// PLOTS COS THETA L
  c->cd(3);
  cosThetaL_sig->SetAxisRange(0,100,"Y");
  cosThetaL_sig->GetYaxis()->SetTitleOffset(1.6);
  cosThetaL_sig->Draw();
  c->cd(4);
  cosThetaL_partrecBKG->SetAxisRange(0,25,"Y");
  cosThetaL_partrecBKG->GetYaxis()->SetTitleOffset(1.6) ;
  cosThetaL_partrecBKG->Draw();

// PLOTS PHI
  c->cd(5);
  phi_sig->SetAxisRange(0,100,"Y");
  phi_sig->GetYaxis()->SetTitleOffset(1.6);
  phi_sig->Draw();
  c->cd(6);
  phi_partrecBKG->SetAxisRange(0,25,"Y");
  phi_partrecBKG->GetYaxis()->SetTitleOffset(1.6) ;
  phi_partrecBKG->Draw();


}
