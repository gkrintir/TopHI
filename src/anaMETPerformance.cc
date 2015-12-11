#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <list>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include <math.h>
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"


#include "UserCode/diall/interface/anaMETPerformance.h"
#include "UserCode/diall/interface/diParticle.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwJetContainer.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/pfParticle.h"


#include "TLorentzVector.h"

#include "TClass.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooVoigtian.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>


using namespace RooFit;
using namespace std;                       


RooRealVar x ("x", "x", -400, 400); // changed the axis range, we only need 800 for zll_pt going to 500 GeV.
RooRealVar g_w ("g_w", "width Gaus", 10., 0., 100., "GeV");	//40
RooRealVar gamma_Z0 ("gamma_Z0_U", "Z0 width", 2.3, 0, 100, "GeV");	//20
RooRealVar v_m ("v_m", "v_m",0,-10.,10.);

RooVoigtian *voigt;
RooFitResult *result;
RooAddPdf *model;

double f;
double efwhm;
//TCanvas *c1 = new TCanvas ("c1", "c1", 800, 800);

anaMETPerformance::anaMETPerformance(const char *name, const char *title) 
:anaBaseTask(name,title),
  fCheckPid(kFALSE),
  fMetType(),      
  fMinPt(0.),      
  fMuonsName(""),  
  fMuons(0x0),     
  fParticles(0x0), 
  fParticlesName(),
  fZs(0x0),
  fZsName(""),
  fh1NMuons(),
  fh2MetCent(),  
  fh2SumEtCent(),         
  fh3PtEtaPhi()       
{

  for(Int_t j = 0; j<10; ++j)
    fh2MetCentPtMin[j] = 0;
  
}

void
anaMETPerformance::ConstructModel(RooDataHist Hist,RooDataHist *bkg_hist, bool BKGSubtract) {

  f=0;
  efwhm=0;

  v_m.setVal(Hist.mean(x) );
  v_m.setRange(Hist.mean(x) - Hist.sigma(x), Hist.mean(x) + Hist.sigma(x));

  voigt =new RooVoigtian ("voigt", "Voigtian", x, v_m, gamma_Z0, g_w);
  
  
  if(BKGSubtract) {
    RooHistPdf *bkg_pdf = new RooHistPdf("bkg_pdf","bkg_pdf",RooArgSet(x),*bkg_hist);
    RooRealVar *lAbkgFrac =new RooRealVar("AbkgFrac","AbkgFrac",0.5,0.,1.);
    RooFormulaVar * sigbkgFrac= new RooFormulaVar("bkgfrac","@0",RooArgSet(*lAbkgFrac));
    model = new RooAddPdf("modelSB","modelSB",*voigt,*bkg_pdf,*sigbkgFrac);
    result = model->fitTo (Hist, RooFit::Minimizer("Minuit2","migrad"),RooFit::Strategy(2), RooFit::SumW2Error (kFALSE), RooFit::Save (kTRUE), RooFit::PrintLevel (-1));	// -1 verbose
                         
  } else {
      result = voigt->fitTo (Hist, RooFit::Minimizer("Minuit2","migrad"),RooFit::Strategy(2), RooFit::SumW2Error (kFALSE), RooFit::Save (kTRUE), RooFit::PrintLevel (-1));	// -1 verbose //
  }

  //if(result->status()!=0) voigt=0;



  //Get the FWHM
  double sigma = g_w.getVal ();
  double gamma = gamma_Z0.getVal ();
  double esigma = g_w.getError ();
  double egamma = gamma_Z0.getError ();

  double Vsg = result->correlation (g_w, gamma_Z0);
  double Vgs = result->correlation (gamma_Z0, g_w);
  double Vss = result->correlation (g_w, g_w);
  double Vgg = result->correlation (gamma_Z0, gamma_Z0);
  cout << "correlacion Vgs " << Vgs << " y correlacion Vsg" << Vsg << endl;
  f = FWHM (sigma, gamma);
  efwhm = FWHMError_fixed (sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg);

  return;

}

//----------------------------------------------------------
void anaMETPerformance::Exec(Option_t * /*option*/)
{

   diParticle *pPart(0);
   TLorentzVector met;
   int count = 0; 
   bool flag_charge = true;
   //printf("anaMETPerformance executing\n");
   //if(!SelectEvent()) return;

   if(!fInitOutput) CreateOutputObjects();

   //Get objects from event

   //Get particles from which MET will be calculated
   if(!fMuons && !fMuonsName.IsNull()) {
     fMuons = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fMuonsName.Data()));
   }
   if(!fMuons) return;
   
   //Make array for Z candidates
   if(!fEventObjects->FindObject(fZsName) && !fZs) {
      fZs = new TClonesArray("diParticle");
      fZs->SetName(fZsName);
      fEventObjects->Add(fZs);
    }
   if(fZs) fZs->Delete();

   Double_t cent = 5.;//fHiEvent->GetCentrality();
   //printf("%f!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", cent);
   Int_t nmuons = fMuons->GetEntriesFast();
   Printf("nmuons: %d",nmuons);
   fh1NMuons->Fill(nmuons);
   if(nmuons<2) return;


   for (int i = 0; i < fMuons->GetEntriesFast(); i++) {
     particleBase *mu1 = static_cast<particleBase*>(fMuons->At(i));
     if(!mu1) {
       Printf("%s ERROR: couldn't get muon",GetName());
       continue;
     }
     if(fCheckPid)
       if(!CheckPid(mu1)) continue;
     count = 0;
     for (int j = i+1; j < fMuons->GetEntriesFast(); j++) {
       particleBase *mu2 = static_cast<particleBase*>(fMuons->At(j));
       if(!mu2) {
         Printf("%s ERROR: couldn't get muon",GetName());
         continue;
       }
       TLorentzVector l1 = mu1->GetLorentzVector();
       TLorentzVector l2 = mu2->GetLorentzVector();
       TLorentzVector dimu = l1 + l2;
       
       //muons should be of opposite sign
       if(mu1->GetCharge()*mu2->GetCharge()<0) {

         if(fCheckPid)
           if(!CheckPid(mu2)) continue;
                   
         //fh3CentPtInvMass->Fill(cent,dimu.Pt(),dimu.M());
         
         //Store Z candidates in event
         if(fZs) {
           pPart = new ((*fZs)[count])
             diParticle(dimu.Pt(),
                        dimu.Eta(),
                        dimu.Phi(),
                        dimu.M(),
                        11);
           pPart->SetCharge(0);
           pPart->AddParticle(mu1);
           pPart->AddParticle(mu2);
           ++count;
         }
       } else {
         //fh3CentPtInvMassSC->Fill(cent,dimu.Pt(),dimu.M());
	 Printf("muon loop prin %d\n", fMuons->GetEntriesFast()-1);
	 flag_charge = false;
       }
             
     }//muon 2 loop
   }//muon 1 loop
   //Get particles from which MET will be calculated
   if(!fParticles && !fParticlesName.IsNull()) {
     //fEventObjects->Print();
     fParticles = dynamic_cast<TClonesArray*>(fEventObjects->FindObject(fParticlesName.Data()));
     if(!fParticles) {
       //check if in jet branch
       lwJetContainer *jetsCont = dynamic_cast<lwJetContainer*>(fEventObjects->FindObject(fParticlesName.Data()));
       if(jetsCont) fParticles = jetsCont->GetJets();
     }
   }
   
   if(!fParticles) {
     Printf("%s: WARNING: Couldn't locate %s branch",GetName(),fParticlesName.Data());
     return;
   }
   

   //Double_t cent = 5.;//fHiEvent->GetCentrality();
   TLorentzVector p4(0.,0.,0.,0.);
   Double_t sumEt = 0.;

   const Int_t nptmin = 10;
   Double_t ptarr[nptmin] {0.,1.,2.,3.,10.,20.,30.,40.,50.,60.};
   TLorentzVector r4[nptmin];
   for(Int_t j = 0; j<nptmin; ++j)
     r4[j].SetPtEtaPhiM(0.,0.,0.,0.);
   for (int i = 0; i < fParticles->GetEntriesFast(); i++) {
     particleBase *p = static_cast<particleBase*>(fParticles->At(i));
     if(!p) {
       Printf("%s ERROR: couldn't get particle",GetName());
       continue;
     }
     TLorentzVector l;
     if(fMetType==kGen || fMetType==kPFRaw) {
     //  if(p->Pt() < fMinPt)
       //  continue;
       l = p->GetLorentzVector();
     }
     else if(fMetType==kVS) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       l.SetPtEtaPhiM(pf->PtVS(),pf->Eta(),pf->Phi(),pf->M());
     }
     else if(fMetType==kPuppi) {
       pfParticle *pf = dynamic_cast<pfParticle*>(p);
       if(!pf) {
         Printf("%s ERROR: couldn't cast particle to pfParticle",GetName());
         return;
       }
       if (pf->GetPuppiWeight()!=0) printf("weight %f\n", pf->GetPuppiWeight());
       l = pf->GetPuppiWeight()*p->GetLorentzVector();
     }
     
     for(Int_t j = 0; j<nptmin; ++j) {
       if(l.Pt()>ptarr[j])
         r4[j]+=l;
     }
    
     if(l.Pt() < fMinPt) continue;
     fh3PtEtaPhi->Fill(l.Pt(),l.Eta(),l.Phi());
     p4+=l;
     sumEt+=l.Et();
 
   }//particle loop
   Printf("muon loop %d\n", fParticles->GetEntriesFast());
	
   met = -p4;
   //printf("!! mpika edp %f %f\n", met.Pt(), p4.Pt());
   fh2MetCent->Fill(cent,met.Pt());
   fh2SumEtCent->Fill(cent,sumEt);

   Int_t nhists = TMath::Min(nptmin,10);
   for(Int_t j = 0; j<nhists; ++j) {
     TLorentzVector met2 = -r4[j];
     fh2MetCentPtMin[j]->Fill(cent,met2.Pt());
   }

   //int error;
   fh1uParaZllPt->Fill(sumEt);
   if (flag_charge){
     //if (pPart->Pt()>20 && pPart->M()>60 && pPart->M()<120)
     //std::pair<double, double> u = compHadronicRecoilProjU(pPart,met, error, count);
     if (pPart->Pt()>48 && pPart->Pt()<60)
       //fh1uParaZllPt->Fill(std::get<1>(u)+pPart->Pt());
       printf("%f", pPart->Pt());
     
   }
   

  /*


  TString samplephys14; TString variablename; TString xvariable; TString tchannel; bool drawchi2=false; bool WantBKGSubtract=false;
  TString variablenamepng=variablename;
  variablenamepng.ReplaceAll("/","over");

  TString DestFolder; // Different folder for background subtracted or not subtracted files
  if(WantBKGSubtract){
    DestFolder="BKG_Subtraction";
  }
  else {
    DestFolder="Not_BKG_Subtraction";
  }

  TH1::SetDefaultSumw2() ;

   
cout << "sample phys " << samplephys14 << endl;
  TString folder = "DY";
  if (samplephys14.Contains ("TT"))
    folder = "TTbar";
  if (samplephys14.Contains ("GJet"))
    folder = "Gamma";
  if(samplephys14.Contains ("QCD"))
    folder = "QCD";
  if(samplephys14.Contains("Data"))
   folder ="Data"; 
   if (samplephys14.Contains("Pseudo"))
   folder="Pseudo";


cout << " folder   " << folder << "  -   " << "Destfolder " << DestFolder << endl;
  
  
  TString titley = "";
  if (variablename == "pfmetx")
    titley = "#sigma(MET_{x}) GeV";
  if (variablename == "pfmety")
    titley = "#sigma(MET_{y}) GeV";




  TCanvas *c1 = new TCanvas ("c1", "c1", 800, 800);

  c1->SetTickx ();
  c1->SetTicky ();
  c1->SetFillColor (kWhite);
  c1->SetFillStyle (0);
  c1->SetRightMargin (0.05);
  c1->SetTopMargin (0.08);




  TFile filephys14 (samplephys14);



  TTree *treephys14 = (TTree *) filephys14.Get ("METtree");




  std::vector < TH1F * >resolution;
  resolution.clear ();




  int limitdown (0), limitup (0);
  TString strlimitup = "0";
  TString strlimitdown = "0";

  int tempsizexarray = 0;

  if (xvariable == "nVert")
    tempsizexarray = 6;
  if (xvariable == "met_sumEt")
    tempsizexarray = 6;
  if (xvariable == "zll_pt")
    tempsizexarray = 10; //6


  const int sizexarray=tempsizexarray;
  TH1F *histonvertex=new TH1F("histonvertex","histonvertex",50,0,50); // added later
  TH1F *histozll_pt=new TH1F("histozll_pt","histozll_pt",100,0,1200); // added later
  TH1F *histomet_sumEt=new TH1F("histomet_sumEt","histomet_sumEt",100,0,4); // added later
  TH1F *histomet_uPara_zllzll_pt=new TH1F("histomet_uPara_zllzll_pt","histomet_uPara_zllzll_pt",50,-300,300); // added later
  TH1F *histomet_uPerp_zll=new TH1F("histomet_uPerp_zll","histomet_uPerp_zll",50,-300,300); // added later
  //TH1F *met_uPara_zllresponse1=new TH1F("met_uPara_zllresponse1","met_uPara_zll",50,-300,300); // added later
  //TH1F *zll_ptresponse1=new TH1F("zll_ptresponse1","zll_pt",100,0,100); // added later
  TH1F *histoscale=new TH1F("histoscale","histoscale",50,-50,50); // added later
  
  
  TString dileptonch="";
  if (tchannel=="MuMu")  dileptonch="1";
  if (tchannel=="EE" || tchannel=="Gamma") dileptonch="0";
  
  
  // Plot inclusive distributions of the main variables
  TString condition="(xsec)*(puWeight)";//"(weighttotal)*(channel=="+dileptonch +")";
  if (samplephys14.Contains("Data")) condition="";
  //cout << condition.Data() << endl;
  if (!samplephys14.Contains("raw") && !samplephys14.Contains("jes_")){ treephys14->Draw ("nVert >> histonvertex", condition.Data());
  treephys14->Draw ("zll_pt >> histozll_pt", condition.Data());
  treephys14->Draw ("met_sumEt/1000 >> histomet_sumEt", condition.Data());
  treephys14->Draw("met_uPara_zll/zll_pt >> histoscale" ) ;
  treephys14->Draw ("met_uPara_zll+zll_pt >> histomet_uPara_zllzll_pt", condition.Data());
  treephys14->Draw ("met_uPerp_zll >> histomet_uPerp_zll",condition.Data());//condition.Data());//, condition.Data());
   }
          
  histonvertex->Draw();
  histonvertex->GetXaxis ()->SetTitle ("Number of Vertices");
  if (!samplephys14.Contains("jes_"))  c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/nVert_inclusive_.png");
  c1->SetLogy();
  histozll_pt->Draw();
  histozll_pt->GetXaxis ()->SetTitle ("zll_pt [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/"+ folder + "/"+tchannel +"/zll_pt_inclusive_.png");
  histomet_sumEt->Draw();
  histomet_sumEt->GetXaxis ()->SetTitle("sumE_{T} [TeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_sumEt_inclusive_.png");
  histomet_uPara_zllzll_pt->Draw();
  histomet_uPara_zllzll_pt->GetXaxis ()->SetTitle ("u_{||}+zll_pt [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/"+ folder + "/"+tchannel +"/met_uPara_zllzll_pt_inclusive_.png");
  histomet_uPerp_zll->Draw();
  histomet_uPerp_zll->GetXaxis ()->SetTitle ("u_{#perp}   [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPerp_zll_inclusive_.png");
  c1->SetLogy(0);
  histoscale->Draw();
  histoscale->GetXaxis ()->SetTitle ("u_{#perp}   [GeV]");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPara_zlloverzll_pt.png");


  Double_t tgraphx[sizexarray], tgraphy[sizexarray], etgraphy[sizexarray],
    etgraphx[sizexarray], tgraphxchi2[sizexarray], tgraphychi2[sizexarray];
    
    



  for (int index = 0; index < sizexarray; index++)
    {
      if (xvariable == "nVert")
	limitup = (index + 1) * 5;
      if (xvariable == "met_sumEt")
	limitup = (index + 1) * 200;
      if (xvariable == "zll_pt")
	limitup = (index + 1) * 12;
      strlimitup = Form ("%d", limitup);
cout << variablename << endl;

cout << "index " << index << endl;
cout << "limit down " << limitdown << endl;
cout << "limit up " << limitup << endl;
      if(variablenamepng.Contains("over"))
        resolution.push_back (new TH1F (Form ("resx%d", index), " ", 200, -200, 200)); // changing the x axis for met_uPara_zlloverzll_pt, so that models are more visible
      else
        resolution.push_back (new TH1F (Form ("resx%d", index), " ", 200, -400, 400));
      


     TString totalnevents="1";
     TString lumi="40.3";
     TH1F * h0 = (TH1F*)filephys14.Get("Count");
     totalnevents=NToString(h0->Integral());
         
      TString conditionbkg="";
      

      if (xvariable == "nVert") condition = "(" + xvariable + "==" + strlimitup +")";
      else condition = "(" + xvariable + "<" + strlimitup + ")*(" + xvariable + ">" + strlimitdown + ")";
      
      conditionbkg=condition;
      if(!samplephys14.Contains("Data")) condition=condition+"*((xsec)*(puWeight)*("+lumi+"))/"+totalnevents;   
      
      cout << "condition data"  << condition.Data() << endl;
      treephys14->Draw (variablename + ">>" + TString (resolution[index]->GetName ()),			condition.Data(), "sames");

    //  c1->Print("~/www/"+TString (resolution[index]->GetName ())+".png");     
      //double m =  resolution[index]->GetMean (); (unused)
      //double um = resolution[index]->GetMean () - resolution[index]->GetRMS ();
      //double uM = resolution[index]->GetMean () + resolution[index]->GetRMS ();



      ////////
      
      RooDataHist Hist ("Hist", "Hist", x,			(TH1 *) resolution[index]->Clone ());

      RooDataHist *bkg_histogram=0;

      if(WantBKGSubtract) {
	TFile *file_ ;
	
	if (tchannel=="Gamma") file_=TFile::Open("QCD_BKG_Train.root");
	else file_=TFile::Open("TTJets13TeV.root");
	
	TString totalnbkg="1";
	TH1F * h1 = (TH1F*)file_->Get("Count");
		totalnbkg=NToString(h1->Integral());
  conditionbkg=conditionbkg+"*((xsec)*(puWeight)*("+lumi+"))/("+totalnbkg+")"; 
	TTree *treephys14bkg = (TTree *) file_->Get ("METtree");
	int bkgbin(0);
	if(variablenamepng.Contains("over"))bkgbin = 5;
	else bkgbin = 400;
	TH1F *h_ = new TH1F("h_"," ", 200, -bkgbin, bkgbin);
	cout << "condition bkg " << conditionbkg << endl;
	treephys14bkg->Draw (variablename + ">>" + TString (h_->GetName ()),	conditionbkg.Data(),"sames");
	bkg_histogram= new RooDataHist("bkg_histogram","bkg_histogram",x,h_);

      }

      // construct the voightian model
      // fit the Hist Dataset also
      // fill f and efwhm that are the parameter of the voightian
      ConstructModel(Hist, bkg_histogram, WantBKGSubtract);
                       
      //if (f/2.3 < 5) continue;
      
      RooPlot *xFrame=x.frame();
      Hist.plotOn (xFrame);

      TString titlexfit = "";
      if (variablename == "pfmetx")
	titlexfit = "MET_{x} [GeV]";
      if (variablename == "pfmety")
	titlexfit = "MET_{y} [GeV]";
      xFrame->SetXTitle (titlexfit);

      int color=kBlack;
      if (xvariable == "nVert")
	color = kRed;
      if (xvariable == "met_sumEt")
	color = kGreen+3;
      if (xvariable == "zll_pt")
	color = kBlue;
      

      cout << "plot made " << endl;
      c1->cd();
      
      //Hist.plotOn(xFrame2);
      //model->plotOn(xFrame2,RooFit::LineColor(kBlack));
      if ( WantBKGSubtract  ){
        model->plotOn(xFrame);
        model->plotOn(xFrame, Components("bkg_pdf"), LineColor(kRed), LineStyle(kDashed), FillColor(kRed), DrawOption("F"));
        model->plotOn(xFrame, Components("voigt"), LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen+1), DrawOption("L"));
      }  
      else {
	Hist.plotOn(xFrame);
        voigt->plotOn(xFrame, RooFit::FillColor(kGray), VisualizeError(*result,1), RooFit::Components(*voigt)); // 1 sigma band in gray
        voigt->plotOn(xFrame, RooFit::LineColor(color));
      }                             
      TString histoname = resolution[index]->GetName ();
      xFrame->GetYaxis()->SetRangeUser(1,200);
      xFrame->GetXaxis() ->SetRangeUser(-200,200);
      xFrame->Draw();
      c1->SetLogy();
      if (!samplephys14.Contains("jes_"))      c1->Print ("~/www/"+DestFolder+"/METModel/" + folder + "/" + tchannel +"/" + histoname + "_" +	variablenamepng + "_vs_" + xvariable + ".png");
      c1->SetLogy(0);

      //c1->Print ("~/www/"+DestFolder+"/METFits/" + folder + "/" + tchannel +"/" + histoname + "_" +	variablenamepng + "_vs_" + xvariable + ".png");

      //Print chi2/dof value

      Double_t chi2 = xFrame->chiSquare ();	//"voigt", "Hist", 3);
      //cout << "chi2 = " << chi2 << endl;

      tgraphx[index] = index;

      if (xvariable == "nVert")
	tgraphx[index] = limitup;
      if (xvariable == "met_sumEt")
	tgraphx[index] = limitup * 0.001;	//For the x axis to be in TEV
      if (xvariable == "zll_pt")
	tgraphx[index] = limitup;

      tgraphxchi2[index] = tgraphx[index];
      tgraphychi2[index] = chi2;

      if (chi2 != chi2 || chi2 >= 100)
    	tgraphychi2[index] = -0.2;
      tgraphy[index] = f / 2.3548;
      if ((variablename == "met_uPara_zll_raw/zll_pt") || (variablename =="met_uPara_zll/zll_pt")|| (variablename =="met_uPara_zll_down/zll_pt")|| (variablename =="met_uPara_zll_up/zll_pt")){
	    tgraphy[index] = -resolution[index]->GetMean ();
	    

// 	    if(tgraphy[index] > 1.0){
// 	       TString condition2="(weighttotal)*(channel=="+dileptonch +")";
//                treephys14->Draw ("met_uPara_zll >> met_uPara_zllresponse1", condition.Data());
//                treephys14->Draw ("zll_pt >> zll_ptresponse1", condition.Data(), "sames");
// 	       TString NumberStr;          // string which will contain the result
//                ostringstream convert;   // stream used for the conversion
//                convert << index;      // insert the textual representation of 'Number' in the characters in the stream
//                NumberStr = convert.str(); 
// 	       met_uPara_zllresponse1->Draw("hist");
//                met_uPara_zllresponse1->GetXaxis ()->SetTitle ("u_{|| }   [GeV]");
//                c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/met_uPara_zll_response1_"+NumberStr+".png");
//                c1->SetLogy(0);
// 	       zll_ptresponse1->Draw("hist");
//                zll_ptresponse1->GetXaxis ()->SetTitle ("zll_pt   [GeV]");
//                c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/"+tchannel +"/zll_pt_response1_"+NumberStr+".png");
//                c1->SetLogy(0);
// 	    }
	    //cout << index << "  and mean: " << -resolution[index]->GetMean () << endl;
      }
      etgraphy[index] = efwhm / 2.3548;
      if (variablename == "met_uPara_zll/zll_pt" || variablename == "met_uPara_zll_raw/zll_pt"|| variablename == "met_uPara_zll_up/zll_pt"|| variablename == "met_uPara_zll_down/zll_pt")
	    {
	    etgraphy[index] = resolution[index]->GetMeanError ();
      etgraphx[index] = 0;}


      //Set limit down
      limitdown = limitup;
      strlimitdown = Form ("%d", limitdown);


    }




  TGraph *gr =  new TGraphErrors (sizexarray, tgraphx, tgraphy, etgraphx, etgraphy);
  gr->SetMarkerColor (4);
  gr->SetMarkerStyle (21);

  TGraph *grchi2 = new TGraph (sizexarray, tgraphxchi2, tgraphychi2);
  grchi2->SetMarkerColor (2);
  grchi2->SetMarkerStyle (34);
  
  
  if (xvariable == "met_sumEt")
    {
      gr->GetXaxis ()->SetTitle ("sumE_{T} [TeV]");
      grchi2->GetXaxis ()->SetTitle ("sumE_{T} [TeV]");
    }
  if (xvariable == "nVert")
    {
      gr->GetXaxis ()->SetTitle ("Number of Vertices");
      grchi2->GetXaxis ()->SetTitle ("Number of Vertices");
    }
  if (xvariable == "zll_pt")
    {
      gr->GetXaxis ()->SetTitle ("zll_pt [GeV]");
      grchi2->GetXaxis ()->SetTitle ("zll_pt [GeV]");
    }



  gr->GetYaxis ()->SetTitle (titley);
  if (variablename == "met_uPara_zll" )
    gr->GetYaxis ()->SetTitle ("#sigma(u_{||}) [GeV]");
  if (variablename == "met_uPara_zll_raw")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{|| raw}) [GeV]");
  if (variablename == "met_uPara_zll+zll_pt")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{||} +zll_pt) [GeV]");
  if (variablename == "met_uPara_zll/zll_pt")
    gr->GetYaxis ()->SetTitle ("-<u_{||}> /zll_pt ");
  if (variablename == "met_uPara_zll_raw/zll_pt")
    gr->GetYaxis ()->SetTitle ("-<u_{|| raw}> /zll_pt ");
  if (variablename == "met_uPara_zll_raw+zll_pt")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{|| raw} +zll_pt) [GeV]");
  if (variablename == "met_uPerp_zll")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{#perp}  ) [GeV]");
  if (variablename == "met_uPerp_zll_raw")
    gr->GetYaxis ()->SetTitle ("#sigma(u_{#perp raw}  ) [GeV]");
  

  if (drawchi2) {
  grchi2->GetYaxis ()->SetTitle ("#Chi^{2}");
  grchi2->Draw ("AP");
  if (!samplephys14.Contains("jes_")) c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/" + tchannel + "/" +  variablenamepng  + "_vs_" +	     xvariable + "_chi2.png");
  c1->Clear (); }
  
  TString direction="";
  if (variablename.Contains("_up")) direction="_up_";
  if (variablename.Contains("_down")) direction="_down_";
  TFile f2 (DestFolder+folder+ direction+"tgraphs_"+samplephys14, "UPDATE");
  gr->Write (tchannel+"_"+variablenamepng + "_vs_" + xvariable);


  c1->Update ();

  TLegend *leg;
  leg = new TLegend (0.60, 0.6, 0.91, 0.81);
  leg->SetFillStyle (0);
  leg->SetBorderSize (0);
  leg->SetTextSize (0.04);
  leg->SetTextFont (42);



  leg->SetFillColor (0);





  // leg->Draw ();
  TLatex l1;
  l1.SetTextAlign (12);
  l1.SetTextSize (0.04);
  l1.SetNDC ();
  l1.DrawLatex (0.155, 0.98, "CMS Preliminary, #sqrt{s} = 13 TeV");

  
  

  gr->Draw ("AP");
  if (variablename!="met_uPara_zll/zll_pt" && variablename!="met_uPara_zll_raw/zll_pt"&& variablename!="met_uPara_zll_down/zll_pt"&& variablename!="met_uPara_zll_up/zll_pt") gr->GetYaxis()->SetRangeUser(0,70);
  else gr->GetYaxis()->SetRangeUser(0.8,1.2);
  c1->Update();
  
  
  if (variablename == "met_uPara_zll/zll_pt" || variablename=="met_uPara_zll_raw/zll_pt" || variablename=="met_uPara_zll_up/zll_pt" || variablename=="met_uPara_zll_down/zll_pt")
    {
      TLine *lineR =  new TLine ( gr->GetHistogram ()->GetXaxis ()->GetXmin (), 1, gr->GetHistogram ()->GetXaxis ()->GetXmax (), 1);
      lineR->SetLineColor (kBlue + 1);
      lineR->SetLineWidth (2);
      lineR->SetLineStyle (2);
      lineR->Draw ();
    }

  else
    {
      
      TLine *lineR =	new TLine (gr->GetHistogram ()->GetXaxis ()->GetXmin (), 0,  gr->GetHistogram ()->GetXaxis ()->GetXmax (), 0);
      lineR->SetLineColor (kBlue + 1);
      lineR->SetLineWidth (2);
      lineR->SetLineStyle (2);
      lineR->Draw ();


    }

  

   if (!samplephys14.Contains("jes_"))  c1->Print ("~/www/"+DestFolder+"/METResolution/" + folder + "/" + tchannel + "/" + variablenamepng  + "_vs_" +	     xvariable + ".png");

  */

}

//----------------------------------------------------------                                                                         
bool anaMETPerformance::CheckPid(particleBase *p) {
  //check if generated particle is muon                                                                                              
  genParticle *gp = dynamic_cast<genParticle*>(p);
  if(!gp) return kFALSE;
  if(abs(gp->GetPID())==13) return kTRUE;
  return kFALSE;
}

std::pair<double, double> 
anaMETPerformance::compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  double uX = -metPx;
  double uY = -metPy;
  uX -= qX;
  uY -= qY;
  
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

std::pair<double, double> 
anaMETPerformance::compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count)
{
  if ( zP4->Pt() == 0. ) {
    Warning ("compMEtProjU", " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!");
    errorFlag = 1;
    return std::pair<double, double>(0., 0.);
  }
  
  double qX = zP4->Px();
  double qY = zP4->Py();
  double qT = TMath::Sqrt(qX*qX + qY*qY);
  
  TLorentzVector zBoson = zP4->GetLorentzVector();
  TLorentzVector hadronicRecoil = -(MET+zBoson);

  double uX = hadronicRecoil.Px();
  double uY = hadronicRecoil.Py();
  uX -= qX;
  uY -= qY;
  
  double u1 = (uX*qX + uY*qY)/qT;
  double u2 = (uX*qY - uY*qX)/qT;

  return std::pair<double, double>(u1,u2);
}

//----------------------------------------------------------
void 
anaMETPerformance::CreateOutputObjects() {

  anaBaseTask::CreateOutputObjects();

  if(!fOutput) {
    Printf("anaMETPerformance: fOutput not present");
    return;
  }

  fh1NMuons = new TH1F("fh1NMuons","fh1NMuons;#it{N}_{muons}",101,-0.5,100.5);
  fOutput->Add(fh1NMuons);

  fh1uParaZllPt =  new TH1F("fh1uParaZllPt","ffh1uParaZllPt;upara+qt",100,0.,200.);
  fOutput->Add(fh1uParaZllPt);

  fh2MetCent = new TH2F("fh2MetCent","fh2MetCent;centrality;MET",100,0,100,500,0,1000.);
  fOutput->Add(fh2MetCent);

  fh2SumEtCent = new TH2F("fh2SumEtCent","fh2SumEtCent;centrality;sumET",100,0,100,500,0,10000.);
  fOutput->Add(fh2SumEtCent);

  fh3PtEtaPhi = new TH3F("fh3PtEtaPhi","fh3PtEtaPhi;pt;eta;phi",500,0,500,100,-5,5,72,-TMath::Pi(),TMath::Pi());
  fOutput->Add(fh3PtEtaPhi);

  for(Int_t j = 0; j<10; ++j) {
    fh2MetCentPtMin[j] = new TH2F(Form("fh2MetCentPtMin%d",j),"fh2MetCent;centrality;MET",100,0,100,500,0,1000.);
    fOutput->Add(fh2MetCentPtMin[j]);
  }
    
  
}



double
anaMETPerformance::FWHM (double sigma, double gamma)
{

  double f_g = 2 * sigma * sqrt (2 * log (2));
  double f_l = 2 * gamma;

  return 0.5346 * 2 * gamma + sqrt (0.2166 * f_l * f_l + f_g * f_g);
}



double
anaMETPerformance::FWHMError (double sigma, double gamma, double esigma, double egamma,
	   double Vss, double Vsg, double Vgs, double Vgg)
{


  double a = 0.5346;
  double b = 0.2166;
  double ef_g = 2 * esigma * sqrt (2 * log (2));
  double ef_l = 2 * egamma;

  double dg =
    2 * a + 4 * b * gamma / sqrt (4 * b * pow (gamma, 2) +
				  4 * pow (sigma, 2) * log (2));

  double ds =
    (sigma * log (4)) / sqrt (b * pow (gamma, 2) + pow (sigma, 2) * log (2));
  
  double p1 = ef_l * ef_l * Vgg * dg;
  double p2 = ef_g * ef_l * Vsg * dg * ds;	//identical (should be)
  double p3 = ef_g * ef_l * Vgs * dg * ds;
  double p4 = ef_g * ef_g * Vss * ds;

  return sqrt (abs (p1) + abs (p2) + abs (p3) + abs (p4));

}


double
anaMETPerformance::FWHMError_fixed (double sigma, double gamma, double esigma, double egamma,
	   double Vss, double Vsg, double Vgs, double Vgg)
{
  // Vss = correlation(sigma, sigma)
  // Vsg = correlation(sigma, gamma)
  // etc
  double a = 0.5346;
  double b = 0.2166;
  double c = 2 * sqrt( 2*log(2) );
  double f_g = c * sigma;
  double f_l = 2 * gamma;
  double sq = sqrt( b * pow(f_l, 2) + pow(f_g, 2) );
  
  // Partial derivatives of f_voigtian w.r.t sigma and gamma
  // f = a * f_l + sqrt( b * f_l^2 + f_g^2 )
  double dfds = c * ( f_g / sq ) ;
  double dfdg = 2 * ( a + b * f_l / sq ) ;
  
  // esigma * esigma * pow( Vss, 2 ) gives covariance(sigma, sigma) etc
  double p1 = dfds * dfds * esigma * esigma * pow( Vss, 2 );
  double p2 = dfds * dfdg * esigma * egamma * pow( Vsg, 2 );
  double p3 = dfdg * dfds * egamma * esigma * pow( Vgs, 2 );
  double p4 = dfdg * dfdg * egamma * egamma * pow( Vgg, 2 );

  return sqrt ( p1 + p2 + p3 + p4 );
}
TString
anaMETPerformance::NToString(Float_t type) 
{ 
   TString name; name.Form("%f",type); 
      cout << name << endl;   
         return name;
         
            } 
            
            
