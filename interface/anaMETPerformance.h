#ifndef anaMETPerformance_h
#define anaMETPerformance_h

#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TClonesArray.h"

#include "UserCode/diall/interface/anaBaseTask.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/diParticle.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooVoigtian.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <utility>
//
// geometrical matching of muons to other pfParticleBase object in fMatch
// match info is stored in muons (fMuons)
//

class anaMETPerformance : public anaBaseTask {
   
public:
  enum metType {
    kGen   = 0,
    kGenEm = 1,
    kPFRaw = 2,
    kVS    = 3,
    kPuppi = 4,
    kCS    = 5
  };
  
   anaMETPerformance() {;}
   anaMETPerformance(const char *name, const char *title);
   virtual ~anaMETPerformance() {;}
   void Exec(Option_t *option="");
   void ConstructModel(RooDataHist Hist, RooDataHist *bkg_hist, bool BKGSubtract);
   std::pair<double, double> compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count);
   std::pair<double, double> compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count);
   void CreateOutputObjects();
   double FWHM (double sigma, double gamma);
   double FWHMError (double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg);
   double FWHMError_fixed (double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg);
   TString NToString(Float_t type);
   void SetCheckPid(Bool_t b)          { fCheckPid = b; }
   void SetMetType(metType t)          { fMetType = t; }
   void SetMinPt(Float_t m)            { fMinPt = m; }
   void SetMuonsName(TString name)     { fMuonsName = name ; }
   void SetParticlesName(TString name) { fParticlesName = name ; }
   
 protected:
   Bool_t            CheckPid(particleBase *p);

   Bool_t            fCheckPid;             //check if candidates are really muons (for simulation)
   metType           fMetType;              //matching type (defines where to store)
   Float_t           fMinPt;                //minimum pT of particles
   TString           fMuonsName;            //name of particles
   TClonesArray     *fMuons;                //!muon array
   TClonesArray     *fParticles;            //!muon array
   TString           fParticlesName;        //name of particles
   TClonesArray     *fZs;                   //!Z candidates container
   TString           fZsName;               // name of Z candidates

   TH1F             *fh1NMuons;             //!# selected muons in event
   TH1F             *fh1uParaZllPt;         //!# Upara (Hadronic recoil paralllel to qT) plus qT (trans. momemtum of dilepton)
   TH2F             *fh2MetCent;            //!MET vs centrality
   TH2F             *fh2SumEtCent;          //!SumEt vs centrality
   TH3F             *fh3PtEtaPhi;           //!particle pt vs eta vs phi
   TH2F             *fh2MetCentPtMin[10];   //!MET vs centrality for various min pt cuts
   
   //ClassDef(anaMETPerformance,2)
};
#endif
