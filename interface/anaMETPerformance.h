#ifndef anaMETPerformance_h
#define anaMETPerformance_h

#include "TString.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "UserCode/diall/interface/anaBaseTask.h"
#include "UserCode/diall/interface/particleBase.h"
#include "UserCode/diall/interface/diParticle.h"

#include <utility>

const Int_t MAXMEC = 300;

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
   std::pair<double, double> compMETProjU(diParticle* zP4, double metPx, double metPy, int& errorFlag, int count);
   std::pair<double, double> compHadronicRecoilProjU(diParticle* zP4, TLorentzVector MET, int& errorFlag, int count);
   void CreateOutputObjects();

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
   TString           fParticlesName;        //name of particles
   TClonesArray     *fMuons;                //!muon array
   TClonesArray     *fParticles;            //!muon array
   TClonesArray     *fZs;                   //!Z candidates container
   TString           fZsName;               // name of Z candidates

   TTree            *fMETPerformanceInfo;
   Int_t             fNZCands;
   Float_t           fMassZCands;
   Float_t           fPtZCands;
   Float_t           fEtaZCands;
   Float_t           fMET;
   Float_t           fuParaZll;
   Float_t           fuParaZllPt;
   Float_t           fuPerpZll;
   //TH2F             *fh2MetCentPtMin[10];   //!MET vs centrality for various min pt cuts
   
   ClassDef(anaMETPerformance,1)
};
#endif
