#ifndef lwElectronProducer_h
#define lwElectronProducer_h

//
// muon candidate producer
//

#include <TNamed.h>
#include <TChain.h>
#include <TString.h>
#include <TClonesArray.h>

#include "UserCode/diall/interface/inputBase.h"
#include <UserCode/diall/interface/ForestElectrons.h>

class lwElectronProducer : public inputBase {
 public:
  lwElectronProducer();
  lwElectronProducer(const char *name);
  virtual ~lwElectronProducer() {;}

  Bool_t   Init();
  Long64_t LoadTree(Long64_t entry);
  Bool_t   InitEventObjects();
  Bool_t   Run(Long64_t entry);
  
  void     SetlwElectronsRecoName(TString n) { flwElectronsRecoName = n; }
  void     SetlwElectronsGeneName(TString n) { flwElectronsGeneName = n; }

  const char* GetlwElectronsRecoName() const { return flwElectronsRecoName.Data() ; }
  const char* GetlwElectronsGeneName() const { return flwElectronsGeneName.Data() ; }
  
 protected:
  Bool_t   AcceptElectron(Int_t i);
  
  TString                      flwElectronsRecoName;// name of reco muons
  TClonesArray                *flwElectronsReco;    //!reco muons
  TString                      flwElectronsGeneName;// name of gene muons
  TClonesArray                *flwElectronsGene;    //!gene muons
  ForestElectrons                  fElectrons;          //! Electrons in forest tree
  Float_t                      fPtMin;          // minimum pT
  Float_t                      fMaxEtaAbs;      // max eta
  Float_t                      fMaxTrkChi2;     // max chi2
  Float_t                      fMaxGlbChi2;     // max chi2
  Int_t                        fMinNMuHits;     // min muon hits
  Int_t                        fMinMS;          // #matched stations
  Float_t                      fMaxDxy;         // max dxy
  Float_t                      fMaxDz;          // max dz
  Float_t                      fMaxtrkDxy;      // max innerTrack dxy
  Float_t                      fMaxtrkDz;       // max innerTrack dz
  Int_t                        fMinNPixHits;    // min pixel hits
  Int_t                        fMinTrkLWM;      // min tracker layer hits

 private:
  lwElectronProducer(const lwElectronProducer& obj); // copy constructor
  lwElectronProducer& operator=(const lwElectronProducer& other); // assignment
  
  //ClassDef(lwElectronProducer,2)
};
#endif
