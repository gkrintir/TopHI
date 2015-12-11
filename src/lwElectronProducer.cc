//
// producer of muon candidates
//

#include "UserCode/diall/interface/lwElectronProducer.h"
#include "UserCode/diall/interface/genParticle.h"
#include "UserCode/diall/interface/lwElectron.h"


//ClassImp(lwElectronProducer)

//__________________________________________________________
lwElectronProducer::lwElectronProducer() :
inputBase("lwElectronProducer"),
  flwElectronsRecoName("lwElectronsReco"),
  flwElectronsReco(0x0),
  flwElectronsGeneName("lwElectronsGene"),
  flwElectronsGene(0x0),
  fElectrons(),
  fPtMin(10.),
  fMaxEtaAbs(2.1),
  fMaxTrkChi2(4.),
  fMaxGlbChi2(10.),
  fMinNMuHits(0),
  fMinMS(1),
  fMaxDxy(0.2),//3.),
  fMaxDz(0.5),//15.),
  fMaxtrkDxy(0.3),//3.),
  fMaxtrkDz(20.),
  fMinNPixHits(0),
  fMinTrkLWM(5)
{
  //default constructor
}

//__________________________________________________________
lwElectronProducer::lwElectronProducer(const char *name) :
  inputBase(name),
  flwElectronsRecoName("lwElectronsReco"),
  flwElectronsReco(0x0),
  flwElectronsGeneName("lwElectronsGene"),
  flwElectronsGene(0x0),
  fElectrons(),
  fPtMin(16.),
  fMaxEtaAbs(2.1),
  fMaxTrkChi2(4.),
  fMaxGlbChi2(10.),
  fMinNMuHits(0),
  fMinMS(1),
  fMaxDxy(0.2),//3.),
  fMaxDz(0.5),//15.),
  fMaxtrkDxy(0.3),//3.),
  fMaxtrkDz(20.),
  fMinNPixHits(0),
  fMinTrkLWM(5)
{
  //standard constructor
}

//__________________________________________________________
Bool_t lwElectronProducer::Init() {

  if(!inputBase::Init()) return kFALSE;
  
  if(fInputMode==hiForest) {
    // Gen Info
    if (fChain->GetBranch("Gen_nptl"))
      fChain->SetBranchAddress("Gen_nptl", &fElectrons.Gen_nptl, &fElectrons.b_Gen_nptl);
    if (fChain->GetBranch("Gen_pid"))
      fChain->SetBranchAddress("Gen_pid", &fElectrons.Gen_pid, &fElectrons.b_Gen_pid);
    if (fChain->GetBranch("Gen_mom"))
      fChain->SetBranchAddress("Gen_mom", &fElectrons.Gen_mom, &fElectrons.b_Gen_mom);
    if (fChain->GetBranch("Gen_pt"))
      fChain->SetBranchAddress("Gen_pt", &fElectrons.Gen_pt, &fElectrons.b_Gen_pt);
    if (fChain->GetBranch("Gen_eta"))
      fChain->SetBranchAddress("Gen_eta", &fElectrons.Gen_eta, &fElectrons.b_Gen_eta);
    if (fChain->GetBranch("Gen_phi"))
      fChain->SetBranchAddress("Gen_phi", &fElectrons.Gen_phi, &fElectrons.b_Gen_phi);
    // Reco Info
    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("nMu", 1);
    fChain->SetBranchStatus("mu*", 1);
    if (fChain->GetBranch("nMu"))
      fChain->SetBranchAddress("nMu", &fElectrons.Glb_nptl, &fElectrons.b_Glb_nptl);
    if (fChain->GetBranch("muCharge"))
      fChain->SetBranchAddress("muCharge", &fElectrons.Glb_charge, &fElectrons.b_Glb_charge);
    if (fChain->GetBranch("muPt"))
      fChain->SetBranchAddress("muPt", &fElectrons.Glb_pt, &fElectrons.b_Glb_pt);
    if (fChain->GetBranch("muEta"))
      fChain->SetBranchAddress("muEta", &fElectrons.Glb_eta, &fElectrons.b_Glb_eta);
    if (fChain->GetBranch("muPhi"))
      fChain->SetBranchAddress("muPhi", &fElectrons.Glb_phi, &fElectrons.b_Glb_phi);
    if (fChain->GetBranch("muD0"))
      fChain->SetBranchAddress("muD0", &fElectrons.Glb_dxy, &fElectrons.b_Glb_dxy);
    if (fChain->GetBranch("muDz"))
      fChain->SetBranchAddress("muDz", &fElectrons.Glb_dz, &fElectrons.b_Glb_dz);
    if (fChain->GetBranch("muPixelHits"))
      fChain->SetBranchAddress("muPixelHits", &fElectrons.Glb_nValPixHits, &fElectrons.b_Glb_nValPixHits);
    //if (fChain->GetBranch("Glb_nValTrkHits"))
    //fChain->SetBranchAddress("Glb_nValTrkHits", &fElectrons.Glb_nValTrkHits, &fElectrons.b_Glb_nValTrkHits);
    if (fChain->GetBranch("muElectronHits"))
      fChain->SetBranchAddress("muElectronHits", &fElectrons.Glb_nValMuHits, &fElectrons.b_Glb_nValMuHits);
    if (fChain->GetBranch("muTrkQuality"))
      fChain->SetBranchAddress("muTrkQuality", &fElectrons.Glb_trkQuality, &fElectrons.b_Glb_trkQuality);
    if (fChain->GetBranch("muIsGood"))
      fChain->SetBranchAddress("muIsGood", &fElectrons.Glb_isGood, &fElectrons.b_Glb_isGood);
    if (fChain->GetBranch("muChi2NDF"))
      fChain->SetBranchAddress("muChi2NDF", &fElectrons.Glb_glbChi2_ndof, &fElectrons.b_Glb_glbChi2_ndof);
    if (fChain->GetBranch("muStations"))
      fChain->SetBranchAddress("muStations", &fElectrons.Glb_nMatchedStations, &fElectrons.b_Glb_nMatchedStations);
    if (fChain->GetBranch("muInnerD0"))
      fChain->SetBranchAddress("muInnerD0", &fElectrons.Glb_trkDxy, &fElectrons.b_Glb_trkDxy);
    if (fChain->GetBranch("muInnerDz"))
      fChain->SetBranchAddress("muInnerDz", &fElectrons.Glb_trkDz, &fElectrons.b_Glb_trkDz);
    if (fChain->GetBranch("muPixelLayers"))
      fChain->SetBranchAddress("muPixelLayers", &fElectrons.Glb_pixLayerWMeas, &fElectrons.b_Glb_pixLayerWMeas);
    if (fChain->GetBranch("muTrkLayers"))
      fChain->SetBranchAddress("muTrkLayers", &fElectrons.Glb_trkLayerWMeas, &fElectrons.b_Glb_trkLayerWMeas);
    if (fChain->GetBranch("muPFChIso"))
      fChain->SetBranchAddress("muPFChIso", &fElectrons.Glb_pfChIso, &fElectrons.b_Glb_pfChIso);
    if (fChain->GetBranch("muPFPhoIso"))
      fChain->SetBranchAddress("muPFPhoIso", &fElectrons.Glb_pfPhoIso, &fElectrons.b_Glb_pfPhoIso);
    if (fChain->GetBranch("muPFNeuIso"))
      fChain->SetBranchAddress("muPFNeuIso", &fElectrons.Glb_pfNeuIso, &fElectrons.b_Glb_pfNeuIso);
    if (fChain->GetBranch("muPFPUIso"))
      fChain->SetBranchAddress("muPFPUIso", &fElectrons.Glb_pfPUIso, &fElectrons.b_Glb_pfPUIso);
    
    fInit = kTRUE;
  }
  return kTRUE;
}

//__________________________________________________________
Bool_t lwElectronProducer::InitEventObjects() {

  //Create event objects
  if(!fEventObjects) {
    Printf("%s: fEventObjects does not exist. Cannot store output",GetName());
    return kFALSE;
  } else {
    if(!fEventObjects->FindObject(flwElectronsRecoName)) {
      flwElectronsReco = new TClonesArray("lwElectron");
      flwElectronsReco->SetName(flwElectronsRecoName);
      fEventObjects->Add(flwElectronsReco);
    }
    if(!fEventObjects->FindObject(flwElectronsGeneName) && !flwElectronsGeneName.IsNull()) {
      flwElectronsGene = new TClonesArray("genParticle");
      flwElectronsGene->SetName(flwElectronsGeneName);
      fEventObjects->Add(flwElectronsGene);
    }
  }
  
  return kTRUE;
}

//__________________________________________________________
Bool_t lwElectronProducer::Run(Long64_t entry) {

  //overloaded run funtion
  Long64_t centry = LoadTree(entry);
  if(centry<0) return kFALSE;

  if(!InitEventObjects()) return kFALSE;
  
  //clear arrays
  flwElectronsReco->Delete();
  if(flwElectronsGene) flwElectronsGene->Delete();

  //reconstructed muons
  Int_t muCount = 0;
  for(Int_t i = 0; i<fElectrons.Glb_nptl; i++) {
    if(!AcceptElectron(i)) continue;
    lwElectron *mu = new lwElectron(fElectrons.Glb_pt->at(i),
                            fElectrons.Glb_eta->at(i),
                            fElectrons.Glb_phi->at(i),
                            0,
                            i);
    mu->SetCharge(fElectrons.Glb_charge->at(i));
    (*flwElectronsReco)[muCount] = mu;
    ++muCount;
  }
  flwElectronsReco->Sort();
  //Printf("%d reconstructed muons",muCount);

  //generated muons
  if(flwElectronsGene) {
    muCount = 0;
    for(Int_t i = 0; i<fElectrons.Gen_nptl; i++) {
      genParticle *mu = new genParticle(fElectrons.Gen_pt[i],
                                        fElectrons.Gen_eta[i],
                                        fElectrons.Gen_phi[i],
                                        0,
                                        i);
      mu->SetCharge(fElectrons.Gen_pid[i]/abs(fElectrons.Gen_pid[i]));
      mu->SetPID(fElectrons.Gen_pid[i]);
      mu->SetPIDMom(fElectrons.Gen_mom[i]);
      (*flwElectronsGene)[muCount] = mu;
      ++muCount;
    }
    flwElectronsGene->Sort();
    //Printf("%d generated muons",muCount);
  }
  
  return kTRUE;
}

//__________________________________________________________
Bool_t lwElectronProducer::AcceptElectron(Int_t i) {

  //muon quality selection 
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  //plus https://github.com/CmsHI/quickZMacros/blob/master/ggHistos.C#L4
  if(!(fElectrons.Glb_isGood->at(i)))                       return kFALSE;
  else if((fElectrons.Glb_pt->at(i))<fPtMin)                return kFALSE;
  else if((fabs(fElectrons.Glb_eta->at(i)))>fMaxEtaAbs)     return kFALSE;
  //else if((fElectrons.Glb_trkChi2_ndof->at(i))>fMaxTrkChi2) return kFALSE;
  else if((fElectrons.Glb_glbChi2_ndof->at(i))>fMaxGlbChi2) return kFALSE;
  else if((fElectrons.Glb_nValMuHits->at(i))<fMinNMuHits)   return kFALSE;
  else if((fElectrons.Glb_nMatchedStations->at(i))<fMinMS)  return kFALSE;
  //else if((fElectrons.Glb_dxy->at(i))>fMaxDxy)            return kFALSE;
  //else if((fElectrons.Glb_dz->at(i))>fMaxDz)              return kFALSE;
  else if((fabs(fElectrons.Glb_trkDxy->at(i)))>fMaxDxy)     return kFALSE;
  else if((fabs(fElectrons.Glb_trkDz->at(i)))>fMaxDz)       return kFALSE;
  else if((fElectrons.Glb_nValPixHits->at(i))<fMinNPixHits) return kFALSE;
  else if((fElectrons.Glb_trkLayerWMeas->at(i))<fMinTrkLWM) return kFALSE;
  else if (!(fElectrons.Glb_trkQuality->at(i)))             return kFALSE; 
  else return kTRUE;
}

//__________________________________________________________
Long64_t lwElectronProducer::LoadTree(Long64_t entry) {

  //overloaded LoadTree function 
  if(!fChain) {
    Printf("fChain doesn't exist");
    return -5;
  }
  
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Init();
    Printf("%lld fCurrent: %d",entry,fCurrent);
  }

  //  fChain->SetMakeClass(1);
 
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) {
    Printf("hiEventProducer: centry smaller than 0");
    return centry;  
  }
  
  fChain->GetEntry(entry);

  return centry;
}
